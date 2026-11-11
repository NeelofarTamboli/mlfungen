# McSplicer notes

Checked 2026-11-11 against `canzarlab/McSplicer` commit `d6a22c3` and
`canzarlab/sigcount` commit `c738d3f` (both cloned into the git-ignored
`external/`). Two kinds of evidence: reading the source, and running things.
Runs: (1) the bundled single-end ENOPH1 example, (2) a small real K562
paired-end subset described in "Validation run" below. Items marked **not
tested** were not run.

Environment: Python 3.13, numpy 2.5, the bundled Linux x86_64 binaries (no
build needed).

## Decisions for this project (2026-11-11)

**1. The first regression target is internal acceptor usage only.** Dataset rows
are McSplicer `s` rows that the transcript-boundary rule (below, "Which `s`/`e`
sites are splice sites") classifies as `acceptor`, K562 only. Transcript starts,
transcript ends and "neither" rows are dropped. Donor (`e`) rows are still
produced and kept in the QC/reporting outputs (site tables, support counts,
usage histograms) so the donor finding is documented, but they are not in the
first regression dataset.

Why donors are excluded: `McSplicer.py` hardcodes `disable_end_sites_computation
= True`, so donor values are not EM estimates but spliced / (spliced + contiguous)
read counts, which is 1 for any ordinary donor. In the real K562 test region
(rep 1, 25 internal donors) 20 are exactly `0.999999`, 4 are `1e-06` (3 with no
reads at all, 1 with 413 contiguous and 0 spliced reads) and 1 is intermediate
(0.375). The 26 internal acceptors in the same region cover the range (10 at
>= 0.99, 10 between 0.05 and 0.95, 3 below 0.05). Turning the donor EM path on
crashes (`ValueError` in `cases_code_new.builder_end_10`), so there is no
working alternative in this release. Details in "How the numbers are computed".

**2. K562 reads are counted strand-aware.** The library is reverse-stranded and
sigcount ignores strand, so each BAM is split by inferred transcript strand
before counting. For each replicate:
1. Split the BAM by fragment strand with `fragment_strand` (read 1 reverse or
   read 2 forward = `+`; otherwise `-`): `src/validation/split_bam_by_strand.py`.
2. Split the original GENCODE GTF by strand (column 7).
3. Per strand: `exonRefine` on that strand's GTF, `bin/sigcount` (no `-s`) on
   that strand's BAM, `McSplicer.py --read_len 100`.
4. Concatenate the two sets of output files (a `gene_id` never spans strands).
Tested on the 60 kb subsets only (see "Strand" and "Validation run"); the full
BAMs have not been processed. Other settings: `bin/sigcount`, `--read_len 100`.

## Pipeline

```
annotation.gtf --exonRefine--> refined.gtf
BAM + refined.gtf --sigcount--> prefix.cnt
refined.gtf + prefix.cnt --McSplicer.py--> <out_dir>/<chr>/<gene_id>[_<prefix>].csv
```

```bash
exonRefine annotation.gtf -p refined
sigcount aligned.bam refined.gtf prefix            # writes prefix.cnt
python3 McSplicer.py --gtf refined.gtf --count_file prefix.cnt \
    --out_dir out --read_len 100 --bootstraps 10 [--gene_id G] [--prefix P]
```

`--gtf --count_file --out_dir --read_len` are required. `--anno y` **not tested**.
BAM must be sorted and indexed per the README; the example BAMs are. I did
not test an unsorted/unindexed BAM.

## Inputs

- **Annotation GTF**: `exon` lines with `gene_id` and `transcript_id`.
  `exonRefine` accepted real GENCODE v29 lines (gene/transcript/CDS/UTR too),
  kept versioned IDs (`ENSG00000116459.10`), and writes `exon` + `subexon` lines.
  `McSplicer.py` reads only `subexon` lines. `exonRefine` splits exons at the
  boundaries of **all** genes in the file, including other overlapping genes on
  the same strand (see "Overlapping genes").
- **`.cnt` file**: one line per read/fragment path,
  `nodes<TAB>totalCount<TAB>uniqueCount`; paired fragments are `left^right`.
  `McSplicer.py` uses only `totalCount`, which includes fractional weights for
  multi-mapped reads (e.g. `920.75` vs unique `919`). Header lines (`#ReadLen`,
  `#FLen`, ...) are skipped by `McSplicer.py`.

## sigcount builds and paired-end data

Tested on the real paired-end subset (both replicates): the plain `sigcount`
and `sigcount_pe` binaries wrote **byte-identical** `.cnt` files, with `left^right`
lines for every fragment, `#ReadLen: 100`, `#FLen: 242` (rep 1) / 257 (rep 2).
On the single-end example, `sigcount` works and `sigcount_pe` writes only the
header. So in these tests `bin/sigcount` handled both single-end and paired data and
`sigcount_pe` only pairs. Use `bin/sigcount`.

This differs from the sigcount GitHub source, where the read mode defaults to
`SINGLE` and the `--mode` option is commented out. The shipped binaries do not
behave like that source. Trust the observed behaviour, and re-test if the
binaries change.

- Fragments whose mate is not in the BAM are skipped with a
  `WARNING: Mate could not be found` (6 in each subset: multi-mapped fragments
  whose other mate lies outside the extracted region; a full BAM should not have
  these).
- `McSplicer.py` splits `left^right` into two independent subpaths and adds the
  fragment count to each (`parse_cnt_file`). The mate pairing, `#FLen`, and the
  PE-aware modules (`compute_frag_length.py`, `populate_pe_ws.py`) are never
  used; `McSplicer.py` does not import them and the `paired_end` argument is never
  passed. So each mate is treated as an independent read (by the code; an
  overlapping pair then counts as two reads over the shared bases, not tested).

## `--read_len`

It is the length of **one mate** (100 for ENCODE 2x100), not the fragment (200).
- Source: `read_length` is used only in `preprocess_read_counts`
  (`subexon_length - (read_length - 1)` = number of positions where a contiguous
  read of that length fits inside a segment) and `get_length_list` (same count per
  segment). Both describe a single contiguous read, and the mates are processed
  as separate reads (above).
- sigcount's own header for the paired-end data says `#ReadLen: 100`.
- Controlled test, rep 1 subset, `--read_len` 50 / 100 / 200 (internal acceptors,
  n = 26): |Δ usage| vs 100 has mean 0.008 (50) and 0.019 (200), max 0.14 and
  0.32. Donor values did not change at all (they do not use it, see below).
  Transcript-start `s` sites changed by up to 0.07 / 0.11. So a wrong value
  mostly moves acceptor estimates on low-coverage or alternative sites.
- STAR QC's "average input read length 200" is the two mates summed. It is not
  the value to use.

## Strand

- The K562 libraries are reverse-stranded and the BAM FLAGs agree: for a
  fragment, read 1 reverse **or** read 2 forward means a `+` strand transcript;
  read 1 forward or read 2 reverse means `-`. Measured on uniquely mapped
  alignments fully inside exons that overlap no gene on the other strand: ATP5PB
  (`+`) 13,345 / 17,927 alignments (rep 1 / rep 2) all `+`; WDR77 (`-`) 5,961 /
  9,969 all `-`. In the 350 bp where their exons overlap, 99.3% (rep 1) and
  100% (rep 2) of fragments are `-`.
- sigcount does **not** use FLAGs or library orientation. `-s` only reads
  TopHat-style `XS:A` tags; the ENCODE BAMs have no `XS` (STAR was run without
  `--outSAMstrandField`), so `-s` has nothing to read. Without `-s` every
  fragment is tested against both strands' subexons. `McSplicer.py` never looks
  at strand in the reads.
- Effect, measured by comparing the normal run with a strand-split run (BAM split
  by the rule above, each half run against only its own strand's genes, outputs
  merged): total count 22,290 vs 22,268 (22 count units counted on both
  strands); on the 84 sites present in both, mean |Δ usage| 0.0004, max 0.023 (rep
  1) / 0.014 (rep 2). All differences are at sites in the ATP5PB/WDR77 antisense
  overlap zone. One site exists only in the strand-unaware run: the `+` pseudogene
  UBE2FP3 start site gets usage 0.998 from reads whose FLAGs say `-`; it has no counts at all in the
  strand-split run.
- So reverse-stranded paired data are not handled by the tools; they are
  tolerated when genes do not overlap antisense. The project uses the strand
  split described under "Decisions" (step 1 can also be done with `samtools
  view`: `+` = `-f 0x80 -F 0x10` plus `-f 0x50`, `-` = `-f 0x40 -F 0x10` plus
  `-f 0x90`; not run here, samtools is not installed).

## Output

One CSV per gene: `<out_dir>/<chr>/<gene_id>[_<prefix>].csv`, header
`step,index,strand,chr,splice site locus,usage estimate`. `chr` is the GTF
seqname (`chr1` with the UCSC-named GENCODE file). No read-count column.
Bootstrap: step 0 is the estimate from the observed counts; steps 1..N repeat
the estimate on counts resampled multinomially (same total, `np.random.seed(42)`
at import). Only the replicate values are written, no SD; compute it from steps
1..N. Rep 1 subset, 10 bootstraps, internal acceptors: SD 0.000-0.035 for the
two well covered genes, but up to 0.17 in C1orf162 (5-6 spliced reads per site).

Real step-0 rows from the bundled single-end example (control sample):

```
0,s3,+,4,83369073,0.9577765854586201
0,s4,+,4,83372196,0.9524496589503981
0,e0,+,4,83352037,0.08170033162496232
0,e1,+,4,83352066,0.999999
0,e6,+,4,83382208,1e-06
```

and from the real K562 paired-end subset (rep 1, `read_len 100`; abridged,
values rounded, comments added):

```
step,index,strand,chr,splice site locus,usage estimate
0,s7,+,chr1,111449837,0.9792        # ATP5PB acceptor
0,s13,+,chr1,111460917,0.1490       # ATP5PB acceptor, alternative
0,e1,+,chr1,111449581,0.999999      # ATP5PB donor
0,s6,-,chr1,111447144,0.9998        # WDR77 acceptor
0,e6,-,chr1,111443349,1e-06         # WDR77 transcript end, not a splice site
```

## Coordinates (verified against the GRCh38 sequence)

`splice site locus` L is the GTF coordinate of the subexon boundary, 1-based
inclusive; McSplicer does not convert. It is always an **exonic** base:
- `s` (start site): first exonic base in transcription direction. `+`: the
  smaller coordinate of the exon; `-`: the larger.
- `e` (end site): last exonic base in transcription direction. `+`: larger
  coordinate; `-`: smaller.

Using the sequence (`seq` = chromosome as a 0-based Python string, so 1-based
position p is `seq[p-1]`), `j` = 0-based index of the first base to the right
(higher coordinate) of the exon/intron junction:

| strand | McSplicer site | junction type | j | two intron bases next to the junction |
|---|---|---|---|---|
| `+` | `s` | acceptor | `L-1` | `seq[L-3:L-1]` should be `AG` |
| `+` | `e` | donor | `L` | `seq[L:L+2]` should be `GT` |
| `-` | `s` | acceptor | `L` | `revcomp(seq[L:L+2])` should be `AG` |
| `-` | `e` | donor | `L-1` | `revcomp(seq[L-3:L-1])` should be `GT` |

Window of W bases centred on the junction: `w = seq[j-W//2 : j+W//2]`,
reverse-complemented if the strand is `-`. After that, `w[:W//2]` is upstream
of the junction and `w[W//2:]` downstream, on both strands. Then for a donor
`w[W//2:W//2+2]` is the intron start (GT) and `w[W//2-1]` is the last exonic
base; for an acceptor `w[W//2-2:W//2]` is the intron end (AG) and `w[W//2]` is
the first exonic base.

Checked on every annotated internal junction of chr1 in GENCODE v29 (28,669
donors, 28,109 acceptors, both strands) against the GRCh38 FASTA: at these
loci donors read GT in 96.4% (`+`) / 96.1% (`-`) and acceptors read AG in
98.1% / 98.2% (most non-GT donors read GC or AT). Moving the locus by ±1 or ±2
drops GT to 0.1-1.6% and AG to 0.1-0.6%, except acceptors shifted by +2 on the
`+` strand and -2 on the `-` strand, which still read AG 6% of the time (not
investigated). The
window formula gives the same GT/AG test result for W = 80, 400, 2000, 10000
(97.2% of 56,778 junctions; none run off the chromosome end on chr1). Also
checked on the McSplicer outputs themselves (below): all 26 internal acceptors
read AG and all 25 internal donors read GT, on both strands.

## Which `s`/`e` sites are splice sites

`s`/`e` mark every subexon boundary labelled L/R/B, so they include transcript
starts and ends and boundaries from other overlapping genes. Rule, using the
original GENCODE transcripts (not the refined GTF), per McSplicer output file
(= per `gene_id`): order each transcript's exons in transcription direction
(descending coordinate on `-`). Then
- acceptor set = start (in transcription direction) of every exon except the
  first; donor set = end of every exon except the last;
- TSS set = start of the first exon; TES set = end of the last exon.

Keep an `s` row if its locus is in the gene's acceptor set, and an `e` row if
its locus is in the gene's donor set. For the first regression dataset only the
kept `s` rows are used (see "Decisions"); kept `e` rows go to QC only. Drop rows that are only in TSS/TES, and
rows in neither set. A locus that is internal in one transcript and a TSS/TES in
another is kept (none occurred in the test region, 0 of 85).

Test region result (rep 1, 85 rows from 6 genes): 26 acceptors + 25 donors kept
(all canonical AG/GT); 17 TSS, 11 TES and 6 "neither" dropped (only 8 of those
34 happen to read GT/AG).

"Neither" rows are boundaries of another gene on the same strand that
`exonRefine` pulled into this gene (e.g. WDR77 rows at loci that belong to the
sense-overlapping lncRNA AL390195.1). The same locus appears again under its own
gene. They have to be dropped, not matched by locus alone.

## Overlapping genes

`WDR77` and the sense-overlapping `AL390195.1` produced 1,250 and 96
`Incompatible read ... separated by start site` warnings (a read path that
cannot be explained by the model; its weight is set to 1e-10). `ATP5PB`,
`C1orf162` and the single-end example gave 0. So this follows gene structure and
not the paired-end input. Sites in such genes need a flag.

## How the numbers are computed (and a limit)

`McSplicer.py` hardcodes `disable_end_sites_computation = True`:
- `s` (acceptors): EM estimate. Real values cover the whole range (26 internal
  acceptors in the test region: 10 at >= 0.99, 10 between 0.05 and 0.95, 3 below 0.05).
- `e` (donors): not EM. Value = spliced / (spliced + contiguous) over reads that
  cover the subexon, clamped to [1e-6, 1-1e-6]. An exon that ends in an intron
  has no contiguous continuation, so ordinary donors come out as exactly
  `0.999999`. In the test region 20 of 25 internal donors are `0.999999`, 4 are
  `1e-06` (3 with no reads, 1 with 413 contiguous and 0 spliced reads), 1 is
  intermediate (0.375, an alternative donor). **Donor labels from this version are
  close to binary**, which is why donors are not in the first regression dataset.
- Flipping the flag to `False` on a scratch copy crashes
  (`ValueError: 'e2' is not in list`, `cases_code_new.builder_end_10`), so the
  donor EM path is not usable as shipped.

## Read support and the `1e-06` ambiguity

Not in the output, but recoverable from the `.cnt`:
- Donor: `get_q_vector_from_junction_reads` returns `[site, key, prob, red,
  black]`; red = spliced count, black = contiguous count. Using it:
  `red = 0, black = 0` -> `1e-06` with **no reads**; `red = 0, black > 0` -> real
  zero spliced (transcript ends show 1-3996 black); `red > 0, black = 0` ->
  `0.999999` with 5-4157 supporting reads in the test region (the clamped 1).
- Acceptor: sum the `.cnt` counts of subpaths that enter the acceptor's subexon
  from a non-adjacent node (spliced-in) and from the adjacent node (contiguous).
  Test region: spliced-in counts range 0-4158 per acceptor (0-6 in the lowest
  expressed gene).
- Counts are mate-level (each fragment is counted once per mate) and use
  fractional multi-mapper weights, so they are not read-pair counts.

## Validation run (what was actually tested)

Region `chr1:111,420,000-111,480,000` (GRCh38), chosen because it holds ATP5PB
(`+`, 647 TPM in the ENCODE V29 quantification) and WDR77 (`-`, 73 TPM) with
overlapping exons, plus C1orf162 (`+`) and four small genes. ENCODE BAMs have no
`.bai`, so each replicate's region was located by binary search over the BGZF
blocks with HTTP range requests and the records re-written as a new BAM
(16.5 MB downloaded per replicate; BAM 1: ENCFF864OUS, 50,961 records; BAM 2:
ENCFF696CSX, 77,540 records; all paired, 100 bp, reads in the region only).
Annotation: the 7 GENCODE v29 genes lying fully inside the region. Timing on
this subset: sigcount 0.4 s, McSplicer (6 genes, no bootstraps) 1.3 s.

The checks are in `src/validation/` (run from the repo root with
`python -m src.validation.<name>`; they need `pysam`, and `curl`):

| module | what it does |
|---|---|
| `bam_region` | region of a remote coordinate-sorted BAM via range requests (no index needed) |
| `stream_fasta` | one pass over the gzipped FASTA: contig table + chosen contigs as plain text |
| `check_coordinates` | locus convention vs sequence (GT/AG at the locus, shifted controls, window helper) |
| `check_strand` | `fragment_strand` rule; strand concordance per gene and in antisense overlaps |
| `split_bam_by_strand` | the strand split used for preprocessing |
| `check_mcsplicer_output` | site classes, dinucleotides, donor/acceptor read support, per-site table |
| `compare_runs` | step-0 differences between two sets of McSplicer outputs |

Re-running them on the same data reproduced the numbers in this note (e.g. 51
internal sites all canonical, the strand split gives 37,701 `+` / 13,260 `-`
alignments in rep 1, max difference 0.0229 between the two runs). They were run on the subset only, not on a
full BAM.

## Open

- Donors: reporting only for now. If donor usage is wanted later it needs a different
  definition or a patched McSplicer; not attempted.
- The strand split on the full BAMs costs one sequential pass per BAM (not run).
- Still undecided, on purpose: coverage threshold, how to combine replicates,
  validation chromosomes. Test chromosomes are provisional.
- Runtime and memory of sigcount/McSplicer on a full 10-13 GB BAM: not measured.
- Which count to use for multi-mapped reads: McSplicer uses `totalCount` only.
- `--bootstraps` cost and how many are needed: not decided.
