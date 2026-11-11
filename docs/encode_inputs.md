# ENCODE inputs for K562 (ENCSR000AEO)

Checked 2026-11-11 from the ENCODE portal JSON (`/experiments/ENCSR000AEO/`,
`/files/<acc>/`) and GEO (`GSM2072366`). Nothing large was downloaded. For BAM
headers I range-requested the first 2 MB of each candidate file; for the FASTA
and GTF, the first 64 KB / 256 KB. Everything below is what the portal or those
bytes said; where I did not check something it says so.

## Experiment

- ENCSR000AEO, polyA plus RNA-seq, K562, status released, lab Brenton
  Graveley (UConn). GEO series GSE78559.
- Reads: 2 x 100, paired-ended, Illumina HiSeq 2000. Four FASTQs, two per
  replicate.
- Both libraries are annotated `strand_specificity: reverse`.

| bio rep | library | GEO | biosample | FASTQs |
|---|---|---|---|---|
| 1 | ENCLB059ZZZ | GSM2072366 | ENCBS087RNA | ENCFF001RDE, ENCFF001RCW |
| 2 | ENCLB060ZZZ | GSM2072367 | ENCBS088RNA | ENCFF001RDD, ENCFF001RCV |

GEO `GSM2072366` title: "polyA mRNA RNA-seq from K562 (ENCLB059ZZZ)", replicate 1.
`GSM2072367` was taken from the ENCODE library record only (GEO page not opened).

## Candidate genomic alignment BAMs

Only released `alignments` BAMs are listed. The experiment has 10 more
`alignments` BAMs with status `archived` (older runs), and separate
`transcriptome alignments` BAMs; neither is a candidate. ENCODE lists no `.bai`
files for these; we would index them ourselves.

| file | rep | assembly | annotation | size | pipeline | analysis status |
|---|---|---|---|---|---|---|
| ENCFF864OUS | 1 | GRCh38 | V29 | 10.39 GB | bulk-rna-seq alignment, pipeline 1.2.1, STAR 2.5.1b, samtools 1.9 | released (ENCAN703TGC) |
| ENCFF696CSX | 2 | GRCh38 | V29 | 13.35 GB | same | released (ENCAN703TGC) |
| ENCFF846WOV | 1 | GRCh38 | V24 | 10.39 GB | lrna-align-star-pe 2.1, STAR 2.5.1b | archived (ENCAN800IJQ) |
| ENCFF588YLF | 2 | GRCh38 | V24 | 13.35 GB | same | archived (ENCAN800IJQ) |
| ENCFF832OIM | 1 | hg19 | V19 | 10.30 GB | lrna-align-star-pe 2.0, STAR 2.4.0k | archived (ENCAN969HWK) |
| ENCFF831GPG | 2 | hg19 | V19 | 13.17 GB | same | archived (ENCAN969HWK) |

All six: file status `released`, file format BAM, output type `alignments`,
mapped run type paired-ended, mapped read length 100. md5 sums are in the portal
records (V29: `c941a6b1...` rep 1, `0222c73b...` rep 2).

Both V29 BAMs came from one workflow run (same
`caper_out_v02_21/rna/e16e2d68-...` path, shard-0 = rep 1, shard-1 = rep 2) and
list the same STAR index `ENCFF598IDH` as input.

STAR QC (V29): rep 1 57.2M input reads, 74.7% uniquely mapped; rep 2 72.1M input
reads, 83.6% uniquely mapped. Rep 1 has fewer reads and lower mapping rate than
rep 2.

### BAM header facts (the four GRCh38 BAMs)

- `@HD ... SO:coordinate`: coordinate-sorted.
- 292 reference sequences: `chr1` (248956422) ... plus `chrM`, other `chr*`
  contigs (195 in total start with `chr`), and `ERCC-*`, `phiX174` and a few
  others spike-in/extra contigs (not enumerated).
- So the BAMs use `chr1`-style names. Not the Ensembl `1`-style.
- STAR attributes are `NH HI AS NM MD`. No `XS` tag (`--outSAMstrandField`
  was not set), so sigcount `-s` has nothing to read.
- The V24 BAM comments say `REFID:GRCh38_no_alt_analysis_set_GCA_000001405.15`
  and `ANNID:gencodeV24pri-tRNAs-ERCC-phiX`. The V29 BAM header has no
  REFID/ANNID comment lines; its reference comes from the listed STAR index
  `ENCFF598IDH`, whose inputs are the V29 GTF, V29 tRNAs, the no_alt analysis
  set FASTA and the spike-ins.

## Reference files that match the V29 BAMs

| role | ENCODE file | size | note |
|---|---|---|---|
| gene annotation | `gencode.v29.primary_assembly.annotation_UCSC_names` (.gtf.gz) | 42.7 MB | seqnames are `chr1` style; `gene_id` has versions, e.g. `ENSG00000223972.5`; header says GENCODE v29 / Ensembl 94 |
| genome FASTA | `GRCh38_no_alt_analysis_set_GCA_000001405.15` (.fasta.gz) | 873 MB | first record `>chr1 ... LN:248956422 ... AS:GRCh38`; md5 `a08035b6...` |

- Contig naming agrees across BAM, GTF and FASTA for chr1 (the only one
  checked in the GTF/FASTA). I did not compare the full contig lists, and I do not
  have a `.fai`.
- The STAR index also used the GENCODE v29 tRNA GTF and the spike-ins. That
  tRNA GTF is not needed for McSplicer unless we decide otherwise.
- For V24 the matching GTF is `gencode.v24.primary_assembly.annotation`
  (38.7 MB); its seqname style was not checked. For hg19 the index lists
  `female.hg19` and `gencode.v19.annotation`; I did not look at those records.

## Inputs selected for the project

| item | value |
|---|---|
| replicate 1 BAM | ENCFF864OUS (library ENCLB059ZZZ, GEO GSM2072366), 10.39 GB, md5 `c941a6b12cce8d8dd4fd3fd1ff29e492` |
| replicate 2 BAM | ENCFF696CSX (library ENCLB060ZZZ, GEO GSM2072367), 13.35 GB, md5 `0222c73b629831008bf0d278529e6ef2` |
| assembly | GRCh38 |
| annotation | GENCODE v29, `gencode.v29.primary_assembly.annotation_UCSC_names` (.gtf.gz, 42,720,023 bytes, md5 `7f3e84122de9f5a56ca4a17b0a44a33c`) |
| FASTA | `GRCh38_no_alt_analysis_set_GCA_000001405.15` (.fasta.gz, 872,949,833 bytes, md5 `a08035b6a6e31780e96a34008ff21bd6`) |

The md5 values are the ones ENCODE lists; I re-computed only the GTF's (the FASTA
was checked by size and per-contig `M5`, the BAMs not at all since they are not
downloaded). This combination passed the compatibility checks below (contig names, lengths
and order; matching headers and pipeline for the two BAMs; GTF md5 and FASTA
`M5` values; the same tools ran on both BAMs). Total download is about 24 GB of
BAM plus 0.9 GB of reference; the full BAMs have not been downloaded yet.

Why not the others: the V24 and hg19 BAMs belong to `archived` analyses; V29 is
the only released analysis. hg19 is the only other clean pair, if GRCh38 ever
becomes a problem. These choices are recorded in `config.yaml`.

## Compatibility check: BAM, GTF, FASTA (2026-11-11)

Done with the full FASTA (streamed once, 872,949,833 bytes = the size ENCODE
lists, nothing kept except chr1/chr2), the full GTF (md5 `7f3e8412...` matches
ENCODE), and the headers of ENCFF864OUS and ENCFF696CSX.

- **FASTA**: 195 contigs. Headers carry `AS:GRCh38`, `LN:` and `M5:`. The md5 of
  the streamed chr1 and chr2 sequences equals the `M5` in their headers. chr1 and
  chr2 have no lower-case bases, but contain IUPAC codes besides N (chr1: M, R;
  chr2: K, M, R, W, Y), which the one-hot encoder must handle. Other chromosomes
  not checked.
- **BAM vs FASTA**: the BAM header has 292 contigs. The 195 FASTA contigs are all
  in it with identical names, lengths and order: 25 primary (chr1-22, X, Y, M),
  42 `_random`, 127 `chrUn_*` and `chrEBV`. The other 97 are spike-ins: 96 `ERCC-*`
  and `phiX174`. They have no counterpart in the GTF, so no reads are counted
  there.
- **GTF vs FASTA/BAM**: the GTF uses 47 seqnames, all present in both the FASTA
  and the BAM.
- **Naming**: `chr1` style everywhere. No renaming needed.
- **Build**: chr1 length and the FASTA headers say GRCh38; GENCODE v29 is
  GRCh38-based. I did not compare every one of the 195 `M5` values with an
  independent source.
- **Rep 1 vs rep 2 BAMs**: identical contig lists, and the non-`@SQ` header lines
  are identical after masking the input paths: same STAR 2.5.1b command line
  (`--outSAMattributes NH HI AS NM MD`, `--outSAMtype BAM SortedByCoordinate`,
  `--outFilterType BySJout`, `--quantMode TranscriptomeSAM`, 32 threads), same
  workflow run, same STAR index (ENCFF598IDH).
- **Handled identically**: I extracted the same chr1 region (60 kb) from both with the
  same tools. Both are 100 bp, every record paired and flagged proper pair, same
  tags, same reverse-stranded orientation (100% concordant on two genes), no `XS`
  tag. The same sigcount and McSplicer commands ran on both without changes.
  Differences are only depth and mapping rate (rep 1: 57.2M reads, 74.7% unique;
  rep 2: 72.1M, 83.6%).

## Not resolved here

- The BAMs carry no `XS` tag and sigcount ignores library strand: the project
  splits each BAM by fragment strand before counting (`mcsplicer_notes.md`,
  "Decisions").
- Runtime and memory of sigcount on a 10-13 GB BAM were not measured (the 51k-record
  subset took 0.4 s).
- Donor labels are not usable as a regression target in this McSplicer release, so
  the first dataset is internal acceptors only (`mcsplicer_notes.md`).
- The hg19 and GRCh38 V24 reference files were not inspected.

## Assumptions in the repo before this check

| where | assumption | status |
|---|---|---|
| `implementation_plan.md`, Phase 1 | use ENCODE's processed GRCh38 BAMs | GRCh38 BAMs exist for both replicates (released V29 pair) and passed the compatibility checks; `genome_build: GRCh38` is now set in `config.yaml` |
| `implementation_plan.md` (Phases 1-5) | HepG2 as a second cell line / two output heads / K562 vs HepG2 analysis | out of scope for now; not touched in the plan text beyond a banner |
| `implementation_plan.md` Phase 1.6, `config.yaml` `splits` | test chromosomes 1,3,5,7,9 (SpliceAI/Pangolin) | provisional. I did not check which chromosomes SpliceAI or Pangolin held out. We also do not know how many usable sites each chromosome will have |
| `implementation_plan.md` Phase 1.3 | McSplicer gives "a usage value per donor and per acceptor site" | McSplicer reports `s`/`e` boundaries incl. transcript starts/ends, and has no read-support column; donor values are near-binary. First regression target is internal acceptors only (`mcsplicer_notes.md`) |
| `implementation_plan.md` Phase 1.4 | e.g. 10-20 junction reads minimum | example only; no cutoff chosen |
| `implementation_plan.md` risks | pool the two replicates' reads | no decision on combining replicates |
| `config.yaml` | validation chromosomes unset | still unset |
