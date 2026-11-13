# Acceptor usage labels (K562)

How the label table is made and what its columns mean. Results are in the last section.
McSplicer details (formats, coordinates, why strand-aware) are in `mcsplicer_notes.md`.

## Run it

```bash
# once: McSplicer (git-ignored), and the annotation
git clone https://github.com/canzarlab/McSplicer external/McSplicer     # tested at d6a22c3
python -m src.download_data gtf

# per replicate (rep1 = ENCFF864OUS, rep2 = ENCFF696CSX); ~10-13 GB each, md5-checked
python -m src.download_data bam rep1
python -c "import pysam; pysam.index('data/raw/ENCFF864OUS.bam')"      # ENCODE ships no .bai
python -m src.run_mcsplicer --bam data/raw/ENCFF864OUS.bam --out data/interim/mcsplicer/rep1 --bootstraps 0 --jobs 3
python -m src.parse_mcsplicer --rep-dir data/interim/mcsplicer/rep1    # -> data/processed/acceptors_rep1.tsv

# after both replicates
python -m src.build_dataset                                            # -> data/processed/acceptor_labels.tsv
python -m analysis.eda --sites-all data/interim/mcsplicer/rep1/sites_all.tsv data/interim/mcsplicer/rep2/sites_all.tsv
python -m analysis.coverage
```

`run_mcsplicer` is resumable (a chromosome with a `done.json` is skipped) and keeps all raw McSplicer
output under `data/interim/` (git-ignored). A small region can be tried first:
`python -m src.validation.bam_region URL chr1 111420000 111480000 region.bam --index` and then
`run_mcsplicer --bam region.bam --region chr1:111420000-111480000`.

## What the pipeline does

For each chromosome: split reads by fragment strand (read 1 reverse or read 2 forward = `+` gene),
then per strand `exonRefine` -> `bin/sigcount` -> `McSplicer.py --read_len 100` on that strand's GENCODE v29
genes. Genes are run in chunks of 100 with a trimmed `.cnt`, because McSplicer re-reads the whole `.cnt`
for every gene (a test with 2-gene chunks gave identical step-0 values). The parser keeps only `s`
rows that are internal acceptors by the GTF rule in `mcsplicer_notes.md` and recovers read support.

## Label table columns (`data/processed/acceptor_labels.tsv`)

One row per internal acceptor; sites seen in only one replicate are kept and flagged.

| column | meaning |
|---|---|
| `site_id` | `chrom:locus:strand:gene_id` |
| `gene_id`, `chrom`, `locus`, `strand` | `locus` is McSplicer's 1-based first exonic base (see notes for sequence conversion) |
| `site_class` | `acceptor` |
| `also_transcript_start` | internal acceptor that is also the first-exon start of another transcript (kept, flagged) |
| `usage_rep1/2` | McSplicer step-0 estimate for that replicate |
| `bootstrap_sd_rep1/2`, `n_bootstraps_rep1/2` | sample SD over bootstrap steps; **empty where no bootstraps were run** (see below) |
| `read_support_rep1/2` | `spliced_in + contiguous_in` |
| `spliced_in_*`, `contiguous_in_*` | `.cnt` counts entering the acceptor's subexon from a non-adjacent / adjacent subexon |
| `n_incompat_warnings_*` | McSplicer "Incompatible read" warnings for the gene (gene-level flag) |
| `has_read_support_rep1/2` | `read_support > 0` |
| `in_rep1`, `in_rep2`, `in_both` | which replicates have the site |
| `in_primary_qc` | in both replicates, nonzero support in both, not `also_transcript_start`. The main QC set; no other threshold |
| `index_*`, `raw_csv_*` | McSplicer's `s<k>` and the raw per-gene CSV the row came from |

Read support counts are per mate (a read pair adds to both mates) and include the fractional weights
McSplicer gets from multi-mapped reads, so they are not read-pair counts.

Replicates are **not averaged** and nothing is filtered on read support.

## Things to know before using the labels

- **Usage without support is not a measurement.** An acceptor in a gene that has reads but with no reads
  at the acceptor still gets a value; it sits at about 0.5 (McSplicer's starting value). On chr3, replicate 1,
  29% of acceptors had zero support and their median usage was 0.4997.
- **`also_transcript_start` sites are kept but flagged.** A locus that is an internal acceptor in one transcript
  and the first exon start of another gets one McSplicer value that mixes both events, so the primary QC leaves
  them out for now. On chr1-4 there are 814 (589 with nonzero support in both replicates); their usage and
  replicate agreement look like the other acceptors (see below). Nothing has been removed from the data.
- **McSplicer crashes on some genes** (an `IndexError` in its own code). The runner catches this, skips the gene,
  carries on with the rest of the chunk and writes `failed_genes.json` (per chunk, and merged in the replicate
  directory); it also prints them. Reruns resume without recomputing finished genes. Known so far: replicate 1
  ENSG00000183474.15 (chr5) and ENSG00000273066.5 (chr9) crash it. None failed on replicate 2 chr1-4. Failed genes
  have no rows in the tables.
- **Bootstraps are only on part of the genome.** Ten bootstraps cost about 11x the time and the heavy genes
  dominate (top 20 genes of chr22 = 78% of the time), so the genome-wide run used none. `bootstrap_sd` is
  only filled for chromosomes that were run with bootstraps (`n_bootstraps_*` says which).
- Genes with many "Incompatible read" warnings (usually same-strand overlapping genes) should be treated with care.

## Interim QC, chr1-4 (both replicates, no bootstraps)

Run: `python -m src.build_dataset --chroms chr1 chr2 chr3 chr4 --out data/processed/acceptor_labels_chr1-4.tsv`, then
`python -m analysis.eda --labels data/processed/acceptor_labels_chr1-4.tsv --tag chr1-4 --bootstrap-table ...`.
Replicate 2 chr1-4 took 24 min on 3 cores. Numbers are in `results/replicate_qc_chr1-4.json`, figures A-D in `figures/qc_[A-D]_*_chr1-4.png` (E, the bootstrap plot, is exploratory and kept local).

| | sites |
|---|---|
| internal acceptors, union of replicates | 62,557 (rep1 59,384, rep2 60,060) |
| matched in both replicates | 56,887 |
| zero read support: rep1 / rep2 | 17,112 / 16,921 |
| matched, nonzero support in both | 39,971 |
| flagged `also_transcript_start` (all / matched / nonzero in both) | 814 / 749 / 589 |
| primary QC set | 39,382 |

Replicate agreement on the primary set: Pearson 0.905, Spearman 0.890, median |diff| 0.029.

This is a replicate-reproducibility (noise) benchmark: it says how consistent the labels are between two biological
replicates of the same cells. It is a useful reference to report next to model results, but it is not a strict
ceiling on model performance.

- Zero-support sites do not lower the correlation, they inflate it and add noise. 11,723 matched sites have zero reads
  in both replicates and sit at about 0.5 in both (they "agree"); 5,033 have zero in only one (Pearson 0.64).
- Agreement rises steadily with read support (median |diff| 0.058 at support 1, 0.029 at 64-127, 0.008 at >=2048;
  share of sites differing by more than 0.1: 40%, 19%, 3%) with no plateau (`results/support_bins_chr1-4.tsv`).
- Bootstrap SD was only run on a 60 kb region of chr1 (23 primary sites). Plot E is exploratory and is not committed.
- 13% of primary sites are in genes with "Incompatible read" warnings; they agree as well as the rest.

## Coverage threshold: not chosen

chr1-4 is enough to show that higher read support gives more reproducible labels (the support-bin table and
plots C and D). It is not used to pick the final support threshold:

- chr1 and chr3 are provisional test chromosomes, so a cutoff tuned on them would leak held-out information.
- The final minimum read support will be selected using training/validation chromosomes only, never the held-out
  test chromosomes. `labels.min_read_support` in `config.yaml` stays `null` until then.
- How the cutoff is chosen (a target level of replicate agreement, or weighting by support instead of a hard cut)
  is still open.

## Handoff to the sequence-extraction step (Person B)

`results/person_b_acceptors_chr1-4.tsv` is the input for the sequence windows: the 39,382 chr1-4 sites in the
primary QC set (both replicates, nonzero support in both, not `also_transcript_start`). Columns: `site_id, gene_id,
chrom, locus, strand, usage_rep1, usage_rep2, read_support_rep1, read_support_rep2, also_transcript_start,
in_primary_qc`. Usage is rounded to 6 decimals and support to 2. Rebuild it with
`python -m src.build_dataset --chroms chr1 chr2 chr3 chr4 --out data/processed/acceptor_labels_chr1-4.tsv --handoff results/person_b_acceptors_chr1-4.tsv`.

- `locus` is McSplicer's 1-based first exonic base of the acceptor. On `+` the junction index is `j = locus - 1`,
  on `-` it is `j = locus` (see "Coordinates" in `mcsplicer_notes.md` for the window rule and checks).
- No read-support threshold has been applied, and the replicates are not averaged; keep both columns.
- chr1 and chr3 are provisional test chromosomes.
