# mlfungen

COMS W4762 (Machine Learning for Functional Genomics, Columbia) class project:
predicting quantitative splice-site usage from genomic sequence.

## Biological question

Can the genomic sequence around a splice site predict how much that site is
used (quantitative donor/acceptor usage) in K562 cells? And how much sequence
context is needed for an accurate prediction?

## Planned pipeline

1. **Data.** K562 RNA-seq from ENCODE experiment **ENCSR000AEO**, two
   biological replicates (processed alignments), plus the matching reference
   genome FASTA and gene annotation GTF.
2. **Labels.** Run [McSplicer](https://github.com/canzarlab/McSplicer)
   on each replicate to estimate usage for each donor and acceptor site.
3. **Replicate QC.** Measure how well the replicates agree before deciding
   whether and how to combine them. This agreement is a replicate-reproducibility
   (noise) benchmark for the labels, to report next to model results; it is not a
   strict limit on model performance.
4. **Dataset.** Filter by read support and build a site table
   (chrom, position, strand, site type, usage, coverage, gene).
5. **Sequences.** Extract windows around each site for each context
   length: 80 bp, 400 bp, 2 kb, 10 kb (the same lengths SpliceAI compared).
6. **Splits.** Chromosome-based train/validation/test split. The planned test
   chromosomes are chr1, 3, 5, 7, 9, the ones SpliceAI and Pangolin held out.
7. **Baselines first.** Mean predictor, then local-feature regression
   (MaxEntScan score, k-mers, GC content, exon/intron lengths).
8. **CNNs.** One model per context length, from the same architecture family.
9. **Comparisons (if practical).** MaxEntScan, Pangolin (and SpliceAI).
10. **Context analysis.** Which splice sites benefit from long-range
    sequence, and why.
11. *Possible later extension:* AlphaGenome. Not required now.

More detail is in [`implementation_plan.md`](implementation_plan.md).

## Current status

- [x] Repository structure, config, `.gitignore`, and requirements
- [x] McSplicer checked from its source and run on its example and on a real K562
      paired-end subset (`docs/mcsplicer_notes.md`)
- [x] ENCODE inputs chosen and checked: BAMs ENCFF864OUS / ENCFF696CSX, GRCh38,
      GENCODE v29, `GRCh38_no_alt_analysis_set_GCA_000001405.15`
      (`docs/encode_inputs.md`, `config.yaml`)
- [x] Decided: first regression target is internal **acceptor** usage only (McSplicer
      donor values are near-binary); K562 reads are counted strand-aware
- [x] Strand-aware McSplicer runner, internal-acceptor parser with read support, and label
      table builder (`docs/usage_labels.md`)
- [x] Both BAMs downloaded; chr1-4 processed in both replicates; interim replicate QC
      (`results/replicate_qc_chr1-4.*`, `figures/qc_[A-D]_*_chr1-4.png`)
- [ ] The rest of the genome for both replicates
- [ ] Read-support threshold and replicate combination (deliberately undecided; the threshold
      will be chosen on training/validation chromosomes only)
- [x] GRCh38 sequence extraction for the 39,382 chr1-4 handoff sites, including
      strand orientation, four centred contexts, N handling and lazy one-hot encoding
- [x] Chromosome splitter, MaxEntScan/GC/k-mer baselines, and test-set regression metrics
- [ ] Final full-genome splits and baseline results (test chromosomes are provisional), CNNs

**Only chr1-4 have been processed in both replicates (39,382 clean matched acceptors), and there
are no final full-genome model results yet.** The current chr1-4 baseline is preliminary: it uses
chr2 for training, chr4 for validation, chr1/chr3 for testing, and the mean of the two replicate
usage estimates as an explicitly selected target. The final analysis must be rerun after the remaining
chromosomes, support threshold and replicate-combination rule are decided.

## Repository structure

```
mlfungen/
├── config.yaml              # project settings (experiment, contexts, splits, paths)
├── requirements.txt
├── implementation_plan.md   # detailed plan incl. instructor feedback
├── docs/
│   ├── mcsplicer_notes.md   # verified McSplicer behaviour + project decisions
│   └── encode_inputs.md     # ENCODE files, reference files, compatibility checks
├── data/                    # NOT committed (only the empty folders are)
│   ├── raw/                 #   downloaded BAM / FASTA / GTF
│   ├── interim/             #   McSplicer intermediate and raw outputs
│   └── processed/           #   final site tables, sequences, splits
├── src/
│   ├── config.py            # loads config.yaml
│   ├── download_data.py     # fetch ENCODE BAMs + reference files
│   ├── inspect_data.py      # sanity-check headers, contigs, read length
│   ├── run_mcsplicer.py     # exonRefine -> sigcount -> McSplicer
│   ├── parse_mcsplicer.py   # McSplicer CSV -> tidy per-site table
│   ├── build_dataset.py     # replicate agreement, filtering, labels
│   ├── extract_sequences.py # windows around each site, one-hot encoding
│   ├── make_splits.py       # chromosome-based splits
│   └── validation/          # small checks used to verify McSplicer and the inputs
├── models/
│   ├── baseline.py          # mean / MaxEntScan / local-feature regression
│   └── cnn.py               # CNNs per context length (not implemented yet)
├── analysis/
│   ├── eda.py               # exploratory analysis and QC plots
│   └── evaluate.py          # metrics vs replicate benchmark
├── figures/
├── results/
└── notebooks/
```

## Setup for collaborators

```bash
git clone https://github.com/NeelofarTamboli/mlfungen.git
cd mlfungen
```

```bash
python3 -m venv .venv
```

```bash
source .venv/bin/activate
```

```bash
pip install -r requirements.txt
```

`src/validation` also needs `curl`. Check that the config loads:

```bash
python -c "from src.config import load_config; print(load_config()['dataset'])"
```

Run scripts from the repository root as modules, e.g.
`python -m src.download_data`. Some later-stage scripts still stop with
`NotImplementedError` until their required inputs or decisions are available.

### Person B: sequences and preliminary baselines

```bash
python -m src.run_part_b
```

That single, resumable command downloads the reference if needed, extracts the
sequences, makes the chromosome split, fits the local baselines and writes the
evaluation metrics. Existing stage outputs are skipped; use `--force` to rerun
them. To run selected stages only, use, for example,
`python -m src.run_part_b --steps split baseline --force`.

When the full-genome handoff arrives, use a new tag so the preliminary chr1-4
artifacts are not overwritten:

```bash
python -m src.run_part_b \
  --sites results/person_b_acceptors_all.tsv \
  --tag all \
  --validation chr2 chr4 \
  --target mean
```

Train the minimal short-context CNNs after the sequence and split artifacts exist:

```bash
python -m models.train_cnn --context 80 --target mean
python -m models.train_cnn --context 400 --target mean
```

Each command selects the best checkpoint by validation MAE and writes the model
to `results/cnn_<context>bp_<target>.pt` plus a JSON validation/test summary.

The extractor stores only the largest 10 kb context as compact integer nucleotide
codes. The 80, 400 and 2,000 bp contexts are centred views, and one-hot arrays are
created per batch.

- **Push to GitHub:** source code, tests, README, and compact files such as
  `results/sequence_qc_*.json`, `results/split_counts_*.tsv`, and
  `results/baseline_metrics_*.tsv`.
- **Do not push:** the FASTA, prepared chromosomes, sequence arrays, feature
  matrices, split tables, or full prediction tables. They are reproducible,
  stored under `data/`, and already ignored by Git.

External tools (not pip packages) needed later: McSplicer (clone it into
`external/`, which is git-ignored) and samtools.

## Data rule: never commit biological data to Git

- BAM/CRAM/SAM, FASTQ, FASTA, GTF/GFF, index files, genomic tracks, large
  arrays, and model checkpoints are all in `.gitignore`. Everything under
  `data/` is ignored except the empty folder placeholders.
- Keep data in `data/raw`, `data/interim`, `data/processed`. Share it through
  the download scripts (with recorded accessions and checksums) or shared
  storage, not through Git.
- Never commit credentials, API keys, or tokens.
- Before committing, run `git status` and make sure no data files are staged.
  Do not use `git add -f` on data.
