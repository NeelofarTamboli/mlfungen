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
- [ ] Sequences, splits (test chromosomes are provisional), baselines, CNNs

**Only chr1-4 have been processed in both replicates (39,382 clean matched acceptors), and there
are no model results yet.** `src/extract_sequences.py`, `src/make_splits.py`, `src/inspect_data.py`
and `models/` are still stubs with docstrings and TODOs.

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
`python -m src.download_data`. Right now every script stops with
`NotImplementedError` and a note on what has to be verified first.

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
