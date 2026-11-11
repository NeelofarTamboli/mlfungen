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
   whether and how to combine them. This agreement is roughly the best score
   a model can reach.
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
- [ ] ENCODE file accessions chosen for replicate 1 and replicate 2
- [ ] Genome build and annotation version confirmed (deliberately **unset**
      in `config.yaml` until checked against the ENCODE file metadata)
- [ ] McSplicer installed and tested on its examples and on one chromosome
- [ ] McSplicer output format verified
- [ ] Everything after that (dataset, sequences, splits, baselines, CNNs)

**No data has been downloaded and there are no results yet.** The Python
modules contain only docstrings, TODOs, and `NotImplementedError` stubs.
They will be implemented once the input formats are checked against real
files.

## Repository structure

```
mlfungen/
├── config.yaml              # project settings (experiment, contexts, splits, paths)
├── requirements.txt
├── implementation_plan.md   # detailed plan incl. instructor feedback
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
│   └── make_splits.py       # chromosome-based splits
├── models/
│   ├── baseline.py          # mean / MaxEntScan / local-feature regression
│   └── cnn.py               # CNNs per context length (not implemented yet)
├── analysis/
│   ├── eda.py               # exploratory analysis and QC plots
│   └── evaluate.py          # metrics vs replicate ceiling
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

Check that the config loads:

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
