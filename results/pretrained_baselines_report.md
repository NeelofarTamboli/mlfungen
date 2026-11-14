# Pretrained baselines: SpliceAI and Pangolin (test chr1 + chr3)

**Status (14 Nov 2026, 03:35 UTC):** SpliceAI is finished and evaluated on every test site.
Pangolin scoring is still running (started 03:32 UTC, expected to finish around 06:00–06:15 UTC).
Its rows will be added to this file and to `pretrained_metrics_test_chr1_chr3.tsv` in a follow-up commit.

## Setup

- **Sites:** the same test rows as the existing baselines, all acceptors on chr1 and chr3 in
  `results/person_b_acceptors_chr1-4.tsv`. That is **23,474 sites** (chr1 14,844, chr3 8,630), of which
  **11,344** are alternative sites (target between 0.1 and 0.9).
- **Target:** mean McSplicer usage of the two K562 replicates, `(usage_rep1 + usage_rep2) / 2`,
  the same as `models/baseline.py --target mean`.
- **Metrics:** `analysis.evaluate.evaluate_predictions` (Pearson, Spearman, R², MAE, after dropping
  non-finite pairs), for all sites and for the alternative subset.
- **No leakage:** SpliceAI (Jaganathan et al. 2019) and Pangolin (Zeng & Li 2022) were both trained
  with chromosomes 1, 3, 5, 7 and 9 held out, so chr1 and chr3 are unseen test data for them too.
  The scorer refuses any other chromosome.
- **Pretrained models are used as is.** Nothing is fitted to our labels.

## What each score means

- **SpliceAI:** probability that the position is an acceptor splice site. It was trained on annotated
  (GENCODE) splice sites, so it measures how "splice-site-like" the site is, not how much it is
  used. It is not on the usage scale, so only Pearson and Spearman are reported (R² and MAE would mostly
  measure the calibration mismatch).
- **Pangolin usage (heart, liver, brain, testis):** predicted splice-site usage in that tissue, trained on
  RNA-seq from four species. It is on the same 0–1 scale as our target, so all four metrics apply.
  The mean of the four tissues is also reported. None of the tissues is K562.

Both models use 5,000 bases on each side of the site. Both peak exactly at our acceptor locus, the first
exonic base: the mean SpliceAI score is 0.816 at offset 0 and ≤ 0.004 at offsets −2..+2.

## Results

All sites (n = 23,474):

| Model | Score | Pearson | Spearman | R² | MAE |
|---|---|---|---|---|---|
| Training mean | usage | – | – | −0.001 | 0.271 |
| Ridge regression (local features) | usage | 0.252 | 0.170 | 0.063 | 0.261 |
| Gradient boosting (local features) | usage | 0.286 | 0.188 | 0.081 | 0.257 |
| SpliceAI | splice-site probability | **0.558** | **0.396** | n/a | n/a |
| Pangolin heart / liver / brain / testis / mean | usage | pending | pending | pending | pending |

Alternative sites, usage 0.1–0.9 (n = 11,344):

| Model | Score | Pearson | Spearman | R² | MAE |
|---|---|---|---|---|---|
| Training mean | usage | – | – | −0.196 | 0.195 |
| Ridge regression (local features) | usage | 0.188 | 0.143 | −0.156 | 0.192 |
| Gradient boosting (local features) | usage | 0.208 | 0.165 | −0.161 | 0.191 |
| SpliceAI | splice-site probability | **0.348** | **0.285** | n/a | n/a |
| Pangolin heart / liver / brain / testis / mean | usage | pending | pending | pending | pending |

The local baselines (MaxEntScan, GC, 2- and 3-mers in 100 bp) were trained on chr2 and validated on chr4.

## Interpretation so far

- SpliceAI's splice-site probability, with 10 kb of context and no K562 training, tracks K562 usage about
  twice as well as the local-feature models (Pearson 0.56 vs 0.29 on all sites).
- It also does better on the alternative sites (0.35 vs 0.21), but this subset is clearly harder for every model.
- Spearman is lower than Pearson for every model. Many sites have usage near 1, which compresses the ranks.

## How to compare the CNN results (for Ahn)

The `models.train_cnn` commands in the README (`--context 80/400 --target mean`) match this setup, but
`models/train_cnn.py` is not on GitHub yet, so it could not be checked directly. For a like-for-like comparison:

1. Train on chr2 and select the checkpoint on chr4, as in `results/split_settings_chr1-4.json`.
   The test chromosomes are chr1 and chr3.
2. Use `--target mean`, i.e. the mean of `usage_rep1` and `usage_rep2`.
3. Report test metrics with `analysis.evaluate.evaluate_predictions(predictions, [your_column], target_column="target", split="test")`.
   This gives Pearson, Spearman, R² and MAE for `all` (should be n = 23,474) and
   `alternative_0.1_to_0.9` (should be n = 11,344). If the counts differ, the test set is different.
4. Compare R² and MAE only with the usage models (the baselines and Pangolin usage). Compare Pearson and
   Spearman with all of them, including SpliceAI.
5. The pretrained models were trained on far more data (many tissues, full genome), so a CNN trained only
   on chr2 is not expected to beat them. Their scores are a reference for how much sequence-only models can explain.

## Files

- `results/pretrained_metrics_test_chr1_chr3.tsv`: the raw metrics (same columns as the tables above).
- `models/pretrained_baselines.py`: scorer and evaluation.
  ```
  python -m models.pretrained_baselines score --model spliceai --chrom chr1 chr3
  python -m models.pretrained_baselines score --model pangolin --chrom chr1 chr3
  python -m models.pretrained_baselines evaluate
  ```
- Per-site predictions are in `data/processed/pretrained/` (git-ignored, not committed).
- Runtime on 4 CPU cores in bfloat16: SpliceAI 0.19 s/site (75 min for all test sites); Pangolin about 0.4 s/site.
