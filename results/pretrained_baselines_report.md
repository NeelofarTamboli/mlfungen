# Pretrained baselines: SpliceAI and Pangolin (test chr1 + chr3)

Both models are scored on every test site. SpliceAI results are unchanged from the first commit;
the Pangolin rows were added once its run finished.

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
exonic base: the mean SpliceAI score is 0.816 at offset 0 and ≤ 0.004 at offsets −2..+2, and the mean
Pangolin usage is 0.73–0.75 at offset 0 and ≤ 0.001 at the other offsets
(`results/pretrained_offsets_test_chr1_chr3.tsv`). The offset was picked from the mean score only, without the labels.

## Results

All sites (n = 23,474):

| Model | Score | Pearson | Spearman | R² | MAE |
|---|---|---|---|---|---|
| Training mean | usage | – | – | −0.001 | 0.271 |
| Ridge regression (local features) | usage | 0.252 | 0.170 | 0.063 | 0.261 |
| Gradient boosting (local features) | usage | 0.286 | 0.188 | 0.081 | 0.257 |
| SpliceAI | splice-site probability | **0.558** | **0.396** | n/a | n/a |
| Pangolin heart | usage | 0.595 | 0.435 | 0.183 | 0.195 |
| Pangolin liver | usage | 0.594 | 0.430 | 0.218 | 0.195 |
| Pangolin brain | usage | 0.599 | 0.436 | 0.195 | 0.195 |
| Pangolin testis | usage | 0.598 | 0.438 | **0.231** | **0.193** |
| Pangolin mean of 4 tissues | usage | **0.601** | **0.442** | 0.223 | **0.193** |

Alternative sites, usage 0.1–0.9 (n = 11,344):

| Model | Score | Pearson | Spearman | R² | MAE |
|---|---|---|---|---|---|
| Training mean | usage | – | – | −0.196 | 0.195 |
| Ridge regression (local features) | usage | 0.188 | 0.143 | −0.156 | 0.192 |
| Gradient boosting (local features) | usage | 0.208 | 0.165 | −0.161 | 0.191 |
| SpliceAI | splice-site probability | **0.348** | **0.285** | n/a | n/a |
| Pangolin heart | usage | 0.381 | 0.330 | −1.126 | 0.268 |
| Pangolin liver | usage | 0.378 | 0.325 | −0.941 | 0.250 |
| Pangolin brain | usage | 0.386 | 0.333 | −1.085 | 0.264 |
| Pangolin testis | usage | 0.385 | 0.334 | −0.940 | 0.251 |
| Pangolin mean of 4 tissues | usage | **0.387** | **0.338** | −0.986 | 0.255 |

The local baselines (MaxEntScan, GC, 2- and 3-mers in 100 bp) were trained on chr2 and validated on chr4.

## Interpretation

- Both pretrained models, with 10 kb of context and no K562 training, track K562 usage about twice as well
  as the local-feature models (Pearson 0.56 for SpliceAI and 0.60 for Pangolin vs 0.29 for gradient boosting on all sites).
- Pangolin usage is slightly better than SpliceAI on every correlation. The four tissues are within 0.01 of each other,
  and their mean is the best by a small margin, so the choice of tissue hardly matters for K562.
- On all sites Pangolin also has the best R² (0.18–0.23) and MAE (0.193–0.195, vs 0.257 for gradient boosting).
- The alternative sites are harder for every model (Pearson 0.35–0.39 for the pretrained models, 0.21 for gradient boosting).
- On the alternative sites Pangolin ranks sites better but has a much worse R² (about −1) and MAE (0.25–0.27) than the
  training mean (0.195). Its predictions for these sites still spread over the whole 0–1 range (SD 0.29) and average 0.75,
  while the targets average 0.61 with SD 0.23. So its usage values are not calibrated to K562 at intermediate sites.
  R² on this subset is measured against the subset's own variance, which is why the local models are also below 0.
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
   Spearman with all of them, including SpliceAI. The main reference rows are gradient boosting (best local model)
   and Pangolin mean of 4 tissues (best pretrained model).
5. The pretrained models were trained on far more data (many tissues, full genome), so a CNN trained only
   on chr2 is not expected to beat them. Their scores are a reference for how much sequence-only models can explain.

## Files

- `results/pretrained_metrics_test_chr1_chr3.tsv`: the raw metrics (same columns as the tables above).
- `results/pretrained_offsets_test_chr1_chr3.tsv`: mean score of each model output at offsets −2..+2 around the locus.
- `models/pretrained_baselines.py`: scorer and evaluation.
  ```
  python -m models.pretrained_baselines score --model spliceai --chrom chr1 chr3
  python -m models.pretrained_baselines score --model pangolin --chrom chr1 chr3
  python -m models.pretrained_baselines evaluate
  ```
- Per-site predictions are in `data/processed/pretrained/` (git-ignored, not committed).
- Runtime on 4 CPU cores in bfloat16: SpliceAI 0.19 s/site (75 min for all test sites); Pangolin about 0.5 s/site (about 3.3 h for all test sites).
