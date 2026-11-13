# Implementation Plan: Predicting Quantitative Splice-Site Usage from Sequence

> **Note (2026-11-11):** this plan was written before any input was checked.
> For now the project is K562 only, so every HepG2 item below is out of scope.
> GRCh38, the test chromosomes, the McSplicer steps and the coverage cutoff are
> still unverified or undecided. See `docs/mcsplicer_notes.md` and
> `docs/encode_inputs.md` for what has been checked.
>
> Decisions since then: the first regression target is internal **acceptor** usage
> only (McSplicer donor values are near-binary), and K562 reads are counted
> strand-aware (BAM split by transcript strand before sigcount). Read every
> "donor and acceptor" below with that in mind.

This plan is based on the proposal and the professor's comments. Items marked **[Prof]** come from the professor's feedback and are new compared to the proposal.

---

## Phase 0: Update the framing (do this first, it's quick)

- **[Prof]** Add **Pangolin** (Zeng & Li, *Genome Biology* 2022) to related work. Pangolin already predicts quantitative splice-site usage from sequence, using a dilated CNN like SpliceAI's. So "we predict usage, not just presence" isn't new by itself. What's new is:
  1. Which **types of sites** gain from long context, and which **distant sequence elements** the model uses to get there.
  2. Whether the amount of context needed **differs between cell types** (K562 vs HepG2).
- **[Prof]** Present the window sweep as an extension of **SpliceAI's 80 / 400 / 2k / 10k nt context comparison**. Use the **same window sizes** so the numbers line up directly. Add 5k or ~50 bp only if there is time.

## Phase 1: Data (target: weeks 1–2)

1. **[Prof]** **Download ENCODE's processed GRCh38 BAMs** for ENCSR000AEO (both replicates) instead of re-aligning reads from SRA. Also get the matching GENCODE GTF and the GRCh38 FASTA. They must be the same versions ENCODE used.
2. **[Prof]** **Pick a HepG2 poly(A)+ RNA-seq experiment** on the ENCODE portal from the same lab and protocol as the K562 one, and get its BAMs the same way.
3. **Run McSplicer** on each replicate and cell line. Budget time for its preprocessing; it needs the BAM plus the gene annotation in its own format. Output: a usage value per donor and per acceptor site.
4. **Filter by coverage and quality control:**
   - Choose a minimum read support per site (e.g., at least 10–20 junction reads).
   - Compute how well the two replicates agree (Pearson/Spearman). Replicate agreement provides a reproducibility benchmark for the RNA-seq-derived usage labels, so report it next to every result.
   - **[Prof]** **Report how many sites pass the filter in the interim report**, broken down by donor/acceptor and cell line.
5. **Build the dataset table:** chrom, position, strand, site type, usage per cell line, coverage, gene ID.
6. **Chromosome split:** use **SpliceAI's and Pangolin's test chromosomes (1, 3, 5, 7, 9)**. Otherwise the pretrained baselines will have seen the test sites during their training, which unfairly helps them. Hold out one or two other chromosomes for validation.
7. **Sequence extraction:** a one-hot encoder with strand handling (reverse-complement minus-strand sites) and N-padding at chromosome ends. Write it once so it can cut windows of any size.

## Phase 2: Baselines (target: by the interim)

- **[Prof]** **MaxEntScan scores**: 9-mer for donors, 23-mer for acceptors (the `maxentpy` package does this). Use them as features in the simple local baseline.
- **Local regression baseline:** linear regression or gradient boosting on MaxEntScan score, k-mer counts in roughly ±50 bp, GC content, and exon/intron lengths.
- **[Prof]** **Pretrained SpliceAI:** score each site with 10k context and correlate its splice-site probability with the usage values.
- **[Prof]** **Pretrained Pangolin** (`pip install pangolin`): extract its per-tissue usage outputs. It has no K562/HepG2 output, so try each of its tissues and/or the maximum. Make clear in the writeup that this is a cross-tissue comparison.

## Phase 3: Our models

- **Framework:** PyTorch. Start with a plain CNN on short windows. For 2k and 10k, use SpliceAI-style dilated residual blocks.
- **One output per site:** take the prediction at the center position.
- **Loss:** use binary cross-entropy on usage, like SpliceAI and Pangolin, rather than squared error. Optionally weight each site by its coverage, as the proposal already suggests.
- **Training:** train one model per window size (80 / 400 / 2k / 10k), with the same architecture family and a few random seeds each, so we get error bars.
- **Two cell lines:** a single network with **two output heads (K562, HepG2)**. That way the cell-type difference comes only from the labels, not the input sequence.
- The SSNet transformer comparison is a **stretch goal**. Drop it if time runs short; the professor's suggestions are worth more.

## Phase 4: Evaluation

- Metrics: R², MAE, Pearson, Spearman, reported **separately for donors and acceptors**, and compared against the replicate agreement.
- **Also evaluate only on alternative sites** (e.g., usage between 0.1 and 0.9). Most sites have usage near 1, which inflates correlations.

## Phase 5: Context analysis (the new contribution)

- **[Prof] Which types of sites gain from context:** for each site, measure the error improvement from the 400 nt model to the 10k model. Compare groups:
  - constitutive vs alternative sites;
  - sites with **another site of the same type nearby** (within roughly 50, 200 or 1000 bp) vs isolated sites;
  - short vs long exons and introns;
  - K562 vs HepG2.
- **[Prof] Attribution:** use in-silico mutagenesis (systematically mutating positions and watching the prediction change) and/or Integrated Gradients or DeepLIFT (the `captum` package) on the 10k model.
  - Add up the attribution in distance bins from the site to show how far away the important signal sits.
  - Run TF-MoDISco-lite on distant high-attribution regions to find motifs (e.g., branch points, ESE/ESS/ISE elements, RBP motifs).
  - Check whether the important distant regions are the competing nearby sites themselves.
  - Full mutagenesis over 10k positions is expensive, so run it on a subset of sites or mutate blocks of positions at a time.

---

## Interim report checklist

- [ ] BAMs downloaded (K562 + HepG2), McSplicer run, replicate agreement computed
- [ ] **Count of sites passing the coverage filter**, by site type and cell line
- [ ] Dataset table and chromosome split (test = chr 1, 3, 5, 7, 9)
- [ ] MaxEntScan/local baseline results
- [ ] SpliceAI and Pangolin baseline scores
- [ ] First CNN trained at 80 or 400 nt

## Main risks

- **McSplicer setup:** try it on a single chromosome in week 1.
- **Too few alternative sites passing the filter:** consider pooling the two replicates' reads.
- **GPU time for the 10k models:** get access to Colab, Google Cloud or the course cluster early.
