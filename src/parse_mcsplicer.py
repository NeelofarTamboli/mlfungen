"""Parse McSplicer output into a tidy per-site usage table.

Format (verified on the bundled example and on a real K562 subset): one CSV per gene,
`<out_dir>/<chr>/<gene_id>[_<prefix>].csv`, header
`step,index,strand,chr,splice site locus,usage estimate`. Step 0 is the estimate from the
observed counts; steps 1..N are bootstrap estimates. `index` is `s<k>` (start site, an
acceptor when internal) or `e<k>` (end site, a donor when internal). The locus is the
1-based GTF coordinate of an exonic base (see docs/mcsplicer_notes.md for the exact
conversion on each strand). There is no read-count column.

Planned output (one row per site per replicate): chrom, locus, strand, gene_id, index,
site class, usage (step 0), bootstrap SD, read support, intron dinucleotide.

`src/validation/check_mcsplicer_output.py` already does the site classification,
dinucleotide check and read-support recovery on a subset; reuse it here.

TODO:
    * Move that logic into this module, run it on the full outputs.
    * Keep only site_class == "acceptor" rows for the regression dataset
      (config `labels.target_site_type`); keep donor rows for QC reporting.
    * Flag rows from genes with "Incompatible read" warnings (overlapping
      same-strand genes) and `neither`-class rows (drop).
"""

from __future__ import annotations


def main() -> None:
    """Command-line entry point. Not implemented until the format is known."""
    raise NotImplementedError("Not implemented yet; the format and plan are in the module docstring.")


if __name__ == "__main__":
    main()
