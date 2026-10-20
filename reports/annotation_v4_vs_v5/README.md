# RatGTEx v4 versus v5 reference annotations

Both releases use GRCr8 (GCF_036323735.1). The GTF headers identify v4 as NCBI RefSeq RS_2024_02 (2024-02-09) and v5 as RS_2026_08 (2026-08-17).

Run `python reports/annotation_v4_vs_v5/compare.py` from the repository root to regenerate `counts.tsv`, `gene_changes.tsv`, and `provenance.json`.

Counts use the actual chromosome-filtered `.chr.gtf` files named in `provenance.json`, before expression filtering or gene collapsing. Genes are unique NCBI GeneIDs; transcript counts use unique versioned transcript IDs from transcript features. Exons are distinct chromosome/start/end/strand intervals, pooled across transcripts and genes. Five extra gene feature records in v4 represent repeated GeneIDs at separate loci. The parser checks that repeated GeneIDs have consistent biotypes and repeated transcript IDs have consistent parent GeneIDs.

The shared/only comparisons use NCBI GeneID, rather than symbols or genomic coordinates. A GeneID found in only one annotation is not necessarily a biological gene discovery or deletion: annotation splits, merges, replacement IDs, and filtering can contribute. Renaming counts compare the GTF gene_id labels. Location changes compare the full sets of gene intervals on the same assembly, and include boundary edits. These are annotation counts, not numbers of expressed or QTL-tested genes.

Slide wording: “Updated RefSeq annotation on the same GRCr8 assembly: 49,017 genes (+4.3%) and 138,231 transcripts (+19.1%). Protein-coding gene counts are nearly unchanged (+0.2%), while lncRNA coverage grows by 6.4%.”
