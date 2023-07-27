# Run this from the main directory ('sh test/diagrams.sh')

cp -r test/Test rn7/

snakemake --rulegraph \
    rn7/Test/Test.aFC.txt \
    rn7/Test/Test.cis_qtl_signif.txt.gz \
    rn7/Test/Test.cis_qtl_all_pvals.txt.gz \
    rn7/Test/Test.trans_qtl_pairs.txt.gz \
    rn7/Test/qc/rna_to_geno_summary.tsv \
    rn7/Test/qc/Test.sex_concordance.txt \
    rn7/Test/splice/Test_splice.cis_independent_qtl.txt.gz \
    rn7/Test/splice/Test_splice.cis_qtl_signif.txt.gz \
    rn7/Test/splice/Test_splice.trans_qtl_pairs.txt.gz \
    | dot -Tpng \
    > test/rulegraph.png

snakemake --filegraph \
    rn7/Test/Test.aFC.txt \
    rn7/Test/Test.cis_qtl_signif.txt.gz \
    rn7/Test/Test.cis_qtl_all_pvals.txt.gz \
    rn7/Test/Test.trans_qtl_pairs.txt.gz \
    rn7/Test/qc/rna_to_geno_summary.tsv \
    rn7/Test/qc/Test.sex_concordance.txt \
    rn7/Test/splice/Test_splice.cis_independent_qtl.txt.gz \
    rn7/Test/splice/Test_splice.cis_qtl_signif.txt.gz \
    rn7/Test/splice/Test_splice.trans_qtl_pairs.txt.gz \
    | dot -Tpng \
    > test/filegraph.png

rm -r rn7/Test
