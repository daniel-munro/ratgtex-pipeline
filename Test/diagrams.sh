# Run this from the main directory ('sh Test/diagrams.sh')
snakemake --rulegraph Test/Test.aFC.txt Test/Test.cis_qtl_signif.txt.gz Test/Test.cis_qtl_all_pvals.txt.gz Test/Test.trans_qtl_pairs.txt.gz | dot -Tpng > Test/rulegraph.png
snakemake --filegraph Test/Test.aFC.txt Test/Test.cis_qtl_signif.txt.gz Test/Test.cis_qtl_all_pvals.txt.gz Test/Test.trans_qtl_pairs.txt.gz | dot -Tpng > Test/filegraph.png
