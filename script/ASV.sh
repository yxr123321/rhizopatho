Rscript dada2.r
sed -i 's/"//g' table.rechim.xls
perl dada2normal.pl table.rechim.xls
sample_order.pl asv_table.tmp.xls raw.fq.list asv_table.xls
usearch11 -closed_ref asv_rep.fasta -db silva_MBPD.fasta -strand both -otutabout close/otu_rep.txt -tabbedout close/closed.txt