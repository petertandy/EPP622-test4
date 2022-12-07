#!/bin/bash

for folder in ../0_raw_data/salmon_quants/*; do
	sample_name=${folder##*/}
	sample_name=${sample_name/_quant}
	cp ${folder}/quant.sf ${sample_name}.sf
done
echo "Copied and renamed the quantity files from 0_raw_data/salmon_quants/(sample name)/quant.sf to 2_DESeq_salmon/(sample name).sf"

# make sure the gff3 file is unzipped for use in analysis
if [ -f ../0_raw_data/Ppersica_298_v2.1.gene_exons.gff3 ]; then
	echo "Not unzipping the gff3 file (../0_raw_data/Ppersica_298_v2.1.gene_exons.gff3.gz), as it appears to already be unzipped."
else
	gzip -dk "../0_raw_data/Ppersica_298_v2.1.gene_exons.gff3.gz"
	echo "Unzipped the gff3 file (../0_raw_data/Ppersica_298_v2.1.gene_exons.gff3.gz)."
fi
