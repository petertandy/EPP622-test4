#!/bin/bash

for folder in ../0_raw_data/salmon_quants/*; do
	sample_name=${folder##*/}
	sample_name=${sample_name/_quant}
	cp ${folder}/quant.sf ${sample_name}.sf
done