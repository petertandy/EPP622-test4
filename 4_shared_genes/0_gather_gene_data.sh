#! /bin/bash

for d in ../*/; do
	[[ $d == "../4_shared_genes/" ]] && continue
	for f in ${d}*genes*; do
		[[ -f $f ]] && cp $f ./
	done
done
