#! /bin/bash
# reformat star count data for DESeq2


for fin in ../0_raw_data/star_counts/*.out.tab; do
	fout=${fin##*/}
	# we will use a different name to keep things distinct
	fout=${fout/.tab/.fixed.tab}
	# we only need lines 5 and onward (tail -n +5)
	# and we only want columns 1 and 2 (cut -f 1,2)
	tail -n +5 $fin | cut -f 1,2 > $fout
done
