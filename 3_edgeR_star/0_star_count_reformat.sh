#! /bin/bash
# reformat star count data for edgeR

# setup a temporary staging folder for intermediate operations
mkdir -p temp

# first, extract all of the data from each raw STAR count
for fin in ../0_raw_data/star_counts/*.out.tab; do
	fout=${fin##*/}
	fout=temp/${fout/.tab/.fixed.tab}
	# we only need lines 5 and onward (tail -n +5)
	# and we only want columns 1 and 2 (cut -f 1,2)
	tail -n +5 $fin | cut -f 1,2 > $fout
done

# now we will join all the data into one table for edgeR

files=($(ls ./temp/*.fixed.tab))
tmp=temp/join.tmp

# use join to merge all files into one table
join -a1 -a2 -e 0 -o auto <(sort ${files[0]}) <(sort ${files[1]}) > $tmp
rest=(${files[@]/${files[0]}})
rest=(${rest[@]/${files[1]}})
for file in ${rest[@]}; do
	join -a1 -a2 -e 0 -o auto $tmp <(sort $file) > $tmp.1
	mv $tmp.1 $tmp
done

# convert all spaces to tabs in the temp file
cat $tmp | tr " " "\t" > $tmp.1
mv $tmp.1 $tmp

# create the final output file with a header
# create the header
echo -n "gene" > gene_counts.tsv
for i in ${!files[@]}; do
	name=${files[$i]##*/}
	name=${name/ReadsPerGene*}
	echo -e -n "\t${name}" >> gene_counts.tsv
done
echo "" >> gene_counts.tsv

# append the gene counts to the file
cat $tmp >> gene_counts.tsv

# clean up intermediate files
rm -rf temp
