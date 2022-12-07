# make sure to set the working directory to the location of this file in the git repo
# in RStudio, use the menu at the top of the screen
# Session->Set Working Directory->To Source File Location

# install required libraries if necessary
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", update=FALSE)
BiocManager::install("GenomicFeatures", update=FALSE)
BiocManager::install("tximport", update=FALSE)

# load the libraries
library("DESeq2")
library("GenomicFeatures")
library("tximport")

# load gene data from the gff3
GFF = "../0_raw_data/Ppersica_298_v2.1.gene_exons.gff3"
txdb = makeTxDbFromGFF(GFF, "gff3", "phytosome", "Prunus persica")
k = keys(txdb, keytype="TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")

# gather the quantity data files
files = grep("*.sf", list.files(path="."), value=TRUE)

# import the quantity data
txi = tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM")

# extract slightly more readable sample name designations from the file names
names = sub("(.*Clone[0-9]+).*", '\\1', files)

# extract the cultivar names from the file names
# the data is cast into type factor for downstream analysis
cultivars = factor(sub("(^[A-Za-z]+).*", "\\1", files))

# extract the cold hour information for each sample
# again, the data is cast into type factor for downstream analysis
coldHours = factor(sub(".*([0-9]{3})CH.*", "\\1", files))

# create a dataframe from the data
df = data.frame(name=names, file=files, cultivar=cultivars, coldHours=coldHours)

# perform the differential expression analysis
dds = DESeqDataSetFromTximport(txi, df, ~cultivar+coldHours)

# filter lowly expressed reads
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

# perform the differential expression analysis
dds = DESeq(dds)

# produce a dispersion plot in order to determine if the data is a good match for the DESeq2 model
# for more info, see: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
# the plot should indicate that this data is a good match for the DESeq2 model
# see visualizations/DESeq2_salmon_dispersion_plot.png in this git repo to view without running this script
plotDispEsts(dds)

# collect the results for each of the categories we tested
# contrast arg determines what relationships to examine within the data
# alpha=0.05 sets the adjusted p-value cutoff to p=0.05
resCultivar = results(dds, contrast=c("cultivar", "Bergeron", "Badami"), alpha=0.05)
summary(resCultivar)
# out of 22033 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3318, 15%
# LFC < 0 (down)     : 2546, 12%
# outliers [1]       : 155, 0.7%
# low counts [2]     : 855, 3.9%
# (mean count < 2)

resColdHours = results(dds, contrast=c("coldHours", "400", "800"), alpha=0.05)
summary(resColdHours)
# out of 22033 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3076, 14%
# LFC < 0 (down)     : 3986, 18%
# outliers [1]       : 155, 0.7%
# low counts [2]     : 855, 3.9%
# (mean count < 2)
