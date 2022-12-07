# make sure to set the working directory to the location of this file in the git repo
# in RStudio, use the menu at the top of the screen
# Session->Set Working Directory->To Source File Location

# install required libraries if necessary
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", update=FALSE)

# load the libraries
library("DESeq2")

# gather the count data files
files = grep("out.fixed.tab", list.files(path="."), value=TRUE)

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

# import the data into a DESeq Data Set
# so DESeq2 can give us information on our multiple variables, we define design=~cultivar+coldHours
dds = DESeqDataSetFromHTSeqCount(sampleTable=df, directory=".", design=~cultivar+coldHours)

# filter lowly expressed reads
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

# perform the differential expression analysis
dds = DESeq(dds)

# produce a dispersion plot in order to determine if the data is a good match for the DESeq2 model
# for more info, see: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
# the plot should indicate that this data is a good match for the DESeq2 model
# see visualizations/DESeq2_star_dispersion_plot.png in this git repo to view without running this script
plotDispEsts(dds)

# collect the results for each of the categories we tested
# contrast arg determines what relationships to examine within the data
# alpha=0.05 sets the adjusted p-value cutoff to p=0.05
resCultivar = results(dds, contrast=c("cultivar", "Bergeron", "Badami"), alpha=0.05)
summary(resCultivar)
# out of 22538 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3601, 16%
# LFC < 0 (down)     : 2770, 12%
# outliers [1]       : 159, 0.71%
# low counts [2]     : 874, 3.9%
# (mean count < 2)

resColdHours = results(dds, contrast=c("coldHours", "400", "800"), alpha=0.05)
summary(resColdHours)
# out of 22538 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3348, 15%
# LFC < 0 (down)     : 4224, 19%
# outliers [1]       : 159, 0.71%
# low counts [2]     : 874, 3.9%
# (mean count < 2)
