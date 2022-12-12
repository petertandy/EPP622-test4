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
sampleNames = sub("(.*Clone[0-9]+).*", '\\1', files)

# extract the cultivar names from the file names
# the data is cast into type factor for downstream analysis
cultivars = factor(sub("(^[A-Za-z]+).*", "\\1", files))

# extract the cold hour information for each sample
# again, the data is cast into type factor for downstream analysis
coldHours = factor(sub(".*([0-9]{3})CH.*", "\\1", files))

# create a dataframe from the data
df = data.frame(sampleName=sampleNames, file=files, cultivar=cultivars, coldHours=coldHours)

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
# contrast arg determines what relationship to examine within the data
# the results of this test will show genes upregulated as being upregulated in Begeron vs Badami
resCultivar = results(dds, contrast=c("cultivar", "Bergeron", "Badami"), alpha=0.05)

# first, we'll remove any NAs

# which() filters to obtain only genes that are upregulated (log2FoldChange > 0)
# and that have an adjusted pvalue of < 0.05 (padj < 0.05)
cultivarUp = resCultivar[which(resCultivar$log2FoldChange > 0 & resCultivar$padj < 0.05),]
cultivarDown = resCultivar[which(resCultivar$log2FoldChange < 0 & resCultivar$padj < 0.05),]

# we will collect the counts of these genes separately as well
cultivarUpCount = nrow(cultivarUp) # 3601
cultivarDownCount = nrow(cultivarDown) # 2770

# as a sanity check, make sure summary() agrees with the above counts
summary(resCultivar)
# it does:
##
# out of 22538 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3601, 16%
# LFC < 0 (down)     : 2770, 12%
# outliers [1]       : 159, 0.71%
# low counts [2]     : 874, 3.9%
# (mean count < 2)
##

resCold = results(dds, contrast=c("coldHours", "800", "400"), alpha=0.05)
coldUp = resCold[which(resCold$log2FoldChange > 0 & resCold$padj < 0.05),]
coldDown = resCold[which(resCold$log2FoldChange < 0 & resCold$padj < 0.05),]
coldUpCount = nrow(coldUp) # 4224
coldDownCount = nrow(coldDown) # 3348
summary(resCold)
##
# out of 22538 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 4224, 19%
# LFC < 0 (down)     : 3348, 15%
# outliers [1]       : 159, 0.71%
# low counts [2]     : 874, 3.9%
# (mean count < 2)
##

# write the regulated genes information to disk for use in downstream analysis
# we make sure to select only the rows with p<0.05 for downstream analysis
write.csv(resCultivar[which(resCultivar$padj < 0.05),], "DESeq2_star_genes_cultivar.csv")
write.csv(resCold[which(resCold$padj < 0.05),], "DESeq2_star_genes_cold.csv")

# write a short summary table to disk
tab = matrix(c(coldUpCount, coldDownCount, cultivarUpCount, cultivarDownCount), ncol=2, byrow=TRUE)
colnames(tab) = c("Upregulated", "Downregulated")
rownames(tab) = c("Cold Hours", "Cultivar")
write.csv(tab, "DESeq2_star_summary.csv")
