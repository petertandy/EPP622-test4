# make sure to set the working directory to the location of this file in the git repo
# in RStudio, use the menu at the top of the screen
# Session->Set Working Directory->To Source File Location

# install required libraries if necessary
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", update=FALSE)

# load the libraries
library("edgeR")

# load our count data table
countData = read.delim("gene_counts.tsv", row.names="gene")

# create a DGEList from our data
sampleNames = names(countData)
y = DGEList(counts=countData)
# normalize the library sizes
y = calcNormFactors(y)

# create a design containing our groups to test
cultivars = factor(sub("(^[A-Za-z]+).*", '\\1', sampleNames))
coldHours = factor(sub(".*([0-9]{3})CH.*", '\\1', sampleNames))
design = model.matrix(~cultivars+coldHours)

# perform dispersion analysis
y = estimateDisp(y, design=design)
# fit quasi-likelihood glm model
y = glmQLFit(y, design)

# perform the test to determine associations between the counts and genes
lrtCold = glmQLFTest(y, coef=3)

# obtain the results
# n=Inf gets topTags to extract the total that meet the other criteria
# p.value=0.05 includes only genes with p<=0.05
resCold = topTags(lrtCold, n=Inf, p.value=0.05)

# tease out the meaningful categories of data
coldAll = resCold$table
coldUp = coldAll[coldAll$logFC > 0,]
coldDown = coldAll[coldAll$logFC < 0,]

# write these tables to disk for downstream analysis
write.csv(coldAll, file="edgeR_star_genes_cold.csv")

# now the same as above but associating cultivar and gene regulation
lrtCultivar = glmQLFTest(y, coef=2)
resCultivar = topTags(lrtCultivar, n=Inf, p.value=0.05)
cultivarAll = resCultivar$table
cultivarUp = cultivarAll[cultivarAll$logFC > 0,]
cultivarDown = cultivarAll[cultivarAll$logFC < 0,]
write.csv(cultivarAll, file="edgeR_star_genes_cultivar.csv")

# create a summary table
coldCountUp = nrow(coldUp)
coldCountDown = nrow(coldDown)
cultivarCountUp = nrow(cultivarUp)
cultivarCountDown = nrow(cultivarDown)
tab = matrix(c(coldCountUp, coldCountDown, cultivarCountUp, cultivarCountDown), ncol=2, byrow=TRUE)
colnames(tab) = c("Upregulated", "Downregulated")
rownames(tab) = c("Cold Hours", "Cultivar")

# write the table to disk
write.csv(tab, file="edgeR_star_summary.csv")
