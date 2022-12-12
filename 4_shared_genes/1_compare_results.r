# make sure to set the working directory to the location of this file in the git repo
# in RStudio, use the menu at the top of the screen
# Session->Set Working Directory->To Source File Location


# load all of the data
deseq2StarCold = read.csv("DESeq2_star_genes_cold.csv")
deseq2StarCultivar = read.csv("DESeq2_star_genes_cultivar.csv")
deseq2SalmonCold = read.csv("DESeq2_salmon_genes_cold.csv")
deseq2SalmonCultivar = read.csv("DESeq2_salmon_genes_cultivar.csv")
edgeRCold = read.csv("edgeR_star_genes_cold.csv")
edgeRCultivar = read.csv("edgeR_star_genes_cultivar.csv")


# the star gene names have an extraneous ".v2.1" that needs to be removed to easily compare them to the salmon genes
# we will keep a copy of the original names versions of the dataframes for comparing star to star
originalStarCold = deseq2StarCold
originalStarCultivar = deseq2StarCultivar
for (i in 1:nrow(deseq2StarCold)) {
  oldName = deseq2StarCold$X[i]
  newName = unlist(strsplit(oldName, ".v"))[1]
  deseq2StarCold[i, "X"] = newName
}
for (i in 1:nrow(deseq2StarCultivar)) {
  oldName = deseq2StarCultivar$X[i]
  newName = unlist(strsplit(oldName, ".v"))[1]
  deseq2StarCultivar[i, "X"] = newName
}


# comparing DESeq2 analyses on star & salmon data
# coldHours
# shared upregulated
deseq2StarColdUp = deseq2StarCold[which(deseq2StarCold$log2FoldChange > 0),]
deseq2SalmonColdUp = deseq2SalmonCold[which(deseq2SalmonCold$log2FoldChange > 0),]
sharedDeseq2ColdUp = intersect(deseq2StarColdUp[,"X"], deseq2SalmonColdUp[,"X"])
# shared downregulated
deseq2StarColdDown = deseq2StarCold[which(deseq2StarCold$log2FoldChange < 0),]
deseq2SalmonColdDown = deseq2SalmonCold[which(deseq2SalmonCold$log2FoldChange < 0),]
sharedDeseq2ColdDown = intersect(deseq2StarColdDown[,"X"], deseq2SalmonColdDown[,"X"])

# cultivar
# shared upregulated
deseq2StarCultivarUp = deseq2StarCultivar[which(deseq2StarCultivar$log2FoldChange > 0),]
deseq2SalmonCultivarUp = deseq2SalmonCultivar[which(deseq2SalmonCultivar$log2FoldChange > 0),]
sharedDeseq2CultivarUp = intersect(deseq2StarCultivarUp[,"X"], deseq2SalmonCultivarUp[,"X"])
# shared downregulated
deseq2StarCultivarDown = deseq2StarCultivar[which(deseq2StarCultivar$log2FoldChange < 0),]
deseq2SalmonCultivarDown = deseq2SalmonCultivar[which(deseq2SalmonCultivar$log2FoldChange < 0),]
sharedDeseq2CultivarDown = intersect(deseq2StarCultivarDown[,"X"], deseq2SalmonCultivarDown[,"X"])

# create a summary table
sharedDeseq2Summary = matrix(c(length(sharedDeseq2ColdUp), length(sharedDeseq2ColdDown), length(sharedDeseq2CultivarUp), length(sharedDeseq2CultivarDown)), ncol=2, byrow=TRUE)
colnames(sharedDeseq2Summary) = c("Upregulated", "Downregulated")
rownames(sharedDeseq2Summary) = c("Cold Hours", "Cultivar")
##
#             Upregulated   Downregulated
# Cold Hours  3675          2854
# Cultivar    2868          2241
##

# save the results to disk
write.csv(sharedDeseq2Summary, "shared_DESeq2_summary.csv")


# comparing star-data-based DESeq2 analysis to star-data-based edgeR analysis
# coldHours
# shared upregulated
originalStarColdUp = originalStarCold[which(originalStarCold$log2FoldChange > 0),]
edgeRColdUp = edgeRCold[which(edgeRCold$logFC > 0),]
sharedStarColdUp = intersect(originalStarColdUp[,"X"], edgeRColdUp[,"X"])
# shared downregulated
originalStarColdDown = originalStarCold[which(originalStarCold$log2FoldChange < 0),]
edgeRColdDown = edgeRCold[which(edgeRCold$logFC < 0),]
sharedStarColdDown = intersect(originalStarColdDown[,"X"], edgeRColdDown[,"X"])

# cultivar
# shared upregulated
originalStarCultivarUp = originalStarCultivar[which(originalStarCultivar$log2FoldChange > 0),]
edgeRCultivarUp = edgeRCultivar[which(edgeRCultivar$logFC > 0),]
sharedStarCultivarUp = intersect(originalStarCultivarUp[,"X"], edgeRCultivarUp[,"X"])
# shared downregulated
originalStarCultivarDown = originalStarCultivar[which(originalStarCultivar$log2FoldChange < 0),]
edgeRCultivarDown = edgeRCultivar[which(edgeRCultivar$logFC < 0),]
sharedStarCultivarDown = intersect(originalStarCultivarDown[,"X"], edgeRCultivarDown[,"X"])

# create a summary table
sharedStarSummary = matrix(c(length(sharedStarColdUp), length(sharedStarColdDown), length(sharedStarCultivarUp), length(sharedStarCultivarDown)), ncol=2, byrow=TRUE)
colnames(sharedStarSummary) = c("Upregulated", "Downregulated")
rownames(sharedStarSummary) = c("Cold Hours", "Cultivar")
##
#             Upregulated   Downregulated
# Cold Hours  3437          2615
# Cultivar    2638          2115
##

# save the results to disk
write.csv(sharedStarSummary, "shared_star_summary.csv")
