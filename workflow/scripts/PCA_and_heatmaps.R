#!/usr/bin/env Rscript

#install.packages("pheatmap")
#install.packages("RColorBrewer")
library(pheatmap)
library(RColorBrewer)
library(DESeq2)

#set path to counts matrix and sample table:

args = commandArgs(trailingOnly=TRUE)
#pathtocounts <- "counts_gene.txt"
pathtocounts <- args[1]
pathtosamples <- args[2]
outputdir <- args[3]
#outputdir <- '/Users/millerv2/Desktop/'
#outputdir <- "/data/millerv2/data_viz/"
#pathtosamples <- "samples.tsv"
setwd(outputdir)
sampleheatmap <- 'sampleheatmap.png'
pcaplot <- "PCAsamples.png"
top50heatmap <- "top50heatmap.png"

#read in counts matrix and sample data table to make a DESeqDataSet
count.matrix <- read.table(pathtocounts,sep="\t",header=TRUE, skip = 1,check.names = FALSE)
count.matrix$Chr   <- count.matrix$Start <- count.matrix$End <- count.matrix$Length <- count.matrix$Strand <- NULL
colnames(count.matrix) <- gsub("(.sorted){0,1}.bam$","",basename(colnames(count.matrix)))
count.matrix <- aggregate(count.matrix[,-1],count.matrix["Geneid"],sum)
countTab           <- count.matrix[,-1]
rownames(countTab) <-count.matrix[,1]


################################################
################################################
## READ IN SAMPLE INFORMATION (CONDITIONS)    ##
################################################
################################################
sample_sheet <- pathtosamples
sample_sheet <- data.frame(read.table(sample_sheet,sep="",fill=FALSE,header=TRUE))
sample <- colnames(countTab)
group <- sample_sheet$group
app <- sample_sheet$application
sampInfo <- data.frame(group, row.names = sample)
if (!all(rownames(sampInfo) == colnames(countTab))){
  sampInfo <- sampInfo[match(colnames(countTab), rownames(sampInfo)),]
}

dds <- DESeqDataSetFromMatrix(countData = countTab, colData = sample_sheet, design = ~ group)
#rlog transformation of dagtaset
rld <- rlog(dds, blind=FALSE)

#make a sample distance heatmap:
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$application, rld$group, sep="-")
colnames(sampleDistMatrix) <- paste(rld$application, rld$group, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(sampleheatmap)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
dev.off()


#make a PCA plot of the samples
png(pcaplot)
plotPCA(rld, intgroup = c("application", "group","replicate"))
dev.off()

#heatmap of 50 most variable genes
geneVars <- rowVars(assay(rld))
geneVarsOrdered <- order(geneVars, decreasing = TRUE)
topVarGenes <- head(geneVarsOrdered, 50)
mat <- assay(rld)[ topVarGenes, ]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("application","group")])
clear_col_names <- paste( rld$group, rld$application, sep=".")
png(top50heatmap,width=750,height=700,units="px")
topGenesHeatmap <- pheatmap(mat, annotation_col=df, labels_col = clear_col_names)
dev.off()