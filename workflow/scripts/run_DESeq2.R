#!/usr/bin/env Rscript

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## DIFFERENTIAL GENE EXPRESSION ANALYSIS
## - GENE COUNT QUANTIFICATION INPUTS FROM EITHER STRINGTIE2+FEATURECOUNTS OR BAMBU
## - SAMPLE INFORMATION INCLUDING CONDITIONS
## - THE PACKAGE BELOW NEEDS TO BE AVAILABLE TO LOAD WHEN RUNNING R

################################################
################################################
## LOAD LIBRARY                               ##
################################################
################################################

library(DESeq2)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
args = commandArgs(trailingOnly=TRUE)
#if (length(args) < 3) {
#    stop("Please featureCounts results, the sample information file, and an output file name", call.=FALSE)
#}
# default output file
#infile <- "counts_gene.txt"
pathtocounts <- args[1]
pathtosamples <- args[2]

#infile <- "counts_gene.txt"
#sample_sheet <- "samples.tsv"
outfile <- "/data/millerv2/DESeq2/deseq2.results.txt"


################################################
################################################
## FORMAT GENE COUNT QUANTIFICATION OUTPUT    ##
################################################
################################################

#create a dataframe for all samples

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
sampInfo <- data.frame(group, row.names = sample)
if (!all(rownames(sampInfo) == colnames(countTab))){
    sampInfo <- sampInfo[match(colnames(countTab), rownames(sampInfo)),]
}

################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################

#make DESeqDataset and write results to file
dds <- DESeqDataSetFromMatrix(countData = countTab, colData = sampInfo, design = ~ group)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=outfile)