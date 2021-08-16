# Introduction
**millerv2/NanoseqSnakemake** is a bioinformatics analysis pipeline that can be used to perform basecalling, demultiplexing, QC, mapping and downstream analysis of Nanopore long-read RNA sequencing data.

The pipeline is written using Snakemake, a Python-baased workflow tool that allows tasks to be run across multiple compute infrastructures in a very portable manner.

# Pipeline Overview
1. Sequencing QC ([`NanoPlot`](https://github.com/wdecoster/NanoPlot))
2. Raw read QC (NanoPlot, FastQC)
3. Read Filtering (Nanofilt)
4. Alignment (minimap2)
    * Capable of performing unspliced and spliced alignment. Sensible defaults will be applied automatically, but specific alignment parameters will be run based on a combination of the input data type and user-specified parameters.
    * Convert SAM to co-ordinate sorted BAM and obtain mapping metrics (SAMtools)
5. Create bigWig (BEDTools, bedGraphToBigWig) and bigBed (BEDTools, bedToBigBed) coverage tracks for visualisation
6. RNA-specific downstream analysis:
    * Transcript reconstruction and quantification (StringTie2)
    * Each sample can be processed individually and combined using Stringtie Merge. After which, featureCounts will be used for both gene and transcript      quantification.
    * Differential expression analysis (DESeq2)
    * Principal component analysis of samples and heatmap generation (pheatmap)
7. Present QC for raw read and alignment results (MultiQC)

This workflow is based on and attempts to improve upon the Nextflow-based Nanoseq Bioinformatics pipeline described here (https://github.com/nf-core/nanoseq/blob/master/README.md) 

