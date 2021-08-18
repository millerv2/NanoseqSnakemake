# Introduction
**millerv2/NanoseqSnakemake** is a bioinformatics analysis pipeline that can be used to perform QC, mapping and downstream analysis of Nanopore long-read RNA sequencing data.

The pipeline is written using Snakemake, a Python-based workflow tool that allows tasks to be run across multiple compute infrastructures in a very portable manner.

# Pipeline Overview
1. Raw read QC ([`NanoPlot`](https://github.com/wdecoster/NanoPlot), [`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Read Filtering ([`Nanofilt`](https://github.com/wdecoster/nanofilt))
3. Alignment ([`minimap2`](https://github.com/lh3/minimap2))
    * Capable of performing unspliced and spliced alignment. Sensible defaults will be applied automatically, but specific alignment parameters will be run based on a combination of the input data type and user-specified parameters.
    * Convert SAM to co-ordinate sorted BAM and obtain mapping metrics (SAMtools)
4. Create bigWig ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/)) and bigBed ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedToBigBed`](http://hgdownload.soe.ucsc.edu/admin/exe/)) coverage tracks for visualisation
5. RNA-specific downstream analysis:
    * Transcript reconstruction and quantification ([`StringTie2`](https://ccb.jhu.edu/software/stringtie/))
    * Each sample can be processed individually and combined using Stringtie Merge. After which, featureCounts will be used for both gene and transcript      quantification.
    * Differential expression analysis ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
    * Principal component analysis of samples and heatmap generation ([`pheatmap`](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap/))
6. Present QC for raw read and alignment results ([`MultiQC`](https://multiqc.info/docs/))

This workflow is based on and attempts to improve upon the Nextflow-based Nanoseq Bioinformatics pipeline described here (https://github.com/nf-core/nanoseq/blob/master/README.md) 

# Running the Workflow
