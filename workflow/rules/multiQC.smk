rule multiQC:
    input:
        expand(os.path.join(base_dir,"fastqc_raw","{sample}_fastqc.html"),sample=SAMPLES),
        expand(os.path.join(base_dir,"fastqc_raw","{sample}_fastqc.zip"),sample=SAMPLES),
        expand(os.path.join(base_dir,"trimmed_reads","{sample}_trimmed.fastq.gz"),sample=SAMPLES),
        os.path.join(base_dir,"Nanoplot","NanoPlot-report.html"),
        expand(os.path.join(base_dir,"fastqc_trimmed","{sample}_trimmed_fastqc.html"),sample=SAMPLES),
        expand(os.path.join(base_dir,"fastqc_trimmed","{sample}_trimmed_fastqc.zip"),sample=SAMPLES),
        expand(os.path.join(base_dir,"genome_alignments","{sample}_alignment.sam"),sample=SAMPLES),
        expand(os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam"),sample=SAMPLES),
        expand(os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam.flagstat"),sample=SAMPLES),
        expand(os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam.stats"),sample=SAMPLES),
        expand(os.path.join(base_dir,"bedtools","bed_files","{sample}.bed"),sample=SAMPLES),
        os.path.join(base_dir,"genomes","chrom.sizes"),
        expand(os.path.join(base_dir,"bedtools","big_bed_files","{sample}.bigBed"),sample=SAMPLES),
        expand(os.path.join(base_dir,"bedtools","bedgraph_files","{sample}.bedGraph"),sample=SAMPLES),
        expand(os.path.join(base_dir,"bedtools","bigWig_tracks","{sample}.bigWig"),sample=SAMPLES),
        expand(os.path.join(base_dir,"stringtie","assemblies","{sample}.stringtie.gtf"),sample=SAMPLES),
        os.path.join(base_dir,"stringtie","mergelist","mergelist.txt"),
        os.path.join(base_dir,"stringtie","merged","stringtie.merged.gtf"),
        expand(os.path.join(base_dir,"stringtie","assemblies","{sample}.stringtie_final.gtf"),sample=SAMPLES),
        os.path.join(base_dir,"featureCounts","counts_gene.txt"),
        os.path.join(base_dir,"featureCounts","counts_transcript.txt")
    output:
        os.path.join(base_dir,"multiqc","multiqc_report.html")
    message:
        "Running multiqc starting in directory: {}".format(base_dir)
    envmodules:
        "multiqc/1.11"
    params:
        workdir = base_dir,
        outdir = os.path.join(base_dir,"multiqc")
    shell:'''
    multiqc --ignore '*/.singularity/*' -f --interactive --outdir {params.outdir} {params.workdir}
    '''