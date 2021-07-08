# Short workflow to run fastqc and cutadapt.
SAMPLES = ["SGNex_A549_directRNA_replicate1_run1","SGNex_A549_directRNA_replicate6_run1"]

#orkflow_base = '/data/millerv2'

rule all:
    input:
        expand("/data/millerv2/fastqc_raw/{sample}_fastqc.html",sample=SAMPLES),
        expand("/data/millerv2/fastqc_raw/{sample}_fastqc.zip",sample=SAMPLES),
        "/data/millerv2/Nanoplot/Yield_By_Length.png",
        "/data/millerv2/Nanoplot/Yield_by_Length.html"
        #expand("trimmed_reads/{sample}_trimmed.fastq.gz",sample=SAMPLES),
        #expand("trimmed_reads/{sample}_trimmed.fastq.gz",sample=SAMPLES),
        #expand("fastqc_trimmed/{sample}_trimmed_fastqc.html",sample=SAMPLES),
        #expand("fastqc_trimmed/{sample}_trimmed_fastqc.zip",sample=SAMPLES)

rule fastqc_raw:
    input:
        "/data/millerv2/samples/{sample}.fastq.gz"
    envmodules:
        "fastqc/0.11.9"
    output:
        "/data/millerv2/fastqc_raw/{sample}_fastqc.html",
        "/data/millerv2/fastqc_raw/{sample}_fastqc.zip"
    shell:'''
    fastqc -o /data/millerv2/fastqc_raw {input}
    '''



rule nanoplot:
    container: 
        "docker://continuumio/miniconda3:4.4.10"
    input:
        expand("/data/millerv2/samples/{sample}.fastq.gz",sample=SAMPLES)
        #"/data/millerv2/samples/{sample}.fastq.gz"
        #"/data/millerv2/samples/SGNex_A549_directRNA_replicate1_run1.fastq.gz",
        #"/data/millerv2/samples/SGNex_A549_directRNA_replicate6_run1.fastq.gz"
    conda:
        "envs/nanoplot.yaml"
# This container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
    output:
        "/data/millerv2/Nanoplot/Yield_By_Length.png",
        "/data/millerv2/Nanoplot/Yield_by_Length.html"
    shell:'''
    NanoPlot -o /data/millerv2/Nanoplot/summary-plots-log-transformed 
    --fastq {input} /
    '''

#
""" rule trim_reads:
    input:
        "A549samples/{sample}.fastq.gz"
    envmodules:
        "cutadapt/3.4"
    output:
        "trimmed_reads/{sample}_trimmed.fastq.gz"
    shell:'''
    cutadapt -q 15,10 -o {output} {input}
    '''
 """
#testing
""" rule fastqc_trimmed:
    input:
        "trimmed_reads/{sample}_trimmed.fastq.gz"
    envmodules:
        "fastqc/0.11.9"
    output:
        "fastqc_trimmed/{sample}_trimmed_fastqc.html",
        "fastqc_trimmed/{sample}_trimmed_fastqc.zip"
    shell:'''
    fastqc -o fastqc_trimmed {input}
    ''' """


 