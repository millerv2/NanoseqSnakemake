# Short workflow to run fastqc and cutadapt.
SAMPLES = ["SGNex_A549_directRNA_replicate1_run1"]
# SGNex_A549_directRNA_replicate5_run1, SGNex_A549_directRNA_replicate6_run1

rule all:
    input:
        expand("fastqc_raw/{sample}_fastqc.html",sample=SAMPLES),
        expand("fastqc_raw/{sample}_fastqc.zip",sample=SAMPLES)
        #expand("trimmed_reads/{sample}_trimmed.fastq.gz",sample=SAMPLES),
        #expand("trimmed_reads/{sample}_trimmed.fastq.gz",sample=SAMPLES),
        #expand("fastqc_trimmed/{sample}_trimmed_fastqc.html",sample=SAMPLES),
        #expand("fastqc_trimmed/{sample}_trimmed_fastqc.zip",sample=SAMPLES)

rule fastqc_raw:
    input:
        "samples/{sample}.fastq.gz"
    envmodules:
        "fastqc/0.11.9"
    output:
        "fastqc_raw/{sample}_fastqc.html",
        "fastqc_raw/{sample}_fastqc.zip"
    shell:'''
    fastqc -o fastqc_raw {input}
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

""" rule nanoplot:
    input:
        "A549samples/{sample}.fastq.gz"
    conda:
        "envs/nanoplot.yaml"
    output:
        "Yield_By_Length.png",
        "Yield_by_Length.html"
    shell:'''
    Nanoplot -o summary-plots-log-transformed 
    --fastq {sample}.fastq.gz /
    '''
 """