rule minimap_alignment:
    input:
        ref = config['genome_fasta'],
        fastq = os.path.join(base_dir,"trimmed_reads","{sample}_trimmed.fastq.gz")
    output:
        os.path.join(base_dir,"genome_alignmments","{sample}_alignment.sam")
    message:
        "Running minimap with {input} with {threads} threads"
    envmodules:
        "minimap2/2.20"
    threads:
        32
    shell:'''
    minimap2 -t {threads} -ax splice -uf -k14 {input.ref} {input.fastq} > {output}
    '''

rule sam_to_sorted_bam:
    input:
        os.path.join(base_dir,"genome_alignments","{sample}_alignment.sam")
    output:
        os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam")
    message:
        "Running Samtools with {input}"
    envmodules:
        "samtools/1.13"
    shell:'''
    samtools sort {input} -o {output}
    '''
    
rule samtools_metrics:
    input:
        os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam")
    output:
        flagstat = os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam.flagstat"),
        stats = os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam.stats")
    message:
        "Running Samtools with {input}"
    envmodules:
        "samtools/1.13"
    shell:'''
    samtools flagstat {input} > {output.flagstat}
    samtools stats {input} > {output.stats}
    '''

#get chromosome sizes in preparation for bedtools data tracks
rule samtools_get_chrom_sizes:
    input:
        ref = config['genome_fasta']
    output:
        os.path.join(base_dir,"genomes","chrom.sizes")
    message:
        "Running Samtools_get_chrom_sizes with {input}"
    envmodules:
        "samtools/1.13"
    shell:'''
    samtools faidx {input.ref}
    cut -f1,2 {input.ref}.fai > {output}
    '''



