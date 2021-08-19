rule minimap_alignment:
    input:
        ref = config['genome_fasta'],
        fastq = os.path.join(base_dir,"filtered_reads","{sample}_filtered.fastq.gz")
    output:
        os.path.join(base_dir,"genome_alignments","{sample}_alignment.sam")
    message:
        "Running minimap_alignment with {input} with {threads} threads"
    envmodules:
        "minimap2/2.20"
    threads: getthreads("minimap_alignment")
    params:
        options = lambda wildcards: config["minimap_options"][sample_to_application[wildcards.sample]]
    shell:'''
    minimap2 -t {threads} {params.options} {input.ref} {input.fastq} > {output}
    '''

rule sam_to_sorted_bam:
    input:
        os.path.join(base_dir,"genome_alignments","{sample}_alignment.sam")
    output:
        os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam")
    message:
        "Running sam_to_sorted_bam with {input}"
    envmodules:
        "samtools/1.13"
    threads: getthreads("sam_to_sorted_bam")
    shell:'''
    samtools sort -@ {threads} -O bam -o {output} {input}
    '''
    
rule samtools_metrics:
    input:
        os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam")
    output:
        flagstat = os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam.flagstat"),
        stats = os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam.stats")
    message:
        "Running samtools_metrics with {input}"
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
        "Running samtools_get_chrom_sizes with {input}"
    envmodules:
        "samtools/1.13"
    shell:'''
cd $(dirname {output})
ln -s {input.ref} .
bn=$(basename {input.ref})
samtools faidx $bn
cut -f1,2 ${{bn}}.fai > {output}
'''



