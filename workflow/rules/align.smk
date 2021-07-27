rule minimap_alignment:
    input:
        ref = config['genome_fasta'],
        fastq = os.path.join(base_dir,"samples","{sample}.fastq")
    output:
        os.path.join(base_dir,"genome_alignmments","{sample}_alignment.sam")
    log:
        os.path.join(base_dir,"logs","{sample}_minimap2.log")
    message:
        "Running minimap with {input} with {threads} threads"
    envmodules:
        "minimap2/2.20"
    threads:
        12
    shell:'''
    minimap2 -t {threads} -ax splice -uf -k14 {input.ref} {input.fastq} > {output}
    '''