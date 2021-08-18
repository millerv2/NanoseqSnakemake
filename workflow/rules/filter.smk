rule nanofilt:
    input:
        os.path.join(sample_dir,"{sample}.fastq.gz") 
    output:
        os.path.join(base_dir,"filtered_reads","{sample}_filtered.fastq.gz") 
    container:
        "docker://mcfonsecalab/nanofilt:latest"
    params:
        OUT_DIR = os.path.join(base_dir,"filtered_reads"),
        headcrop = config['nanofilt']['headcrop'],
        length = config['nanofilt']['length']
    log:
        os.path.join(base_dir,"logs","{sample}_Nanofilt.log")
    message:
        "Running Nanofilt with {input}"
    shell:'''
    gunzip -c {input} | NanoFilt -l {params.length} --headcrop {params.headcrop} | gzip - > {output}
    '''  