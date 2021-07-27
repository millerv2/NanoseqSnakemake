rule nanofilt:
    input:
        os.path.join(base_dir,"samples","{sample}.fastq") 
    output:
        os.path.join(base_dir,"trimmed_reads","{sample}_trimmed.fastq") 
    container:
        "docker://mcfonsecalab/nanofilt:latest"
    params:
        OUT_DIR = os.path.join(base_dir,"trimmed_reads"),
        headcrop = config['headcrop'],
        length = config['length']
    log:
        os.path.join(base_dir,"logs","{sample}_Nanofilt.log")
    message:
        "Running Nanofilt with {input}"
    shell:'''
    NanoFilt -l {params.length} --headcrop {params.headcrop} < {input} > {output}
    '''  