rule stringtie_assembly:
    input:
        bam = os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam"),
        ref_annot = config['genome_gtf']
    output:
        os.path.join(base_dir,"stringtie","{sample}","{sample}.stringtie.gtf")
    message:
        "Running stringtie with {input}"
    envmodules:
        "stringtie/2.1.5"
    shell:'''
    stringtie -L -G {input.ref_annot} -o {output} {input.bam}
    '''

rule get_merge_list:
    input:
        expand(os.path.join(base_dir,"stringtie","{sample}","{sample}.stringtie.gtf"),sample=SAMPLES)
    output:
        os.path.join(base_dir,"stringtie","mergelist","mergelist.txt")
    message: 
        "creating list of paths to all individual assembly files (.gtf)"
    shell:'''
    echo {input} > {output}
    '''

rule stringtie_merge:
    input:
        gtfs = os.path.join(base_dir,"stringtie","mergelist","mergelist.txt"),
        ref_annot = config['genome_gtf']
    output:
        os.path.join(base_dir,"stringtie","merged","stringtie.merged.gtf")
    message:
        "Running stringtie with {input}"
    envmodules:
        "stringtie/2.1.5"
    shell:'''
    stringtie --merge {input.gtfs} -G {input.ref_annot} -o {output}
    '''

