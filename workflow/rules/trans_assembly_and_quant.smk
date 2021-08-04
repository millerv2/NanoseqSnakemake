rule stringtie_assembly:
    input:
        bam = os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam"),
        ref_annot = config['genome_gtf']
    output:
        os.path.join(base_dir,"stringtie","assemblies","{sample}.stringtie.gtf")
    message:
        "Running stringtie with {input}"
    envmodules:
        "stringtie/2.1.5"
    shell:'''
    stringtie -L -G {input.ref_annot} -o {output} {input.bam}
    '''

rule get_merge_list:
    input:
        os.path.join(base_dir,"stringtie","assemblies")
    output:
        os.path.join(base_dir,"stringtie","mergelist","mergelist.txt")
    message: 
        "putting list of paths to all individual assembly files (.gtf) into a text file"
    shell:'''
    cd {input}
    find "$(pwd)" -type f > {output}
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

rule stringtie_assembly_final:
    input:
        bam = os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam"),
        merged_gtf = os.path.join(base_dir,"stringtie","merged","stringtie.merged.gtf")
    output:
        os.path.join(base_dir,"stringtie","assemblies","{sample}.stringtie_final.gtf")
    message:
        "Running stringtie_final with {input}"
    envmodules:
        "stringtie/2.1.5"
    shell:'''
    stringtie -L -G {input.merged_gtf} -o {output} {input.bam}
    '''