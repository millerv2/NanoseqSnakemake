rule bedtools_bam_to_bed:
    input:
        os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam")
    output:
        os.path.join(base_dir,"bedtools","bed_files","{sample}.bed")
    message:
        "Running bedtools with {input}"
    envmodules:
        "bedtools/2.30.0"
    shell:'''
    bedtools bamtobed -i {input} > {output}
    '''

#convert BED to BigBed track for visualization 

rule bed_to_bigbed:
    input:
        bed = os.path.join(base_dir,"bedtools","bed_files","{sample}.bed"),
        chrom_sizes = os.path.join(base_dir,"genomes","chrom.sizes")
    output:
        os.path.join(base_dir,"bedtools","bed_files","{sample}.bigBed")
    message:
        "Running bedtools with {input}"
    container:
        "docker://zavolab/bedtobigbed:2.7"
    shell:'''
    bedtoBigBed {input.bed} {input.chrom_sizes} > {output}
    '''