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
        sorted_bed = os.path.join(base_dir,"bedtools","bed_files","{sample}.sorted.bed"),
        big_bed = os.path.join(base_dir,"bedtools","big_bed_files","{sample}.bigBed")
    message:
        "Running bedtools with {input}"
    envmodules:
        "ucsc/416"
    shell:'''
    sort -k1,1 -k2,2n {input.bed} > {output.sorted_bed}
    bedToBigBed {output.sorted_bed} {input.chrom_sizes} {output.big_bed}
    '''