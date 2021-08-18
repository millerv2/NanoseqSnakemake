rule convert_bam:
    input:
        bam = os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam"),
        chrom_sizes = os.path.join(base_dir,"genomes","chrom.sizes")
    output:
        bed = os.path.join(base_dir,"bedtools","bed_files","{sample}.bed"),
        big_bed = os.path.join(base_dir,"bedtools","big_bed_files","{sample}.bigBed"),
        bedGraph = os.path.join(base_dir,"bedtools","bedgraph_files","{sample}.bedGraph"),
        bigWig = os.path.join(base_dir,"bedtools","bigWig_tracks","{sample}.bigWig")
    message:
        "Running bedtools with {input}"
    envmodules:
        "bedtools/2.30.0","ucsc/416"
    shell:'''
    bedtools bamtobed -i {input.bam} > {output.bed}
    bedSort {output.bed} {output.bed}
    bedToBigBed {output.bed} {input.chrom_sizes} {output.big_bed}
    bedtools genomecov -ibam {input.bam} -bg | bedtools sort > {output.bedGraph}
    bedGraphToBigWig {output.bedGraph} {input.chrom_sizes} {output.bigWig}
'''

