#gene and transcript quantification with featureCounts

rule featureCounts:
    input:
        bams = expand(os.path.join(base_dir,"genome_alignments","sorted_bam","{sample}.sorted.bam"),sample=SAMPLES),
        merged_gtf = os.path.join(base_dir,"stringtie","merged","stringtie.merged.gtf")
    output:
        genes = os.path.join(base_dir,"featureCounts","counts_gene.txt"),
        transcripts = os.path.join(base_dir,"featureCounts","counts_transcript.txt")
    message:
        "Running featurecounts with {input}"
    envmodules:
        "subread/2.0.2"
    params:
        clean_path = os.path.join(base_dir,"genome_alignments","sorted_bam",""),
        extension = ".sorted.bam"
    shell:"""
    featureCounts -L -O -g gene_id -t exon -a {input.merged_gtf} -o {output.genes} {input.bams} 
    featureCounts -L -O -f --primary --fraction -F GTF -g transcript_id -t transcript --extraAttributes gene_id -a {input.merged_gtf} -o {output.transcripts} {input.bams}
    sed -i '2 s@{params.clean_path}@@g' {output.genes} 
    sed -i '2 s@{params.clean_path}@@g' {output.transcripts} 
    sed -i '2 s@{params.extension}@@g' {output.genes}
    sed -i '2 s@{params.extension}@@g' {output.transcripts}
    """
 
 #take and define pattern with s (substitute) and @ is delimiter where in between the pattern is nested. g is global replacement
 #2 is only second line of file, below header line
 #-i is edit file
