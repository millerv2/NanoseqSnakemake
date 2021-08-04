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
    shell:'''
    featureCounts -L -O -f -g gene_id -t exon -a {input.merged_gtf} -o {output.genes} {input.bams}
    featureCounts -L -O -f --primary --fraction -F GTF -g transcript_id -t transcript --extraAttributes gene_id -a {input.merged_gtf} -o {output.transcripts} {input.bams}
    '''
 