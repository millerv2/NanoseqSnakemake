rule DESeq:
    input:
        script_dir = os.path.join(base_dir,"NanoseqSnakemake","workflow","scripts"),
        genecounts = os.path.join(base_dir,"featureCounts","counts_gene.txt"),
        samples = "/data/millerv2/NanoseqSnakemake/config/samples.tsv"
        #samples = config['samples']
    output:
        os.path.join(base_dir,"DESeq2","deseq2.results.txt")
    message:
        "Running DESeq2 with {input}"
    envmodules:
        "R/4.1.0"
    shell:"""
    cd {input.script_dir}
    ./run_DESeq2.R {input.genecounts} {input.samples}
    """