rule DESeq:
    input:
        genecounts = os.path.join(base_dir,"featureCounts","counts_gene.txt"),
        samples = config['samples']
    output:
        os.path.join(base_dir,"DESeq2","deseq2.results.txt")
    params:
        script_dir = os.path.join(base_dir,"scripts"),
        script_name = "run_DESeq2.R"
    message:
        "Running DESeq2 with {input}"
    envmodules:
        "R/4.1.0"
    shell:"""
    Rscript {params.script_dir}/{params.script_name} {input.genecounts} {input.samples}
    """
