rule DESeq:
    input:
        os.path.join(base_dir,"featureCounts","counts_gene.txt")
    output:
        os.path.join(base_dir,"DESeq2","deseq2.results.txt")
    message:
        "Running DESeq2 with {input}"
    envmodules:
        "R/4.1.0"
    script:'''
    run_deseq2.r $params.quantification_method {input}
    '''