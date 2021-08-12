rule dataviz:
    input:
        out_dir = os.path.join(base_dir,"downstream"),
        script_dir = os.path.join(base_dir,"NanoseqSnakemake","workflow","scripts"),
        genecounts = os.path.join(base_dir,"featureCounts","counts_gene.txt"),
        samples = config['samples']
        #samples = config['samples']
    output:
        os.path.join(base_dir,"downstream","sampleheatmap.png"),
        os.path.join(base_dir,"downstream","PCAsamples.png"),
        os.path.join(base_dir,"downstream","top50heatmap.png")
    message:
        "Running dataviz script with {input}"
    envmodules:
        "R/4.1.0"
    shell:"""
    cd {input.script_dir}
    ./PCA_and_heatmaps.R {input.genecounts} {input.samples} {input.out_dir}
    """