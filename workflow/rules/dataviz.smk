rule dataviz:
    input:
        genecounts = os.path.join(base_dir,"featureCounts","counts_gene.txt"),
        samples = config['samples']
        #samples = config['samples']
    output:
        os.path.join(base_dir,"downstream","sampleheatmap.png"),
        os.path.join(base_dir,"downstream","PCAsamples.png"),
        os.path.join(base_dir,"downstream","top50heatmap.png")
    message:
        "Running dataviz script with {input}"
    params:
        out_dir = os.path.join(base_dir,"downstream"),
        scripts = os.path.join(base_dir,"NanoseqSnakemake","workflow","scripts","PCA_and_heatmaps.R")
    envmodules:
        "R/4.1.0"
    shell:"""
    mkdir -p {params.out_dir}
    {params.scripts} {input.genecounts} {input.samples} {params.out_dir}
    """