rule dataviz:
    input:
        genecounts = os.path.join(base_dir,"featureCounts","counts_gene.txt"),
        samples = config['samples']
        #samples = config['samples']
    output:
        os.path.join(base_dir,"dataviz","sampleheatmap.png"),
        os.path.join(base_dir,"dataviz","PCAsamples.png"),
        os.path.join(base_dir,"dataviz","top50heatmap.png")
    message:
        "Running dataviz script with {input}"
    params:
        out_dir = os.path.join(base_dir,"dataviz"),
        script_dir = os.path.join(base_dir,"scripts"),
        script_name = "PCA_and_heatmaps.R"
    envmodules:
        "R/4.1.0"
    shell:"""
    if [ ! -f {params.out_dir} ]; then mkdir -p {params.out_dir};fi
    Rscript {params.script_dir}/{params.script_name} {input.genecounts} {input.samples} {params.out_dir}
    """