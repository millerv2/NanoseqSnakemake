rule multiQC:
    input:
        base_dir
    output:
        os.path.join(base_dir,"multiqc","multiqc_report.html")
    message:
        "Running multiqc with recursion starting in directory: {input}"
    envmodules:
        "multiqc/1.11"
    shell:'''
    multiqc {input}
    '''