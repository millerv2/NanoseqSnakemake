rule fastqc_raw:
    input:
        os.path.join(base_dir,"samples","{sample}.fastq")
    output:
        os.path.join(base_dir,"fastqc_raw","{sample}_fastqc.html"),
        os.path.join(base_dir,"fastqc_raw","{sample}_fastqc.zip")
    params:
        OUT_DIR = os.path.join(base_dir,"fastqc_raw")
    log:
        os.path.join(base_dir,"logs","{sample}_fastqc_raw.log")
    message:
        "Running fastqc_raw with {input}"
    envmodules:
        "fastqc/0.11.9"
    shell:
        "fastqc -o {params.OUT_DIR} {input} &> {log} "
    #test

rule nanoplot:
    input:
        expand(os.path.join(base_dir,"samples","{sample}.fastq"),sample=SAMPLES)
    container:
        "docker://staphb/nanoplot:1.33.0"
# This container defines the underlying OS for each job when using the workflow
# with --use-singularity
    output:
        os.path.join(base_dir,"Nanoplot","Dynamic_Histogram_Read_length.html"),
        os.path.join(base_dir,"Nanoplot","Dynamic_Histogram_Read_length.png"),
        os.path.join(base_dir,"Nanoplot","HistogramReadlength.png"),
        os.path.join(base_dir,"Nanoplot","LengthvsQualityScatterPlot_dot.png"),
        os.path.join(base_dir,"Nanoplot","LengthvsQualityScatterPlot_kde.png"),
        os.path.join(base_dir,"Nanoplot","LogTransformed_HistogramReadlength.png"),
        os.path.join(base_dir,"Nanoplot","NanoPlot-report.html"),
        os.path.join(base_dir,"Nanoplot","NanoStats.txt"),
        os.path.join(base_dir,"Nanoplot","Weighted_HistogramReadlength.png"),
        os.path.join(base_dir,"Nanoplot","Weighted_LogTransformed_HistogramReadlength.png"),
        os.path.join(base_dir,"Nanoplot","Yield_By_Length.png")
    params:
        OUT_DIR = os.path.join(base_dir,"Nanoplot")
    log:
        os.path.join(base_dir,"logs","Nanoplot.log")
    message:
        "Running Nanoplot with {input}"
    shell:'''
    NanoPlot -o {params.OUT_DIR} --fastq {input} 
    '''
