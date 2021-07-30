#!/bin/bash
#SBATCH --cpus-per-task=4 
#SBATCH --mem=8g
#SBATCH --time=1-00:00:00 
#SBATCH --parsable 
#SBATCH -J "Nanoseq" 
#SBATCH --mail-type=BEGIN,END,FAIL
set -euo pipefail
module load snakemake
module load singularity
snakemake -s /data/millerv2/NanoseqSnakemake/workflow/Snakefile --directory /data/millerv2 --use-singularity --singularity-args "'-B /data/millerv2,/data/millerv2/NanoseqSnakemake/.tests'" --configfile /data/millerv2/NanoseqSnakemake/config/config.yaml --use-envmodules --printshellcmds --latency-wait 120 --cluster-config /data/millerv2/NanoseqSnakemake/resources/cluster.json --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} " -j 500 --rerun-incomplete --keep-going --stats "/data/millerv2/snakemake.stats" 2>&1 | tee "/data/millerv2/snakemake.log"
