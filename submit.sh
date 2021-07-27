#!/bin/bash
#SBATCH --cpus-per-task=4 
#SBATCH --mem=8g
#SBATCH --time=1-00:00:00 
#SBATCH --parsable 
#SBATCH -J "Nanoseq" 
#SBATCH --mail-type=BEGIN,END,FAIL
module load snakemake
module load singularity
snakemake --directory "/data/millerv2" --use-singularity --singularity-args "'-B /data/millerv2'" --use-envmodules --printshellcmds --latency-wait 120 --cluster-config resources/cluster.json --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} " -j 500 --rerun-incomplete --keep-going --stats "/data/millerv2/snakemake.stats" 2>&1 | tee "/data/millerv2/snakemake.log"
