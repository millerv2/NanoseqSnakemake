#!/bin/bash
WORKDIR=$1
# runmode = cluster or dryrun or unlock
RUNMODE=$2

PIPELINE_HOME="/data/millerv2/NanoseqSnakemake"
SNAKEFILE="$PIPELINE_HOME/workflow/Snakefile"

if [ "$RUNMODE" == "dryrun" ];then
	snakemake -n \
	-s $SNAKEFILE \
	--directory $WORKDIR \
	--printshellcmds \
	--cluster-config $PIPELINE_HOME/resources/cluster.json \
	--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} " \
	-j 500
elif [ "$RUNMODE" == "unlock" ];then
	snakemake --unlock \
	-s $SNAKEFILE \
	--directory $WORKDIR
elif [ "$RUNMODE" == "cluster" ];then
	if [ ! -d $WORKDIR ]; then mkdir -p $WORKDIR;fi
	cat << EOF > submit.sh 
	#!/bin/bash
	#SBATCH --cpus-per-task=4 
	#SBATCH --mem=8g
	#SBATCH --time=1-00:00:00 
	#SBATCH --parsable 
	#SBATCH -J "Nanoseq" 
	#SBATCH --mail-type=BEGIN,END,FAIL
	module load snakemake
	module load singularity
	snakemake \
	-s $SNAKEFILE \
	--directory $WORKDIR \
	--use-singularity \
	--singularity-args "'-B $WORKDIR'" \
	--use-envmodules \
	--printshellcmds \
	--latency-wait 120 \
	--cluster-config $PIPELINE_HOME/resources/cluster.json \
	--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} " \
	-j 500 \
	--rerun-incomplete \
	--keep-going \
	--stats "${WORKDIR}/snakemake.stats" \
	2>&1 | tee "${WORKDIR}/snakemake.log"
	EOF

	sbatch submit.sh
fi