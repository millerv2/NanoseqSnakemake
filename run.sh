#!/bin/bash
# runmode : dryrun or unlock or cluster or init
RUNMODE=$1
# workdir : dir where all the results will go
WORKDIR=$2
# inputdir : dir where all the input fastqs are located
INPUTDIR=$3

PIPELINE_HOME="$(readlink -f $(dirname "${BASH_SOURCE[0]}"))"
SNAKEFILE="$PIPELINE_HOME/workflow/Snakefile"
CONFIGYAML="$WORKDIR/config.yaml"
CLUSTERJSON="$WORKDIR/cluster.json"

echo "Pipeline Dir:$PIPELINE_HOME"
echo "Snakefile:$SNAKEFILE"
echo "Output Dir:$WORKDIR"
echo "Config file:$CONFIGYAML"
echo "Cluster resources file:$CLUSTERJSON"

module load snakemake
#module load snakemake/6.0.5

function check_workdir_exists() {
if [ ! -d $WORKDIR ];then
	echo "Run \"init\" runmode to create $WORKDIR!"
	exit
fi
}

if [ "$RUNMODE" == "dryrun" ];then

	check_workdir_exists 
	snakemake -n \
	-s $SNAKEFILE \
        --configfile $CONFIGYAML \
	--directory $WORKDIR \
	--printshellcmds \
	--cluster-config $CLUSTERJSON \
	--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} " \
	-j 500 | tee ${WORKDIR}/dryrun.log

elif [ "$RUNMODE" == "init" ];then
	
	if [ -d "$WORKDIR" ];then
		echo "$WORKDIR already exists! Cannot re-initialize!"
		exit
	fi
	mkdir $WORKDIR
	mkdir ${WORKDIR}/logs
	sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/config.yaml > $WORKDIR/config.yaml
	sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/samples.tsv > $WORKDIR/samples.tsv
	cp $PIPELINE_HOME/resources/cluster.json $WORKDIR/cluster.json
	cp -r $PIPELINE_HOME/workflow/scripts $WORKDIR/

elif [ "$RUNMODE" == "unlock" ];then

	check_workdir_exists
	snakemake --unlock \
	-s $SNAKEFILE \
	--configfile $CONFIGYAML \
	--directory $WORKDIR

elif [ "$RUNMODE" == "cluster" ];then

	check_workdir_exists
	cd $WORKDIR
	cat << EOF > submit.sh 
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
snakemake \
-s $SNAKEFILE \
--directory $WORKDIR \
--use-singularity \
--singularity-args "'-B $WORKDIR,$INPUTDIR'" \
--configfile $CONFIGYAML \
--use-envmodules \
--printshellcmds \
--latency-wait 120 \
--cluster-config $CLUSTERJSON \
--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} " \
-j 500 \
--rerun-incomplete \
--keep-going \
--stats "${WORKDIR}/snakemake.stats" \
2>&1 | tee "${WORKDIR}/snakemake.log"

if [ "\$?" -eq "0" ];then
  snakemake -s $SNAKEFILE \
  --directory $WORKDIR \
  --report ${WORKDIR}/runslurm_snakemake_report.html \
  --configfile ${WORKDIR}/config.yaml 
fi

bash <(curl https://raw.githubusercontent.com/CCBR/Tools/master/Biowulf/gather_cluster_stats_biowulf.sh 2>/dev/null) ${WORKDIR}/snakemake.log > ${WORKDIR}/snakemake.log.HPC_summary.txt

EOF

	sbatch submit.sh
fi
