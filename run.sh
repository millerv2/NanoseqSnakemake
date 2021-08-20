#!/bin/bash
# runmode : dryrun or unlock or cluster or init
RUNMODE=$1
# workdir : dir where all the results will go
WORKDIR=$2
# inputdir : dir where all the input fastqs are located
INPUTDIR=$3

DEFAULT_SINGULARITY_BINDS="/gpfs/gsfs11/users/$USER,/data/$USER"

PIPELINE_HOME="$(readlink -f $(dirname "${BASH_SOURCE[0]}"))"
SNAKEFILE="$PIPELINE_HOME/workflow/Snakefile"
CONFIGYAML="$WORKDIR/config.yaml"
CLUSTERJSON="$WORKDIR/cluster.json"
SCRIPTNAME="$0"
SCRIPTBASENAME=$(readlink -f $(basename $0))

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


function usage() { 

# This function prints generic usage of the wrapper script.

    cat << EOF
${SCRIPTBASENAME}
--> run NanoSeqPipeline
USAGE:
  bash ${SCRIPTNAME} <RUNMODE> <RESULTSDIR> <SAMPLESDIR>
Required Arguments:
1.  RUNMODE: [Type: String] Valid options:
    *) init : initialize workdir
    *) cluster : run with slurm
    *) dryrun : dry run snakemake to generate DAG
    *) unlock : unlock workdir if locked by snakemake
2.  RESULTSDIR: [Type: String]: 
	Absolute path to the output folder with write permissions.
3.  SAMPLESDIR: [Type: String]: 
	Absolute path to the folder where the sample fastqs are. You only need read permissions to his folder.

EOF

}

if [ $# -ne 3 ]; then usage; exit 1;fi

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
	## Archive previous run files
	if [ -f ${WORKDIR}/snakemake.log ];then 
		modtime=$(stat ${WORKDIR}/snakemake.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
		mv ${WORKDIR}/snakemake.log ${WORKDIR}/logs/snakemake.${modtime}.log
		if [ -f ${WORKDIR}/snakemake.log.HPC_summary.txt ];then 
		mv ${WORKDIR}/snakemake.log.HPC_summary.txt ${WORKDIR}/logs/snakemake.${modtime}.log.HPC_summary.txt
		fi
		if [ -f ${WORKDIR}/snakemake.stats ];then 
		mv ${WORKDIR}/snakemake.stats ${WORKDIR}/logs/snakemake.${modtime}.stats
		fi
	fi
	nslurmouts=$(find ${WORKDIR} -maxdepth 1 -name "slurm-*.out" |wc -l)
	if [ "$nslurmouts" != "0" ];then
		for f in $(ls ${WORKDIR}/slurm-*.out);do mv ${f} ${WORKDIR}/logs/;done
	fi

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
--singularity-args "'-B $DEFAULT_SINGULARITY_BINDS,$WORKDIR,$INPUTDIR'" \
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
