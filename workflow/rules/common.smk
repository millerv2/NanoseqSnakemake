from snakemake.utils import validate
import pandas as pd
import yaml

##### load config and sample sheets #####

#configfile: "/data/millerv2/NanoseqSnakemake/config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")
#print(workflow.overwrite_configfiles)
#CONFIGFILE = str(workflow.overwrite_configfiles[0])
sample_dir=config['sample_dir']
base_dir = config['base_dir']
samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
#validate(samples, schema="../schemas/samples.schema.yaml")
SAMPLES = list(samples['sample'])
application_type = list(samples['application'])
sample_to_application = { SAMPLES[i] : application_type[i] for i in range(len(SAMPLES))}

## Load cluster.json
try:
    CLUSTERJSON = config["clusterjson"]
except KeyError:
    CLUSTERJSON = os.path.join(base_dir,"cluster.json")
# check_readaccess(CLUSTERJSON)
with open(CLUSTERJSON) as json_file:
    CLUSTER = json.load(json_file)

getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER and "threads" in CLUSTER[rname] else int(CLUSTER["__default__"]["threads"])