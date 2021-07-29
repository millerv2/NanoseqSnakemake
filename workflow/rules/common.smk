from snakemake.utils import validate
import pandas as pd

##### load config and sample sheets #####

#configfile: "/data/millerv2/NanoseqSnakemake/config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")
#print(workflow.overwrite_configfiles)
#CONFIGFILE = str(workflow.overwrite_configfiles[0])

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
#validate(samples, schema="../schemas/samples.schema.yaml")
