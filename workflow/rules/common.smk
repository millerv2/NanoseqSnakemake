from snakemake.utils import validate
import pandas as pd

##### load config and sample sheets #####

#configfile: "config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")
CONFIGFILE = str(workflow.overwrite_configfiles[0])

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
#validate(samples, schema="../schemas/samples.schema.yaml")
