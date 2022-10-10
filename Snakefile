import os
import pandas as pd
from snakemake.utils import min_version

min_version("7.2.1")

configfile: "config.json"


##### Config and reference processing #####
#
GLOBAL_REF_PATH = "/mnt/nfs/shared/S3acgt/resources"

reference_directory = os.path.join(GLOBAL_REF_PATH,"organisms",config["organism"],config["reference"])
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")


wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),


##### Target rules #####

rule all:
     input: expand("results/{sample}.mismatch_tab.tsv", sample = sample_tab.sample_name),

##### Modules #####

include: "rules/get_editing_sites.smk"

