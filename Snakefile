import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

#GLOBAL_REF_PATH = "/mnt/references/"
GLOBAL_REF_PATH = "/mnt/ssd/ssd_3/references/"

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

###########################################
# DEFINITION OF VARIABLES
#
cfg = pd.DataFrame(config)

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
wildcard_constraints:
    sample = "|".join(sample_tab.sample_name)

##### Target rules #####

rule all:
     input: editing_sites = expand("per_sample_results/{sample}.editing_sites.tsv", sample = sample_tab.sample_name)

##### Modules #####

include: "rules/get_editing_sites.smk"


