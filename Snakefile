import pandas as pd
import json
import os
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

config = BR.load_organism()

tools = BR.load_tooldir()


#### FOLDERS
wildcard_constraints:
    sample = "|".join(sample_tab.sample_name)

####################################
# SEPARATE RULES
include: "rules/edit_analysis.smk"

####################################
# RULE ALL
rule all:
    input: expand("jacusa_call/{sample}_potential_edit_sites.tsv", sample=sample_tab.sample_name)