import os
import pandas as pd
import json

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]

##### BioRoot utilities #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

read_pair_tags = BR.set_read_pair_qc_tags() # ["SE"] / ["R1", "R2"]
pair_tag = BR.set_read_pair_tags() # [""] / ["_R1", "_R2"]
paired = BR.set_paired_tags() # "SE" / "PE"

config = BR.load_organism()

#### FOLDERS
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

####################################
# SEPARATE RULES
include: "rules/edit_analysis.smk"

####################################
# RULE ALL
rule all:
    input: "reports/editing_analysis_report.html"
