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
