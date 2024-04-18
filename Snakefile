import os
import pandas as pd
import json

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]

##### Config processing #####
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

#### Reference info processing
#### Setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"))
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]

#### FOLDERS
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

####################################
# SEPARATE RULES
include: "rules/edit_analysis.smk"

####################################
# RULE ALL
rule all:
    input: "reports/editing_analysis_report.html"
