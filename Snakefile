import os
import pandas as pd

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

# ####################################
# # VARIALBES FROM CONFIG
# used_SV_callers = []
# used_CNV_callers = []
# if config["use_gatk_cnv"]:
#     used_CNV_callers.append("gatk_cnv")
# if config["use_cnvkit"]:
#     used_CNV_callers.append("cnvkit")
# if config["use_jabCoNtool"]:
#     used_CNV_callers.append("jabCoNtool")
# if config["use_control_freec"]:
#     used_CNV_callers.append("control_freec")
# if config["use_manta"]:
#     used_SV_callers.append("manta")
# if config["use_gridss"]:
#     used_SV_callers.append("gridss")

#
# wildcard_constraints:
#      tumor_normal = "tumor|normal|sample",


####################################
# SEPARATE RULES
include: "rules/edit_analysis.smk"
# include: "rules/vep.smk"



####################################
# RULE ALL
rule all:
    input: "reports/editing_analysis_report.html"
