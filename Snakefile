import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("7.2.1")

configfile: "config.json"


##### Config and reference processing #####
#

f = open(os.path.join(config["globalResources"],"reference_info","reference2.json"),)
reference_dict = json.load(f)
f.close()
config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
if len(config["species_name"].split(" (")) > 1:
    config["species"] = config["species_name"].split(" (")[1].replace(")","")


reference_directory = os.path.join(config["globalResources"],config["organism"],config["reference"])
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")


wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),

##### Target rules #####

rule all:
     input: expand("results/{sample}.mismatch_cov_tab.tsv", sample = sample_tab.sample_name),

##### Modules #####

include: "rules/get_editing_sites.smk"
