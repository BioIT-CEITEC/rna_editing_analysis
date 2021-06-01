import math
import subprocess
import json
import re
import os.path
import pandas as pd
from snakemake.utils import R
from snakemake.utils import report
from os.path import split
from helper import define_variable
from snakemake.utils import min_version

include: "prepare_reference.smk"

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

##### Target rules #####

rule all:
     input: ediing_sites = expand(ADIR + "/per_sample_results/{full_name}.editing_sites.tsv", full_name = cfg[SAMPLE].tolist())

##### Modules #####

include: "rules/get_editing_sites.smk"


