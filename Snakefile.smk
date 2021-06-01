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

GLOBAL_REF_PATH = "/mnt/references/"

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])




###########################################
# DEFINITION OF VARIABLES
#
cfg = pd.DataFrame(config)

REF_DIR = define_variable(cfg, "REF_DIR")
PROJECT_NAME = define_variable(cfg, "PROJECT_NAME")
PROJECT_DIR = define_variable(cfg, "PROJECT_DIR")
INPUTS_DIR = define_variable(cfg, "INPUTS_DIR")
ADIR = define_variable(cfg, "ANALYSIS_DIR")
SAMPLE = "full_name"

##### Target rules #####

rule all:
     input: ediing_sites = expand(ADIR + "/per_sample_results/{full_name}.editing_sites.tsv", full_name = cfg[SAMPLE].tolist())

##### Modules #####

include: "rules/get_editing_sites.smk"


