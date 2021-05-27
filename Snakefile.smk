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

include: "prepare_reference.smk"


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

# REF_DIR = "/mnt/ssd/ssd_3/references"
# PROJECT_NAME = cfg['project'].tolist()[0].replace("/",".")
# PROJECT_DIR = os.path.join(cfg['project_owner'].tolist()[0]
#                     ,"sequencing_results"
#                     ,"projects"
#                     ,cfg['project'].tolist()[0])
# INPUTS_DIR = os.path.join(PROJECT_DIR,"input_files")
# ADIR = os.path.join(PROJECT_DIR,cfg['analysis_name'].tolist()[0])


####################################
# FINAL RESULTING FILES
#



# def bams_input(wildcards):
#     if cfg["material"].tolist()[0] != "RNA":
#         return set(expand(INPUTS_DIR+"/mapped/{full_name}.bam",full_name = cfg["full_name"].tolist()))
#     else:
#         return set(expand(INPUTS_DIR+"/mapped/{full_name}.RNAsplit.bam",full_name = cfg["full_name"].tolist()))
#
# def raw_fastqc_input(wildcards):
#     if cfg.loc[cfg.full_name ==  wildcards.full_name,"material"].min() != "RNA":
#         return INPUTS_DIR+"/mapped/"+wildcards.full_name+".bam"
#     else:
#         return INPUTS_DIR+"/mapped/"+wildcards.full_name+".RNAsplit.bam"
#
# def raw_fastqc_input_index(wildcards):
#     if cfg.loc[cfg.full_name ==  wildcards.full_name,"material"].min() != "RNA":
#         return INPUTS_DIR+"/mapped/"+wildcards.full_name+".bam.bai"
#     else:
#         return INPUTS_DIR+"/mapped/"+wildcards.full_name+".RNAsplit.bam.bai"


 #rule final_RNA_edit_report:
  #   input: reports = ADIR+"/EditingIndex.csv"
 #    output: html = ADIR + "/" + PROJECT_NAME + "_finished"
 #    shell:
  #       "touch {output.html}"

#rule RNA_edit_postprocessing:
   # input: #reports = ADIR+"/EditingIndex.csv",
          # editing_maps = expand(ADIR+"/per_sample_results/{full_name}.edit_tab.Rdata",full_name = cfg['full_name'].tolist())
   # output: html = ADIR + "/" + PROJECT_NAME + "_finished"
   # shell:
     #   "touch {output.html}"

#rule get_full_editing_map:
 #   input: bam = INPUTS_DIR + "/mapped/{full_name}.bam",
          # ref = lambda wildcards:  expand("{dir}/{organism}/{ref}/seq/{ref}.fa",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == #wildcards.full_name,"reference"].min()),
          # region_bed = lambda wildcards:  expand("{dir}/{organism}/{ref}/intervals/RNA/all_As.pos",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc
#[cfg  [SAMPLE] == wildcards.full_name,"reference"].min()),
    #output: rdata = ADIR + "/per_sample_results/{full_name}.edit_tab.Rdata",
   # log:    run = ADIR + "/sample_logs/{full_name}/get_full_editing_map.log"
   # params: res_base = ADIR + "/results/{full_name}",
   # threads: 10
   # conda:  "../wraps/RNA_editing/get_full_editing_map/env.yaml"
   # script: "../wraps/RNA_editing/get_full_editing_map/script.py"
#rule all:
  # input: 
    #   rdata = ADIR + "/per_sample_results/{full_name}.edit_tab.Rdata"
rule all:
     input: ediing_sites = expand(ADIR + "/per_sample_results/{full_name}.editing_sites.tsv", full_name = cfg[SAMPLE].tolist())


rule get_editing_sites:
     input:
       bam = INPUTS_DIR + "/mapped/{full_name}.bam",
       gtf = lambda wildcards:  expand("{dir}/{organism}/{ref}/annot/{ref}.gtf",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"reference"].min())[0],
       known_editing = lambda wildcards:  expand("{dir}/{organism}/{ref}/other/known_editing_sites/Human_AG_all_hg38_v2.csv",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"reference"].min())[0],
       snp = lambda wildcards:  expand("{dir}/{organism}/{ref}/other/snp/{ref}.snp.bed",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"reference"].min())[0],
     output:
       fa = ADIR + "/per_sample_results/{full_name}.seq.fa",
       editing_sites = ADIR + "/per_sample_results/{full_name}.editing_sites.tsv",
     log:
       run = ADIR + "/sample_logs/{full_name}/get_editing_sites.log"
     threads: 20
     conda: "../wraps/RNA_editing/get_editing_sites/env.yaml"
     script: "../wraps/RNA_editing/get_editing_sites/script.py"


#rule get_first_tmp_mileup:
   # input: 
    #   bam = INPUTS_DIR + "/mapped/{full_name}.bam",
    #   ref = lambda wildcards:  expand("{dir}/{organism}/{ref}/seq/{ref}.fa",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"reference"].min()),
     #  region_bed = lambda wildcards:  expand("{dir}/{organism}/{ref}/intervals/RNA/all_As_genes_strand_wise.pos",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"reference"].min())
  #  output:
     #   first_tmp_mileup = ADIR + "/per_sample_results/{full_name}.mpileup.tmp1.tsv"
   # log:
   #   run = ADIR + "/sample_logs/{full_name}/get_first_tmp_mileup.log"
   # threads: 10
   # conda:  "../wraps/RNA_editing/get_first_tmp_mileup/env.yaml"
  #  script: "../wraps/RNA_editing/get_first_tmp_mileup/script.py"

#rule get_first_filtered_positions:
     #input:
    #   first_tmp_mileup = ADIR + "/per_sample_results/{full_name}.mpileup.tmp1.tsv"
    # output:
     #  first_filtered_pos = ADIR + "/per_sample_results/{full_name}.filtered_pos1.tsv"
    # log:
    #  run = ADIR + "/sample_logs/{full_name}/get_first_filtered_positions.log"
    # threads: 10
   #  conda:  "../wraps/RNA_editing/get_first_filtered_positions/env.yaml"
    # script: "../wraps/RNA_editing/get_first_filtered_positions/script.py"
     

#rule get_second_tmp_mileup: 
    #  input:
    #    first_filtered_pos = ADIR + "/per_sample_results/{full_name}.filtered_pos1.tsv",
     #   ref = lambda wildcards:  expand("{dir}/{organism}/{ref}/seq/{ref}.fa",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"reference"].min()),
     #   bam = INPUTS_DIR + "/mapped/{full_name}.bam"
     # output:
     #   second_tmp_mileup =  ADIR + "/per_sample_results/{full_name}.mileup.tmp2.tsv"
     # log:
     #  run = ADIR + "/sample_logs/{full_name}/get_second_tmp_mileup.log"
     # threads: 10
     # conda:  "../wraps/RNA_editing/get_second_tmp_mileup/env.yaml"
    #  script: "../wraps/RNA_editing/get_second_tmp_mileup/script.py"

#rule get_editing_map:
    # input:
     #   second_tmp_mileup =  ADIR + "/per_sample_results/{full_name}.mileup.tmp2.tsv",
     #   bam = INPUTS_DIR + "/mapped/{full_name}.bam"
    # output:
     #   rdata = ADIR + "/per_sample_results/{full_name}.edit_tab.Rdata"     
    # log:    run = ADIR + "/sample_logs/{full_name}/get_editin
#g_map.log"
   #  threads: 10
   #  conda:  "../wraps/RNA_editing/get_editing_map/env.yaml"
   #  script: "../wraps/RNA_editing/get_editing_map/script.py"


#rule RNA_edit_with_a2i:
#    input:  bam = expand(INPUTS_DIR+"/mapped/{full_name}.bam",full_name = cfg['full_name'].tolist()),
            #bai = expand(INPUTS_DIR+"/mapped/{full_name}.bam.bai",full_name = cfg['full_name'].tolist()),
            #refdir = lambda wildcards: expand("{dir}/{organism}/{ref}/",dir = REF_DIR,organism = cfg["organism"].min(),ref = cfg["reference"].min()),
            #fa  = lambda wildcards:  expand("{dir}/{organism}/{ref}/seq/{ref}.fa",dir = REF_DIR,organism = cfg["organism"].min(),ref = cfg["reference"].min()),
            #fai = lambda wildcards:  expand("{dir}/{organism}/{ref}/seq/{ref}.fa.fai",dir = REF_DIR,organism = cfg["organism"].min(),ref = cfg["reference"].min()),
            #snp = lambda wildcards:  expand("{dir}/{organism}/{ref}/other/a2i_data/ensembl/SNPs150.bed.gz",dir = REF_DIR,organism = cfg["organism"].min(),ref = cfg["reference"].min()),
            #exp = lambda wildcards:  expand("{dir}/{organism}/{ref}/other/a2i_data/ensembl/GTExGeneExpression.bed.gz",dir = REF_DIR,organism = cfg["organism"].min(),ref = cfg["reference"].min()),
            #ref = lambda wildcards:  expand("{dir}/{organism}/{ref}/other/a2i_data/ensembl/RefSeqCuratedAnnotation.bed.gz",dir = REF_DIR,organism = cfg["organism"].min(),ref = cfg["reference"].min()),
            #reg = lambda wildcards:  expand("{dir}/{organism}/{ref}/other/a2i_data/ensembl/AluRegions.bed.gz",dir = REF_DIR,organism = cfg["organism"].min(),ref = cfg["reference"].min()),
            # bed = lambda wildcards:  expand("{dir}/{organism}/{ref}/intervals/{kit}/{kit}.bed",\
                                             #dir = REF_DIR,\
                                             #organism = cfg.loc[cfg.full_name == wildcards.full_name,"organism"].min(),\
                                            # ref = cfg.loc[cfg.full_name == wildcards.full_name,"reference"].min(),\
                                             #kit = cfg.loc[cfg.full_name == wildcards.full_name,"library_scope"].min())
    #output: sum = ADIR+"/EditingIndex.csv"
    #log:    run = ADIR+"/EditingIndex.log"
   # params: AF_threshold = 0.01,
           # paired = cfg['paired'].min(),
            #strand = 'unstr' if 'strandness' not in cfg else cfg['strandness'].min()
   # threads: 20,
    #conda:  "../wraps/RNA_editing/RNA_edit_with_a2i/env.yaml"
    #script: "../wraps/RNA_editing/RNA_edit_with_a2i/script.py"
