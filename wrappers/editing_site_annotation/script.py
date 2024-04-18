######################################
# wrapper for rule: variant_annotation
######################################
import os
import subprocess
import re
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: vep annot \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("vep 2>&1 | grep \"ensembl-vep\" | cut -f 2 -d \":\"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: ensembl-vep  :  "+version+"\n")
f.close()

assembly = re.sub("\-.*","",snakemake.params.ref_name)

if snakemake.threads > 1:
    fork_text = "--fork " + str(snakemake.threads)
else:
    fork_text = ""

cache_version = 0

organism = snakemake.params.organism_name

if organism == "homsap":
    organism = "homo_sapiens"

if os.path.isdir(snakemake.params.vep_dir):

    merged = ""
    vep_dir = snakemake.params.vep_dir + "/" + organism
    for dir in os.listdir(vep_dir):
         if re.search("^[0-9]+_",dir):
             num = int(re.sub("_.*","",dir))
             if num > cache_version:
                 cache_version = num

    command = "vep --dir " + snakemake.params.vep_dir + \
              " --everything " + \
              " --fasta " + snakemake.params.ref + \
              " --species " + organism + \
              " --offline --assembly " + assembly + \
              " --cache " +\
              " --cache_version " + str(cache_version) + \
              " --input_file " + snakemake.input.tsv_for_vep + \
              " --output_file " + snakemake.output.annotated + \
              " --force_overwrite " + \
              fork_text + " >> " + log_filename + " 2>&1"

    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)




