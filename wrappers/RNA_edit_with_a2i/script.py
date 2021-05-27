##########################################
# wrapper for rule: RNA_edit_with_a2i #
##########################################
import os
import sys
import math
import random
import string
import subprocess
import re
from snakemake.shell import shell


f = open(snakemake.log.run, 'wt')
f.write("\n##\n## RULE: RNA_edit_with_a2i \n##\n")
f.close()

shell.executable("/bin/bash")

## Don't change this part unless you know what it does
env_command = "export"
env_command+= " TEMPDIR=/mnt/ssd/ssd_1/tmp/"
env_command+= " SINGULARITY_CACHEDIR=/mnt/ssd/ssd_1/tmp/"
env_command+= " SINGULARITY_LOCALCACHEDIR=/mnt/ssd/ssd_1/tmp/"
# env_command+= " SINGULARITY_DISABLE_CACHE=True"

## Here is declaration of used container (docker, singularity, etc.)
image = "/home/325073/RNAEditingIndexer/a2i_editing.simg"
f = open(snakemake.log.run, 'at')
f.write("\n## IMAGE: "+image+"\n##\n")
f.close()

## This part must be tailored to particular rule for all needed files following pattern: path_outside_container:path_inside_container
## where file was defined in the rule using expand() function it must follow [0] (e.g., snakemake.input.bam[0])
# following temporary directory must be specified inside this script to be unique for each instance of the rule
temp_dir = "/mnt/ssd/ssd_1/tmp/a2i_"+''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(10))+"/"
## inputs
outside_input_dir = temp_dir+"in/"
inside_input_dir = "/inputs/"
inside_refdir = "/ref_data/"
#inside_bed = inside_input_dir+os.path.basename(snakemake.input.bed[0])
## outputs
outside_output_dir = temp_dir+"out/"
inside_output_dir = "/workdir/"
inside_sum = inside_output_dir+os.path.basename(snakemake.output.sum)
## binding of files/dirs from outside of container to inside container files/dirs
image_options = " -B "+snakemake.input.refdir[0]+":"+snakemake.input.refdir[0]
# image_options = " -B "+snakemake.input.fa[0]+":"+inside_refdir+os.path.basename(snakemake.input.fa[0])
# image_options = " -B "+snakemake.input.fai[0]+":"+inside_refdir+os.path.basename(snakemake.input.fai[0])
# image_options = " -B "+snakemake.input.snp[0]+":"+inside_refdir+os.path.basename(snakemake.input.snp[0])
# image_options = " -B "+snakemake.input.exp[0]+":"+inside_refdir+os.path.basename(snakemake.input.exp[0])
# image_options = " -B "+snakemake.input.reg[0]+":"+inside_refdir+os.path.basename(snakemake.input.reg[0])
# image_options = " -B "+snakemake.input.ref[0]+":"+inside_refdir+os.path.basename(snakemake.input.ref[0])
image_options+= " -B "+outside_input_dir+":"+inside_input_dir
image_options+= " -B "+outside_output_dir+":"+inside_output_dir
#image_options+= " -B "+snakemake.params.outdir+":"+inside_output_dir
image_options+= " -c "

debug = ""
## Allow this if you want to see DEBUG (VERBOSE) informations
# debug = "--debug"

## This should not need any change
image_command = "singularity "+debug+" exec "+image_options+" "+image

## Here must be specified a command for particular tool in container (used files must be inside container)
tool_command = "RNAEditingIndex"
tool_command+= " -d "+inside_input_dir
tool_command+= " -f "+".bam"
tool_command+= " -l "+inside_output_dir+"/logs/"
tool_command+= " -o "+inside_output_dir+"/cmpileups/"
tool_command+= " -os "+inside_output_dir+"/summary/"
tool_command+= " --genome_fasta "+snakemake.input.fa[0]
tool_command+= " --genome "+"UserProvided"
tool_command+= " --regions "+snakemake.input.reg[0]
tool_command+= " --refseq "+snakemake.input.ref[0]
tool_command+= " --snps "+snakemake.input.snp[0]
tool_command+= " --genes_expression "+snakemake.input.exp[0]
# tool_command+= " --verbose"
if snakemake.params.strand.lower() == 'fwd' or snakemake.params.strand.lower() == 'rev':
  tool_command+= " --stranded"
if snakemake.params.paired.lower() == 'pe':
  tool_command+= " --paired"

## Create temporary directory 
command = "mkdir -p "+outside_input_dir+" "+outside_output_dir+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

## Link input bams inside temporary dir
list_of_inputs = snakemake.input.bam+snakemake.input.bai
for inp in list_of_inputs:
  command = "ln "+inp+" "+outside_input_dir+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

## Main command
command = env_command + " && " + image_command + " " + tool_command + " >> " + snakemake.log.run + " 2>&1"
print(command)
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n\n")
f.close()
shell(command)

## List all the outputs for check
# command = "ls -l "+snakemake.params.outdir+" >> "+snakemake.log.run+" 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

## Move needed output file
command = "mv "+outside_output_dir+"/summary/"+os.path.basename(snakemake.output.sum)+" "+snakemake.output.sum+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

## Move needed output file
# command = "mv "+outside_output_dir+"/summary/"+os.path.basename(snakemake.output.sum)+" "+snakemake.output.sum+" >> "+snakemake.log.run+" 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

## Clean temporary directory
command = "rm -rf "+temp_dir+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
