######################################
# wrapper for rule: control_freec
######################################
import os
from re import sub
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: control freec \n##\n")
f.close()

shell.executable("/bin/bash")

version = sub(" : a.*", "", str(subprocess.Popen("freec 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8'))
f = open(log_filename, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

#not implemented --targeted
# if snakemake.params.library_scope == "wgs":
#     scope = ""
# else:
#     scope = " --exome --callRegions " + snakemake.input.regions_gz + " "

f = open(snakemake.input.config_template, 'at')
# Read in the config template
with open(snakemake.input.config_template, 'r') as file :
  filedata = file.read()

# Replace specific parameters
filedata = filedata.replace('X_window_X', str(snakemake.params.window_size))
filedata = filedata.replace('X_ref_fai_X', snakemake.input.ref_fai)
filedata = filedata.replace('X_GC_profile_file_X', snakemake.input.GC_profile_file)
filedata = filedata.replace('X_output_dir_X', os.path.dirname(snakemake.output.config))
filedata = filedata.replace('X_input_bam_X', snakemake.input.tumor)

if hasattr(snakemake.input, 'normal'):
  filedata = filedata.replace('#normal_X__', '')
  filedata = filedata.replace('X_input_control_bam_X', snakemake.input.normal)

# Write the specific output
with open(snakemake.output.config, 'w') as file:
  file.write(filedata)


# command = "samtools view -b " + snakemake.input.tumor + " {1,10,11,12,13,14,15,16,17,18,19,2,21,22,3,4,5,6,7,8,9,X,Y} > " + snakemake.input.tumor + " 2>&1"
#
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

command = "freec -conf " + snakemake.output.config + " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# remove workspace - lots of files messing copying speeds
shell("cp " + snakemake.params.control_freec_outfile + " " + snakemake.output.CNVs_vcf)






# command = "gatk MergeVcfs " + \
#         " -I "+ snakemake.params.dir + "/results/variants/somatic.indels.vcf.gz "+\
#         " -I "+ snakemake.params.dir + "/results/variants/somatic.snvs.vcf.gz "+\
#         " -R "+ snakemake.input.ref +\
#         " -D " + snakemake.input.dict +\
#         " -O " + snakemake.output.vcf +\
#         " >> " + log_filename + " 2>&1"
#
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)


# command =  "cp " + snakemake.params.dir + "/results/variants/somatic.indels.vcf.gz " + str(snakemake.output.vcf) + ".gz"
#
# f = open(log_filename, 'a+')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
#
# command =  "gunzip " + str(snakemake.output.vcf) + ".gz"
#
# f = open(log_filename, 'a+')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
