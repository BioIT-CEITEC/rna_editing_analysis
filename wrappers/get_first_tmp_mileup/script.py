######################################
# wrapper for rule: get_first_tmp_mileup
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: get_first_tmp_mileup \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("samtools --version 2>&1 | grep samtools",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = "mkdir -p " + os.path.dirname(snakemake.output.first_tmp_mileup)
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
						#+ snakemake.input.region_bed					
command = "samtools mpileup --positions " + snakemake.input.region_bed[0] + \
                                          " -f " + snakemake.input.ref[0] + " " + snakemake.input.bam + \
                                          " 2>> " +snakemake.log.run + " > " + snakemake.output.first_tmp_mileup
f = open(snakemake.log.run, 'at')                                             
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


