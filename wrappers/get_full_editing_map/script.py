######################################
# wrapper for rule: get_full_editing_map
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: get_full_editing_map \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("samtools --version 2>&1 | grep samtools",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = "mkdir -p " + os.path.dirname(snakemake.output.rdata)
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "samtools mpileup --positions " + snakemake.input.region_bed[0] + \
                                          " -f " + snakemake.input.ref[0] + " " + snakemake.input.bam + \
                                          " 2>> " +snakemake.log.run + " > " + snakemake.params.res_base + ".mpileup.tmp.tsv"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/filter_pos.R "+\
            snakemake.params.res_base + ".mpileup.tmp.tsv " +\
            snakemake.params.res_base + ".filtered_pos.tsv" +\
            " >> " + snakemake.log.run + " 2>&1"

f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "samtools mpileup --output-QNAME --positions " + snakemake.params.res_base + \
                                                         ".filtered_pos.tsv -f " + snakemake.input.ref[0] + " " + snakemake.input.bam + \
                                                         " 2>> " +snakemake.log.run + " > " + snakemake.params.res_base + ".mpileup.tmp2.tsv"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/create_tab.R "+\
            snakemake.params.res_base + ".mpileup.tmp2.tsv " +\
            snakemake.input.bam + " " +\
            snakemake.output.rdata + \
            " >> " + snakemake.log.run + " 2>&1"

f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
