######################################
# wrapper for rule: get_editing_map
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: get_editing_map \n##\n")
f.close()



command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/create_tab.R "+\
            snakemake.input.second_tmp_mileup + " " +\
            snakemake.input.bam + " " +\
            snakemake.output.rdata + \
            " >> " + snakemake.log.run + " 2>&1"

f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
