######################################
# wrapper for rule: get_first_filtered_positions
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: get_first_filtered_positions \n##\n")
f.close()


command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/filter_pos.R "+\
            snakemake.input.first_tmp_mileup + " " +\
            snakemake.output.first_filtered_pos + " " +\
            " >> " + snakemake.log.run + " 2>&1"

f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


