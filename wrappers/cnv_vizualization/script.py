#############################################################
# wrapper for rule: cnv_visualization
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'wt')
f.write("\n##\n## RULE: cnv_visualization \n##\n")
f.close()


command = " Rscript  " + os.path.abspath(os.path.dirname(__file__)) + "/cnv_visualization.R " \
                       + " " + snakemake.input.rdata \
                       + " " + snakemake.output.pdf \
                       + " 2>> " + snakemake.log.run
f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
