######################################
# wrapper for rule: get_editing_sites
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

f = open(snakemake.log.run, 'w')
f.write("\n##\n## RULE: get_editing_sites \n##\n")
f.close()

command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/rna_editing_bam.R "+\
            snakemake.input.bam + " " +\
            snakemake.input.gtf + " " +\
            snakemake.input.known_editing + " " +\
            snakemake.input.snp + " " +\
            snakemake.params.chrm + " " +\
            snakemake.log.run + " " +\
            snakemake.output.fa + " " +\
            snakemake.output.editing_sites + " " +\
            " >> " + snakemake.log.run + " 2>&1"

f = open(snakemake.log.run, 'a')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


