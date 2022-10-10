######################################
# wrapper for rule: get_editing_sites
######################################
import os
from snakemake.shell import shell

f = open(snakemake.log.run, 'w')
f.write("\n##\n## RULE: get_editing_sites \n##\n")
f.close()

command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/rna_editing_bam.R "+\
            snakemake.input.bam + " " +\
            snakemake.input.gtf + " " +\
            snakemake.params.chr + " " +\
            snakemake.output.mismatch_tab + " " +\
            " >> " + snakemake.log.run + " 2>&1"

f = open(snakemake.log.run, 'a')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


