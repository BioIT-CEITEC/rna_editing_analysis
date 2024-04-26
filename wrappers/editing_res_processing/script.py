######################################
# wrapper for rule: merge_variant_callers
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: process editing \n##\n")
f.close()


# after_merge_processing
command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/editing_res_processing.R "
            # " >> " + log_filename + " 2>&1"
#
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
# f.close()

shell(command)
