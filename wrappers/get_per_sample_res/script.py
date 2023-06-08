#############################################################
# wrapper for rule: cnv_computation
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: per_sample_results \n##\n")
f.close()

command = "Rscript  " + os.path.abspath(os.path.dirname(__file__)) + "/get_per_sample_res.R" \
                       + " " + snakemake.output.CNV_res \
                       + " " + snakemake.input.all_res_prob_tab\
                       + " " + snakemake.wildcards.sample_name\
                       + " 2>> " + snakemake.log.run

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()
shell(command)

# command = "touch " + snakemake.output.all_res_prob_tab
# shell(command)