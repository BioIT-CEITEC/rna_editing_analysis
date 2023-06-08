#############################################################
# wrapper for rule: get_binned_bed_from_dict
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: get_binned_bed_from_dict \n##\n")
f.close()

shell.executable("/bin/bash")


command = "Rscript  " + os.path.abspath(os.path.dirname(__file__)) + "/get_binned_bed_from_dict.R" \
                        + " " + snakemake.input.ref_dict \
                        + " " + snakemake.output.bed\
                        + " " + str(snakemake.params.window_size) \
                        + " 2>> " + log_filename


f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()
shell(command)

