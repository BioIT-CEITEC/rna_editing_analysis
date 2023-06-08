#############################################################
# wrapper for rule: get_binned_gc_content
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: get_binned_gc_content \n##\n")
f.close()

shell.executable("/bin/bash")
command = "bedtools nuc -fi " + \
          " " + snakemake.input.ref + \
          " -bed " + snakemake.input.region_bed + \
          " > " + snakemake.output.gc_content_file + ".tmp"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)


command = "Rscript  " + os.path.abspath(os.path.dirname(__file__)) + "/get_binned_gc_content.R" \
                        + " " + snakemake.output.gc_content_file + ".tmp" \
                        + " " + snakemake.output.gc_content_file \
                        + " 2>> " + log_filename


f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()
shell(command)

command = "rm " + snakemake.output.gc_content_file + ".tmp"
shell(command)
