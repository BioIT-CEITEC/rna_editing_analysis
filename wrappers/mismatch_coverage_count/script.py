######################################
# wrapper for rule: mismatch_coverage_count
######################################
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: mismatch_coverage_count \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()


command = "bedtools coverage -a " + snakemake.input.bed + " -b " + snakemake.input.bam + " -counts " + snakemake.output.tsv + " >> " + log_filename + " 2>&1"
with open(log_filename,'at') as f:
   f.write("## COMMAND: "+command+"\n")
shell(command)



