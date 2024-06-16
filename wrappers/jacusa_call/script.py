###########################################
###########################################
# wrapper for rule: jacusa_call
###########################################

import subprocess
from os.path import dirname
from snakemake.shell import shell

shell.executable("/bin/bash")
log_filename = str(snakemake.log)

version = str(subprocess.Popen("conda list", shell = True, stdout = subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, "wt")
f.write("\n##\n## CONDA: " + version + "\n")
f.close()

f = open(log_filename, "at")
f.write("\n##\n## RULE: JACUSA2 call \n##\n")
f.close()

command = "java -Xmx130g -jar " + snakemake.params.jacusa + " call-1 -p " + str(snakemake.threads) + \
          " -P " + snakemake.params.strand + " -R " + snakemake.input.ref +" -r " + snakemake.output.jacusa + \
          " -F 1024 -m " + str(snakemake.params.min_mapq)  + " -filterNH " + str(snakemake.params.filternh) + \
          " -filterNM " + str(snakemake.params.filternm) + " " + snakemake.input.bam + " >> " + log_filename + " 2>&1"
f = open(log_filename, "at")
f.write("\n##\n## COMMAND: " + command + "\n")
f.close()
shell(command)