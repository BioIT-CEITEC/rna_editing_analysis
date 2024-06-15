###########################################
# wrapper for rule: jacusa_helper
###########################################

import subprocess
import os
from os.path import dirname
from snakemake.shell import shell

shell.executable("/bin/bash")
log_filename = str(snakemake.log)

version = str(subprocess.Popen("conda list", shell = True, stdout = subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, "at")
f.write("\n##\n## CONDA: " + version + "\n")
f.close()

command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/jacusa_helper.R "+\
            dirname(snakemake.input.raw) + " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)