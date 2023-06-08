######################################
# wrapper for rule: germline_strelka
######################################
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: cnv cnvkit fix and segment \n##\n")
f.close()

shell.executable("/bin/bash")


command = "cnvkit.py fix" +  \
          " " + snakemake.input.targetcoverage + \
          " " + snakemake.input.antitargetcoverage + \
          " " + snakemake.input.cnv_reference + \
          " -o " + snakemake.output.fix + \
          " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

command = "cnvkit.py segment" +  \
          " " + snakemake.output.fix + \
          " -o " + snakemake.output.segments + \
          " >> " + log_filename + " 2>&1"


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)