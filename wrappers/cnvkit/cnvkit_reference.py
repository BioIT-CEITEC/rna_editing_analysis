######################################
# wrapper for rule: germline_strelka
######################################
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: cnv reference \n##\n")
f.close()

shell.executable("/bin/bash")


if snakemake.input.normal_coverage_inputs:
    command = "cnvkit.py reference " + " ".join(snakemake.input.normal_coverage_inputs) + \
              " --fasta " + snakemake.input.reference + \
              " -o " + snakemake.output.reference_cnn + \
              " >> " + log_filename + " 2>&1"
else:
    command = "cnvkit.py reference " + \
              " --fasta " + snakemake.input.reference + \
              " -o " + snakemake.output.reference_cnn + \
              " -t " + snakemake.input.target + \
              " -a " + snakemake.input.antitarget + \
              " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
