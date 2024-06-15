rule bam_preprocessing:
    input:  bam = "mapped/{sample}.bam",
            ref = config["organism_fasta"]
    output: bam = "mapped/{sample}.MD.bam",
    log:    "logs/{sample}/bam_preprocessing.log"
    threads: 30
    resources: mem=10
    conda:  "../wrappers/bam_preprocessing/env.yaml"
    script: "../wrappers/bam_preprocessing/script.py"

rule jacusa_call:
    input:  bam = "mapped/{sample}.MD.bam",
            ref = config["organism_fasta"]
    output: jacusa = "jacusa_call/{sample}.txt",
    log:    "logs/{sample}/jacusa_calling.log"
    threads: 30
    params: jacusa = config["tooldir"] + "/JACUSA2/JACUSA_v2.0.4.jar",
            strand = config["strand"],
            min_mapq = config["min_mapq"],
            filternh = config["filternh"],
            filternm = config["filternm"]
    conda:  "../wrappers/jacusa_call/env.yaml"
    script: "../wrappers/jacusa_call/script.py"

rule jacusa_helper:
    input:  raw = "jacusa_call/{sample}.txt"
    output: processed = "jacusa_call/{sample}_potential_edit_sites.tsv"
    log:    "logs/{sample}/jacusa_helper.log"
    threads: 30
    conda:  "../wrappers/jacusa_helper/env.yaml"
    script: "../wrappers/jacusa_helper/script.py"