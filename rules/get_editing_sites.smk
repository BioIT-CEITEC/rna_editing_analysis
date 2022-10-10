rule get_editing_sites:
     input:
        bam = "/mnt/nfs/shared/RNA_and_Immunity/sequencing_results/projects/M6A_KD_flies_editing/RNA_alignment/mapped/{sample}.bam",
        gtf = expand("{ref_dir}/annot/{rel}/{ref}.gtf",ref_dir = reference_directory,rel = config["release"],ref = config["reference"])[0],
     params:
        chr = config["chromosome"]
     output:
        mismatch_tab = "results/{sample}.mismatch_tab.tsv",
     log:
        run = "logs/{sample}/get_editing_sites.log"
     threads: 20
     conda: "../wrappers/get_editing_sites/env.yaml"
     script: "../wrappers/get_editing_sites/script.py"


