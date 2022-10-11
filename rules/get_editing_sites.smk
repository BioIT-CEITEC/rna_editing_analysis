rule get_editing_sites:
     input:
        bam = "/mnt/nfs/shared/RNA_and_Immunity/sequencing_results/projects/M6A_KD_flies_editing/RNA_alignment/mapped/{sample}.bam",
        gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir = reference_directory,ref = config["reference"])[0],
     output:
        mismatch_tab = "results/{sample}.mismatch_cov_tab.tsv",
     params:
        prefix = "{sample}"
     log:
        run = "logs/{sample}/get_editing_sites.log"
     threads: 20
     conda: "../wrappers/get_editing_sites/env.yaml"
     script: "../wrappers/get_editing_sites/script.py"


