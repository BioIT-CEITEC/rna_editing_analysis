# rule get_editing_sites:
#      input:
#         bam = "/mnt/nfs/shared/RNA_and_Immunity/sequencing_results/projects/M6A_KD_flies_editing/RNA_alignment/mapped/{sample}.bam",
#         gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir = reference_directory,ref = config["reference"])[0],
#      output:
#         mismatch_tab = "results/{sample}.mismatch_cov_tab.tsv",
#      params:
#         prefix = "{sample}"
#      log:
#         run = "logs/{sample}/get_editing_sites.log"
#      threads: 20
#      conda: "../wrappers/get_editing_sites/env.yaml"
#      script: "../wrappers/get_editing_sites/script.py"

rule get_editing_sites:
     input:
        bam = BR.remote("mapped/{sample}.bam"),
        bai = BR.remote("mapped/{sample}.bam.bai"),
        gtf = BR.remote(expand("{ref_dir}/annot/{ref}.gtf",ref_dir = reference_directory,ref = config["reference"])),
     output:
        mismatch_tab = BR.remote("results/{sample}.mismatch_cov_tab.tsv"),
        bed = BR.remote("results/{sample}.positions_tmp.bed")
     log:
        run = BR.remote("logs/{sample}/get_editing_sites.log")
     threads: 20
     conda: "../wrappers/get_editing_sites/env.yaml"
     script: "../wrappers/get_editing_sites/script.py"


#
# rule read_single_chr_bam_as_sam:
#     input:
#         bam = "mapped/{sample}.bam",
#         bed = "results/{sample}.positions.bed",
#     output:
#         cov_tab="results/{sample}.cov_tab.tsv",
#     params:
#         prefix="{sample}"
#     log:
#         run="logs/{sample}/mismatch_coverage_count.log"
#     threads: 20
#     conda: "../wrappers/mismatch_coverage_count/env.yaml"
#     script: "../wrappers/mismatch_coverage_count/script.py"
#
#
# rule create_mismatch_table:
#     input:
#         bam = "mapped/{sample}.bam",
#         bed = "results/{sample}.positions.bed",
#     output:
#         cov_tab="results/{sample}.cov_tab.tsv",
#     params:
#         prefix="{sample}"
#     log:
#         run="logs/{sample}/mismatch_coverage_count.log"
#     threads: 20
#     conda: "../wrappers/mismatch_coverage_count/env.yaml"
#     script: "../wrappers/mismatch_coverage_count/script.py"
#
#
# rule coverage_count:
#     input:
#         bam = "mapped/{sample}.bam",
#         bed = "results/{sample}.positions.bed",
#     output:
#         cov_tab="results/{sample}.cov_tab.tsv",
#     params:
#         prefix="{sample}"
#     log:
#         run="logs/{sample}/mismatch_coverage_count.log"
#     threads: 20
#     conda: "../wrappers/mismatch_coverage_count/env.yaml"
#     script: "../wrappers/mismatch_coverage_count/script.py"
#
#
# rule mismatch_coverage_table:
#     input:
#         bam = "mapped/{sample}.bam",
#         bed = "results/{sample}.positions.bed",
#     output:
#         cov_tab="results/{sample}.cov_tab.tsv",
#     params:
#         prefix="{sample}"
#     log:
#         run="logs/{sample}/mismatch_coverage_count.log"
#     threads: 20
#     conda: "../wrappers/mismatch_coverage_count/env.yaml"
#     script: "../wrappers/mismatch_coverage_count/script.py"