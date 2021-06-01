rule get_editing_sites:
     input:
       bam = INPUTS_DIR + "/mapped/{full_name}.bam",
       gtf = lambda wildcards:  expand("{dir}/{organism}/{ref}/annot/{ref}.gtf",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"reference"].min())[0],
       known_editing = lambda wildcards:  expand("{dir}/{organism}/{ref}/other/known_editing_sites/Human_AG_all_hg38_v2.csv",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"reference"].min())[0],
       snp = lambda wildcards:  expand("{dir}/{organism}/{ref}/other/snp/{ref}.snp.bed",dir = REF_DIR,organism = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"organism"].min(),ref = cfg.loc[cfg[SAMPLE] == wildcards.full_name,"reference"].min())[0],
     output:
       fa = ADIR + "/per_sample_results/{full_name}.seq.fa",
       editing_sites = ADIR + "/per_sample_results/{full_name}.editing_sites.tsv",
     log:
       run = ADIR + "/sample_logs/{full_name}/get_editing_sites.log"
     threads: 20
     conda: "../wraps/RNA_editing/get_editing_sites/env.yaml"
     script: "../wraps/RNA_editing/get_editing_sites/script.py"
