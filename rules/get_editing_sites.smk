rule get_editing_sites:
     input:
        bam = "mapped/{sample}.bam",
        gtf = expand("{ref_dir}/annot/{ref}.gtf", ref_dir = reference_directory, ref = config["reference"])[0],
        known_editing = expand("{ref_dir}/other/known_editing_sites/{ref}.csv", ref_dir = reference_directory, ref = config["reference"])[0],
        snp = expand("{ref_dir}/other/snp/{ref}.snp.bed", ref_dir = reference_directory, ref = config["reference"])[0]
     params:
        chrm = "{chrom}"
     output:
        fa = "per_sample_results/{sample}.seq.{chrom}.fa",
        editing_sites = "per_sample_results/{sample}.editing_sites.{chrom}.tsv"
     log:
        run = "sample_logs/{sample}/get_editing_sites.{chrom}.log"
     threads: 20
     conda: "../wrappers/get_editing_sites/env.yaml"
     script: "../wrappers/get_editing_sites/script.py"
