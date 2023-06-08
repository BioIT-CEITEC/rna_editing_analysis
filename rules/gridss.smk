

rule gridss_computation:
    input:  bams = expand("mapped/{sample_name}.bam",sample_name=sample_tab.sort_values(by=['tumor_normal']).sample_name.tolist()),
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            ref_fai=expand("{ref_dir}/seq/{ref_name}.fa.fai",ref_dir=reference_directory,ref_name=config["reference"])[0],
            bwa_index = expand("{ref_dir}/index/BWA/{ref}.{suffix}", ref_dir=reference_directory, ref=config["reference"], suffix=["amb", "ann", "bwt", "pac", "sa"])
    output: vcf="structural_varcalls/all_samples/gridss/all_sample_variants.vcf"
    params: tumor_normal_paired = config["tumor_normal_paired"]
    log:    "logs/all_samples/gridss/gridss_computation.log",
    threads: 16
    conda:  "../wrappers/gridss/env.yaml"
    script: "../wrappers/gridss/script.py"


rule gridss_get_per_sample_res:
    input:  all_sample_vcf="structural_varcalls/all_samples/gridss/all_sample_variants.vcf"
    output: vcf="structural_varcalls/{sample_name}/gridss/result_SV.vcf",
    log:    "logs/{sample_name}/gridss/get_per_sample_res.log",
    shell:
        "touch {output.vcf}"
