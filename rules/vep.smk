

rule filter_vep:
    input:
        vcf="snv_indels/ensemble_vcf/{sample_name}.ensembled.vep_annotated.vcf", #vstup je ENSEMBLE VCF z nekolika indel calleru
    output:
        vcf="cnv_sv/filter_vep_germline_vcf/{sample_name}.germline.vcf",
    params:
        filter='--filter "DP > 50" --filter "MAX_AF >= 0.001"',
    log:
        "logs/filter_vep/{sample_name}.germline.vcf.log",
    threads: 8
    conda:
        "../wrappers/vep/env.yaml"
    shell:
        "(filter_vep -i {input.vcf} -o {output.vcf} {params.filter}) &> {log}"
