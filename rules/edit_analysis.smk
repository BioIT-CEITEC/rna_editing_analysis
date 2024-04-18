#
# def get_bam_input(wildcards):
#     if config["calling_type"] == "tumor_normal":
#         input_bam_name = sample_tab.loc[(sample_tab["tumor_normal"] == wildcards.tumor_normal) & (sample_tab["donor"]==wildcards.sample_name), "sample_name"]
#         return "mapped/" + input_bam_name + ".bam"
#     else:
#         return expand("mapped/{input_bam}.bam",input_bam=wildcards.sample_name)[0]
#
# def get_region_bed_input(wildcards):
#     if config["lib_ROI"] != "wgs":
#         return expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
#     else:
#         return expand("structural_varcalls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0]
#
# rule jabCoNtool_per_sample_coverage:
#     input:  bam= get_bam_input,
#             region_bed = get_region_bed_input,
#             ref_dict= expand("{ref_dir}/seq/{ref_name}.fa.fai",ref_dir=reference_directory,ref_name=config["reference"])[0],
#     output: cov_tab = "structural_varcalls/{sample_name}/jabCoNtool/{tumor_normal}.region_coverage.tsv",
#     log:    "logs/{sample_name}/jabCoNtool/{tumor_normal}_get_coverage.log"
#     threads: 8
#     resources: mem=10
#     conda:  "../wrappers/jabCoNtool/per_sample_coverage_computing/env.yaml"
#     script: "../wrappers/jabCoNtool/per_sample_coverage_computing/script.py"



# rule editing_site_per_sample_snp_AF:
#     input:  bam = "mapped/{sample_name}.bam",
#             ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#             snp_tsv = expand("{ref_dir}/other/known_editing_sites/{ref_name}.csv",ref_dir=reference_directory,ref_name=config["reference"])[0],
#     output: snp_tab = "editing_analysis/known_editing_sites/{sample_name}.editing_sites_AF.tsv",
#     log:    "logs/{sample_name}/editing_sites_get_AF.log"
#     threads: 8
#     resources: mem=10
#     conda:  "../wrappers/per_sample_snp_AF_computing/env.yaml"
#     script: "../wrappers/per_sample_snp_AF_computing/script.py"


rule editing_site_per_sample_snp_AF:
    input:  bam = "mapped/{sample_name}.bam",
            ref = config["organism_fasta"],
            snp_tsv = "editing_analysis/potential_edit_sites.tsv"
    output: snp_tab = "editing_analysis/potential_edit_site_AFs/{sample_name}.editing_sites_AF.tsv",
    log:    "logs/{sample_name}/editing_sites_get_AF.log"
    threads: 8
    resources: mem=10
    conda:  "../wrappers/editing_site_per_sample_snp_AF/env.yaml"
    script: "../wrappers/editing_site_per_sample_snp_AF/script.py"


rule editing_site_annotation:
    input:  tsv_for_vep = "editing_analysis/potential_edit_sites.tsv"
    output: annotated = "editing_analysis/potential_edit_sites.annotated.tsv"
    log:    "logs/potential_edit_site_annotation.log"
    threads: 20
    resources:
        mem_mb=8000
    params: ref = config["organism_fasta"],
            vep_dir = config["organism_vep_dir"],
            organism_name = config["organism"]
    conda:  "../wrappers/editing_site_annotation/env.yaml"
    script: "../wrappers/editing_site_annotation/script.py"


rule editing_res_processing:
    input:  all_vars_tsv = expand("editing_analysis/known_editing_sites/{sample_name}.known_editing_sites_AF.tsv",sample_name=sample_tab.sample_name.tolist()),
            annotated = "editing_analysis/potential_edit_sites.annotated.tsv"
    output: res_tab = "editing_analysis/res_tab.tsv"
    log:    "logs/editing_res_processing.log"
    threads: 8
    resources: mem=10
    conda:  "../wrappers/editing_res_processing/env.yaml"
    script: "../wrappers/editing_res_processing/script.py"


rule final_editing_report:
    input:  res_tab = "editing_analysis/res_tab.tsv"
    output: html = "reports/editing_analysis_report.html"
    # params: sample_name = sample_tab.sample_name,
    #         config = "./config.json"
    # conda: "../wrappers/final_alignment_report/env.yaml"
    # script: "../wrappers/final_alignment_report/script.Rmd"
    shell:
        "mkdir -p reports; touch {output.html}"

#
# def jabCoNtool_cnv_computation_inputs(wildcards):
#     input_dict = {}
#     if config["calling_type"] == "tumor_normal":
#         input_dict["normal_sample_cov"] = set(expand("structural_varcalls/{sample_name}/jabCoNtool/normal.region_coverage.tsv",sample_name=sample_tab.loc[
#             sample_tab.tumor_normal == "normal", "donor"].tolist()))
#         input_dict["sample_cov"] = set(expand("structural_varcalls/{sample_name}/jabCoNtool/tumor.region_coverage.tsv",sample_name=
#             sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist()))
#         if config["jabCoNtool_use_snps"] == True:
#             input_dict["snp_AF"] = set(expand("structural_varcalls/{sample_name}/jabCoNtool/tumor.snpAF.tsv",sample_name=
#                 sample_tab.loc[sample_tab.tumor_normal == "tumor", "donor"].tolist()))
#             input_dict["normal_snp_AF"] = set(expand("structural_varcalls/{sample_name}/jabCoNtool/normal.snpAF.tsv",sample_name=
#                 sample_tab.sample_name.tolist()))
#     else:
#         input_dict["sample_cov"] = set(expand("structural_varcalls/{sample_name}/jabCoNtool/sample.region_coverage.tsv",sample_name=sample_tab.sample_name.tolist()))
#         if config["jabCoNtool_use_snps"] == True:
#             input_dict["snp_AF"] = set(expand("structural_varcalls/{sample_name}/jabCoNtool/sample.snpAF.tsv",sample_name=sample_tab.sample_name.tolist()))
#     if config["use_cohort_data"] == True:
#         input_dict["cohort_data"] = "cohort_data/cohort_data/jabCoNtool/cohort_info_tab.tsv"
#     if config["lib_ROI"] == "wgs":
#         input_dict["region_bed"] = expand("structural_varcalls/all_samples/binned_genome_{window_size}.bed",window_size=config["wgs_bin_size"])[0]
#         if config["jabCoNtool_normalize_to_GC"] == True:
#             input_dict["GC_profile_file"] = expand("structural_varcalls/all_samples/GC_profile_{window_size}.cnp",window_size=config["wgs_bin_size"])[0]
#         if config["jabCoNtool_remove_centromeres"] == True:
#             input_dict["cytoband_file"] = expand("{ref_dir}/other/cytoband/{ref}.cytoband.tsv",ref_dir=reference_directory,ref = config["reference"])[0]
#     else:
#         input_dict["region_bed"] = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
#     if config["jabCoNtool_use_snps"] == True:
#         input_dict["snp_bed"] = expand("{ref_dir}/other/snp/{lib_ROI}/{lib_ROI}_snps.tsv",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
#     return input_dict
#
#
# rule jabCoNtool_cnv_computation:
#     input: unpack(jabCoNtool_cnv_computation_inputs)
#     output: all_res_prob_tab="structural_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv",
#             cohort_info_tab="structural_varcalls/all_samples/jabCoNtool/cohort_info_tab.tsv"
#     params: jabCoNtool_predict_TL = config["jabCoNtool_predict_TL"],
#             calling_type = config["calling_type"],
#             lib_ROI= config["lib_ROI"],
#             max_CNV_occurance_in_cohort = config["max_CNV_occurance_in_cohort"]
#     log:    "logs/all_samples/jabCoNtool/cnv_computation.log",
#     threads: workflow.cores
#     conda:  "../wrappers/jabCoNtool/cnv_computation/env.yaml"
#     script: "../wrappers/jabCoNtool/cnv_computation/script.py"
#
#     rulefinal_alignment_report:
#     input: all_vars_tsv="final_CNV_results/CNV_variants.tsv",
#




# rule jabCoNtool_get_per_sample_res:
#     input:  all_res_prob_tab="structural_varcalls/all_samples/jabCoNtool/final_CNV_probs.tsv"
#     output: CNV_res="structural_varcalls/{sample_name}/jabCoNtool/CNV_varcalls.tsv",
#     log:    "logs/{sample_name}/jabCoNtool/get_per_sample_res.log",
#     conda:  "../wrappers/jabCoNtool/get_per_sample_res/env.yaml"
#     script: "../wrappers/jabCoNtool/get_per_sample_res/script.py"