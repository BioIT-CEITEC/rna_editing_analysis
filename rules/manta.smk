# def normal_bam_input(wildcards):
#     return expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "normal") & (sample_tab["donor"]==wildcards.donor), "sample_name"])
#
# def normal_bam_bai_input(wildcards):
#     return expand("mapped/{input_bam}.bam.bai",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "normal") & (sample_tab["donor"]==wildcards.donor), "sample_name"])
#
# def tumor_bam_input(wildcards):
#     return expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "tumor") & (sample_tab["donor"]==wildcards.donor), "sample_name"])
#
# def tumor_bam_bai_input(wildcards):
#     return expand("mapped/{input_bam}.bam.bai",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "tumor") & (sample_tab["donor"]==wildcards.donor), "sample_name"])

def bam_inputs(wildcards):
    if config["calling_type"] == "tumor_normal":
        return {'tumor': expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "tumor") & (sample_tab["donor"]==wildcards.sample_name), "sample_name"])[0],
                'normal':expand("mapped/{input_bam}.bam",input_bam=sample_tab.loc[(sample_tab["tumor_normal"] == "normal") & (sample_tab["donor"]==wildcards.sample_name), "sample_name"])[0]}
    else:
        return {'tumor': expand("mapped/{tumor_bam}.bam",tumor_bam=wildcards.sample_name)[0]}

    # if config["material"] != "RNA":
    #     tag = "bam"
    # else:
    #     tag = "RNAsplit.bam"



rule manta:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions_gz=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed.gz",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
        regions_tbi=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed.gz.tbi",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output: vcf="structural_varcalls/{sample_name}/manta/result_SV.vcf",
    log: "logs/{sample_name}/manta/manta.log"
    threads: 5
    resources: mem=6
    params: dir = "structural_varcalls/{sample_name}/manta",
            manta_sv_vcf="structural_varcalls/{sample_name}/manta/results/variants/tumorSV.vcf.gz",
            library_scope = config["lib_ROI"],
            calling_type = config["calling_type"]
    conda:  "../wrappers/manta/env.yaml"
    script: "../wrappers/manta/script.py"

# rule manta:
#     input:
#         unpack(bam_inputs),
#         ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#     output:
#         cand_si_vcf="cnv_sv/manta/{donor}/results/variants/candidateSmallIndels.vcf.gz",
#         cand_si_tbi="cnv_sv/manta/{donor}/results/variants/candidateSmallIndels.vcf.gz.tbi",
#         cand_sv_vcf="cnv_sv/manta/{donor}/results/variants/candidateSV.vcf.gz",
#         cand_sv_tbi="cnv_sv/manta/{donor}/results/variants/candidateSV.vcf.gz.tbi",
#         dipl_sv_vcf="cnv_sv/manta/{donor}/results/variants/diploidSV.vcf.gz",
#         dipl_sv_tbi="cnv_sv/manta/{donor}/results/variants/diploidSV.vcf.gz.tbi",
#         som_sv_vcf="cnv_sv/manta/{donor}/results/variants/somaticSV.vcf.gz",
#         som_sv_tbi="cnv_sv/manta/{donor}/results/variants/somaticSV.vcf.gz.tbi",
#     log:
#         "logs/manta/{donor}/run_workflow.py.log",
#     threads: 8
#     conda:  "../wrappers/manta/env.yaml"
#     script: "../wrappers/manta/script.py"
    # shell:
    #     "configManta.py "
    #     "--tumorBam={input.bam_t} "
    #     "--normalBam={input.bam_n} "
    #     "--referenceFasta={input.ref} "
    #     "--runDir=cnv_sv/manta/{wildcards.donor} &> {log}"

# rule manta_run_workflow_tn:
#     input:
#         bam_t=tumor_bam_input,
#         bai_t=tumor_bam_bai_input,
#         bam_n=normal_bam_input,
#         bai_n=normal_bam_bai_input,
#         ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#         scrpt="cnv_sv/manta/{donor}/runWorkflow.py",
#     output:
#
#         wrk_dir=directory("cnv_sv/manta/{donor}/workspace"),
#     log:
#         "logs/manta/{donor}/manta_tn.log",
#     threads: 8
#     conda:
#         "../wrappers/manta/env.yaml"
#     shell:
#         "{input.scrpt} "
#         "-j {threads} "
#         "-g unlimited &> {log}"
