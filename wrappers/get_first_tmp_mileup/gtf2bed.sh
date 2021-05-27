#gtf2bed
INPUT_DIR=/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/annot
GTF=GRCh38-p10.gtf
SCRATCH=/mnt/ssd/ssd_3/temp/vasek/gtf2bed


source activate rnaseq

mkdir -p ${SCRATCH}
cp ${INPUT_DIR}/${GTF} ${SCRATCH}
cd ${SCRATCH}

awk '{ if ($0 ~ "transcript_id") ; else print $0" transcript_id \"\";"; }' ${GTF} | gtf2bed > ${GTF%.gtf}_whole_genes.bed

cp ${GTF%.gtf}_whole_genes.bed ${INPUT_DIR}
