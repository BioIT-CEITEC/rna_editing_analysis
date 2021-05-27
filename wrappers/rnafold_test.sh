#RNAFold test
INPUT=/mnt/ssd/ssd_3/temp/ailar/seq.fa
SCRATCH=/mnt/ssd/ssd_3/temp/


source activate rna_edit
cp ${INPUT} ${SCRATCH}
cd ${SCRATCH}

############
#analyse
RNAfold --MEA -d2 < seq.fa
