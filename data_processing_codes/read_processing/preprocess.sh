#!/bin/bash

set -x
TCGA_PROJECT=$1
MANIFEST_FILE_PATH=$2

echo TCGA_PROJECT: ${TCGA_PROJECT}
echo MANIFEST_FILE_PATH: ${MANIFEST_FILE_PATH}

cd /TCGA-BRCA/${TCGA_PROJECT}/


# bam to unmapped_bam
mkdir -p unmapped_bam
echo "===Extract unmapped==="
time python /TCGA-BRCA/code/bamTOunmapped_manifest_20180802.py ${MANIFEST_FILE_PATH} \
/TCGA-BRCA/${TCGA_PROJECT}/bam/ \
/TCGA-BRCA/${TCGA_PROJECT}/unmapped_bam/ >bam2unmappedBam_log

# unmapped_bam to fastq
mkdir -p unmapped_fastq
echo "===Turn unmapped bam to fastq==="
time python /TCGA-BRCA/code/bamTofastq_remainUnmapped_manifest_20180802.py ${MANIFEST_FILE_PATH} \
/TCGA-BRCA/${TCGA_PROJECT}/unmapped_bam/ \
/TCGA-BRCA/${TCGA_PROJECT}/unmapped_fastq/ >unmappedBam2fastq_log

# sRNAnalyzer processing
cd unmapped_fastq
echo "===sRNAnalyzer preprocessing==="
~/cancerMicrobiome/sRNAnalyzer/preprocess.pl --config /TCGA-BRCA/TCGA_sRNAnalyzer.yaml 2>&1 | tee /TCGA-BRCA/${TCGA_PROJECT}/preprocessing_log
