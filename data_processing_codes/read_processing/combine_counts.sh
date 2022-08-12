#!/bin/bash

TCGA_PROJECT=$1
WK_DIR=$2
UUID_LIST=$3
CODE_DIR=/work1/chialang1220/reference/kpc/pipeline_code/
LOG_FILE=${WK_DIR}/${TCGA_PROJECT}_combine_count_log

errorCheck() {
	if [ $1 -ne 0 ];then
		echo "Error occurring in taxalevel $2 with step $3"
		exit 1
	else
		echo "Successfully executed in taxalevel $2 with step $3"
	fi
}

echo "TCGA-PROJECT: ${TCGA_PROJECT}" 2>&1 | tee -a ${LOG_FILE}
echo "WORKING_DIRECTORY: ${WK_DIR}" 2>&1 | tee -a ${LOG_FILE}
echo "UUID_LIST: ${UUID_LIST}" 2>&1 | tee -a ${LOG_FILE}
echo "CODE_DIR: ${CODE_DIR}" 2>&1 | tee -a ${LOG_FILE}

cd ${WK_DIR}

echo "===Combine count===" 2>&1 | tee -a ${LOG_FILE}
echo "${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Species" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Species 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "Species" "combining_count" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE} 

echo "${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Genus" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Genus 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "Genus" "combining_count" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE} 

echo "${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Family" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Family 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "Family" "combining_count" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE} 

echo "${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Order" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Order 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "Order" "combining_count" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE} 

echo "${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Class" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Class 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "Class" "combining_count" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE} 

echo "${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Phylum" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/combineSamplesBacteriaCounts.py ${TCGA_PROJECT} ${WK_DIR} ${UUID_LIST} Phylum 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "Phylum" "combining_count" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE} 

echo "===Comine successfully===" 2>&1 | tee -a ${LOG_FILE}
