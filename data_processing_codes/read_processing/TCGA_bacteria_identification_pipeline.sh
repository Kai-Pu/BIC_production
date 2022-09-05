#!/bin/bash

WK_DIR=$1
UUID=$2

export PATH=$PATH:/work1/chialang1220/reference/kpc/bowtie-1.2.1.1
SRNANALYZER_DIR=/work1/chialang1220/reference/kpc/sRNAnalyzer/
SRNANALYZER_YAML_FILE=/work1/chialang1220/reference/kpc/TCGA_sRNAnalyzer.yaml
SRNANALYZER_DB_CONFIG_FILE=/work1/chialang1220/reference/kpc/TCGA_sRNAnalyzer.config
CODE_DIR=/work1/chialang1220/reference/kpc/pipeline_code/
LOG_FILE=${UUID}_exec_log

errorCheck() {
	if [ $1 -ne 0 ];then
		echo "Error occurring in file $2 with step $3"
		exit 1
	else
		echo "Successfully executed in file $2 with step $3"
	fi
}


echo "SRNANALYZER_DIR: ${SRNANALYZER_DIR}" 2>&1 | tee ${LOG_FILE}
echo "SRNANALYZER_YAML_FILE: ${SRNANALYZER_YAML_FILE}" 2>&1 | tee -a ${LOG_FILE}
echo "SRNANALYZER_DB_CONFIG_FILE: ${SRNANALYZER_DB_CONFIG_FILE}" 2>&1 | tee -a ${LOG_FILE}
echo "CODE_DIR: ${CODE_DIR}" 2>&1 | tee -a ${LOG_FILE}
echo "WORKING_DIRECTORY: ${WK_DIR}" 2>&1 | tee -a ${LOG_FILE}
echo "UUID: ${UUID}" 2>&1 | tee -a ${LOG_FILE}
echo "LOG_FILE: ${LOG_FILE}" 2>&1 | tee -a ${LOG_FILE}

echo "Enter processing file directory: ${WK_DIR}/${UUID}" 2>&1 | tee -a ${LOG_FILE}
cd ${WK_DIR}/${UUID}

#echo "sRNAnalyzer preprocessing start" 2>&1 | tee -a ${LOG_FILE}
#echo "${SRNANALYZER_DIR}/preprocess.pl --config ${SRNANALYZER_YAML_FILE}" 2>&1 | tee -a ${LOG_FILE}
#${SRNANALYZER_DIR}/preprocess.pl --config ${SRNANALYZER_YAML_FILE} 2>&1 | tee -a ${LOG_FILE}
#errorCheck "$?" "${UUID}" "preprocessing" 2>&1 | tee -a ${LOG_FILE}

echo 2>&1 | tee -a ${LOG_FILE}
echo "===sRNAnalyzer aligning start===" 2>&1 | tee -a ${LOG_FILE}
echo "${SRNANALYZER_DIR}/align.pl ${WK_DIR}/${UUID} ${SRNANALYZER_YAML_FILE} ${SRNANALYZER_DB_CONFIG_FILE}" 2>&1 | tee -a ${LOG_FILE}
${SRNANALYZER_DIR}/align.pl ${WK_DIR}/${UUID} ${SRNANALYZER_YAML_FILE} ${SRNANALYZER_DB_CONFIG_FILE} 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "alinging" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "===adding header line===" 2>&1 | tee -a ${LOG_FILE}
echo "sed -i '1i \\matchType.strand\treadID\treadSeq\tmatchID\toffsetq\tmatchNum' ${UUID}_Processed.profile" 2>&1 | tee -a ${LOG_FILE}
sed -i '1i \\matchType.strand\treadID\treadSeq\tmatchID\toffsetq\tmatchNum' ${UUID}_Processed.profile
errorCheck "$?" "${UUID}" "adding_header" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "===sRNAnalyzer adding description information start===" 2>&1 | tee -a ${LOG_FILE}
echo "${SRNANALYZER_DIR}/selfDev/desProfile.pl ${UUID}_Processed.profile ${SRNANALYZER_DB_CONFIG_FILE}" 2>&1 | tee -a ${LOG_FILE}
${SRNANALYZER_DIR}/selfDev/desProfile.pl ${UUID}_Processed.profile ${SRNANALYZER_DB_CONFIG_FILE} 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "adding_description" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "===sRNAnalyzer adding taxa information===" 2>&1 | tee -a ${LOG_FILE}
echo "${SRNANALYZER_DIR}/selfDev/taxProfile.pl ${UUID}_Processed_des.profile" 2>&1 | tee -a ${LOG_FILE}
${SRNANALYZER_DIR}/selfDev/taxProfile.pl ${UUID}_Processed_des.profile 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "adding_taxa_information" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "===build OTU table and find seq mapped bacteria number===" 2>&1 | tee -a ${LOG_FILE}
echo "${CODE_DIR}/buildOTUtable_findSeqFreq_20200221_v2.py ${UUID}" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/buildOTUtable_findSeqFreq_20200221_v2.py ${UUID} 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "building_OTU_and_finding_seq_mapped_bacteria_number" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "===finding unique seq matched bacteria===" 2>&1 | tee -a ${LOG_FILE}
echo "${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Species" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Species 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "finding_unique_seq_bacteria_Species" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Genus" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Genus 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "finding_unique_seq_bacteria_Genus" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Family" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Family 2>&1 | tee -a ${LOG_FILE} 
errorCheck "$?" "${UUID}" "finding_unique_seq_bacteria_Family" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Order" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Order 2>&1 | tee -a ${LOG_FILE} 
errorCheck "$?" "${UUID}" "finding_unique_seq_bacteria_Order" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Class" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Class 2>&1 | tee -a ${LOG_FILE} 
errorCheck "$?" "${UUID}" "finding_unique_seq_bacteria_Class" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Phylum" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/findBothBacteriaNumAndMatchedSeq_for_databases_20200223.py ${UUID} Phylum 2>&1 | tee -a ${LOG_FILE} 
errorCheck "$?" "${UUID}" "finding_unique_seq_bacteria_Phylum" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "${CODE_DIR}/findBacteriaWithUniqueSeq_20200223.py ${UUID}" 2>&1 | tee -a ${LOG_FILE}
${CODE_DIR}/findBacteriaWithUniqueSeq_20200223.py ${UUID} 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "finding_unique_seq_bacteria_all_level" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "===remove intermediate files===" 2>&1 | tee -a ${LOG_FILE}
#echo "rm ${UUID}_Processed.fa" 2>&1 | tee -a ${LOG_FILE}
#rm ${UUID}_Processed.fa 2>&1 | tee -a ${LOG_FILE}
#errorCheck "$?" "${UUID}" "removing_Processed.fa" 2>&1 | tee -a ${LOG_FILE}

echo "rm ${UUID}_Processed.feature" 2>&1 | tee -a ${LOG_FILE}
rm ${UUID}_Processed.feature 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "removing_Processed.feature" 2>&1 | tee -a ${LOG_FILE}

echo "rm ${UUID}_Processed.profile" 2>&1 | tee -a ${LOG_FILE}
rm ${UUID}_Processed.profile 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "removing_Processed.profile" 2>&1 | tee -a ${LOG_FILE}

echo "rm ${UUID}_Processed_des.profile" 2>&1 | tee -a ${LOG_FILE}
rm ${UUID}_Processed_des.profile 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "removing_Processed_des.profile" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "===gzip the annot.profile===" 2>&1 | tee -a ${LOG_FILE}
echo "gzip ${UUID}_Processed_anno.profile" 2>&1 | tee -a ${LOG_FILE}
gzip ${UUID}_Processed_anno.profile 2>&1 | tee -a ${LOG_FILE}
errorCheck "$?" "${UUID}" "gzip_Processed_anno.profile" 2>&1 | tee -a ${LOG_FILE}
echo 2>&1 | tee -a ${LOG_FILE}

echo "===Finished processes===" 2>&1 | tee -a ${LOG_FILE}
