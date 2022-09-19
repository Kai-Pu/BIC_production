########################################################################
## usage: python3 find_readnum_of_bam_unmappedFastq_processedFa.py 
## aim: find reads info
## total_reads_tcga_bam, total_reads_tcga_bam_human_mapped, total_reads_human_unmapped_fastq, total_reads_srnanalyzer_processed_fa
## samtools version: 1.3.1
## author: Kai-Pu Chen
## date: 20220907
######################################################

import os

## test: subset => 5 GBM samples
# barcode_uuid_matched_file_path = "/home/kaipu/microbiome_TCGA_miRNA/FFPE_nonFFPE/nonFFPE_PT_STN/TCGA_miRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal_GBM.txt"
# total_reads_tcga_bam_output_path = "/TCGA-BRCA/reads_info/total_reads_tcga_bam_GBM.txt"
# total_reads_tcga_bam_human_mapped_output_path = "/TCGA-BRCA/reads_info/total_reads_tcga_bam_human_mapped_GBM.txt"
# total_reads_human_unmapped_fastq_output_path = "/TCGA-BRCA/reads_info/total_reads_human_unmapped_fastq_GBM.txt"
# total_reads_srnanalyzer_processed_fa_output_path = "/TCGA-BRCA/reads_info/total_reads_srnanalyzer_processed_fa_GBM.txt"

barcode_uuid_matched_file_path = "/home/kaipu/microbiome_TCGA_miRNA/FFPE_nonFFPE/nonFFPE_PT_STN/TCGA_miRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal.txt"
total_reads_tcga_bam_output_path = "/TCGA-BRCA/reads_info/total_reads_tcga_bam.txt"
total_reads_tcga_bam_human_mapped_output_path = "/TCGA-BRCA/reads_info/total_reads_tcga_bam_human_mapped.txt"
total_reads_human_unmapped_fastq_output_path = "/TCGA-BRCA/reads_info/total_reads_human_unmapped_fastq.txt"
total_reads_srnanalyzer_processed_fa_output_path = "/TCGA-BRCA/reads_info/total_reads_srnanalyzer_processed_fa.txt"


print("barcode_uuid_matched_file_path: ", barcode_uuid_matched_file_path, "\n")



barcode_uuid_match_table = []

with open(barcode_uuid_matched_file_path, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        # print(line)
        barcode_uuid_match_table.append(line)


error_time_count = 0

for i in range(1, len(barcode_uuid_match_table)):
    # Patient_barcode Project_name    Sample_barcode  Aliquot_barcode miRNAseq_file_uuid      Sample_type     Filtered_Aliquot_barcode
    # TCGA-06-0678    GBM     TCGA-06-0678-11A        TCGA-06-0678-11A-32R-A36C-13    e45a44b3-a124-4a4c-9fec-3991740ada26    Solid Tissue Normal
    # TCGA-06-0675    GBM     TCGA-06-0675-11A        TCGA-06-0675-11A-32R-A36C-13    7f9878f5-e606-46e9-8b1a-de9e7a617575    Solid Tissue Normal
    project_name = barcode_uuid_match_table[i][1]
    sample_uuid = barcode_uuid_match_table[i][4]
    # /TCGA-BRCA/GBM/bam/e45a44b3-a124-4a4c-9fec-3991740ada26/TCGA-06-0678-11A-32R-A36C-13_mirna_gdc_realn.bam 
    bam_path = "/TCGA-BRCA/" + project_name + "/bam/" + sample_uuid + "/*.bam"
    # /TCGA-BRCA/GBM/unmapped_fastq/e45a44b3-a124-4a4c-9fec-3991740ada26.fastq
    human_unmaped_fastq_path = "/TCGA-BRCA/" + project_name + "/unmapped_fastq/" + sample_uuid + ".fastq"
    # /TCGA-BRCA/GBM/unmapped_fastq/e45a44b3-a124-4a4c-9fec-3991740ada26_Processed.fa
    srnanalyzer_processed_fa_path = "/TCGA-BRCA/" + project_name + "/unmapped_fastq/" + sample_uuid + "_Processed.fa"


    ## 1.find and save bam flagstat
    samtools_flagstat_output_path = "/TCGA-BRCA/" + project_name + "/bam/" + sample_uuid + "/samtools_flagstat.txt"

    cmd = "samtools flagstat " + bam_path + " >" + samtools_flagstat_output_path
    print("===== " + str(i) + " =====")
    # print(cmd)

    cmdStateCode = os.system(cmd)
    stderrMsgSTATE_OK = cmd + ' => Done!'
    stderrMsgSTATE_FAIL = "Error in processing the " + samtools_flagstat_output_path
    if cmdStateCode == 0:
        print(stderrMsgSTATE_OK)
    else:
        print(stderrMsgSTATE_FAIL)


    ## 2. find total reads of tcga bam
    cmd = "samtools view -c " + bam_path + " >>" + total_reads_tcga_bam_output_path
    # print(cmd)

    cmdStateCode = os.system(cmd)
    stderrMsgSTATE_OK = cmd + ' => Done!'
    stderrMsgSTATE_FAIL = "Error in processing the " + total_reads_tcga_bam_output_path
    if cmdStateCode == 0:
        print(stderrMsgSTATE_OK)
    else:
        print(stderrMsgSTATE_FAIL)
        error_time_count += 1


    ## 3. find total reads of tcga bam human mapped
    cmd = "samtools view -c -F 260 " + bam_path + " >>" + total_reads_tcga_bam_human_mapped_output_path
    # print(cmd)

    cmdStateCode = os.system(cmd)
    stderrMsgSTATE_OK = cmd + ' => Done!'
    stderrMsgSTATE_FAIL = "Error in processing the " + samtools_flagstat_output_path
    if cmdStateCode == 0:
        print(stderrMsgSTATE_OK)
    else:
        print(stderrMsgSTATE_FAIL)
        error_time_count += 1


    ## 4. find total reads of human unmapped fastq
    # echo $(cat /TCGA-BRCA/GBM/unmapped_fastq/e45a44b3-a124-4a4c-9fec-3991740ada26.fastq | wc -l)/4|bc
    cmd = "echo $(cat " + human_unmaped_fastq_path + " | wc -l)/4|bc " + " >>" + total_reads_human_unmapped_fastq_output_path
    # print(cmd)

    cmdStateCode = os.system(cmd)
    stderrMsgSTATE_OK = cmd + ' => Done'
    stderrMsgSTATE_FAIL = "Error in processing the " + samtools_flagstat_output_path
    if cmdStateCode == 0:
        print(stderrMsgSTATE_OK)
    else:
        print(stderrMsgSTATE_FAIL)
        error_time_count += 1


    ## 5. find total reads of srnanalyzer Processed.fa
    # >1-3408
    # TCGTTTCCGGCTCGAAGGACCA
    # >2-3158
    # AAGCTGCCAGTTGAAGAACTGTATATCGT
    # >3-3105
    # TCTTTGGTTATCTAGCTGTATGAAACTCGT

    total_count = 0
    with open(srnanalyzer_processed_fa_path, 'r') as f:
        for line in f:
            if line[0] == ">":
                count = int(line.rstrip().split("-")[1])
                total_count += count
    
    print("total counts in sRNAnalyzer Processed.fa: ", total_count, "\n\n")

    with open(total_reads_srnanalyzer_processed_fa_output_path, "a") as f:
        f.write(str(total_count)+"\n")


print("Done!")
print("error_time_count: ", error_time_count)