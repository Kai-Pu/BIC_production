## parse following tables
'''
    patient
    patient_sample_type
'''

aliquot_file_path = "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/TCGA_miRNAseq_and_mRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal_v2.txt"
#cancer_type_file_path = "/home/kaipu/projects/bic/data/cancer_type.txt"
sample_type_file_path = "/home/kaipu/projects/bic/data/sample_type.txt"

patient_output_path = "/home/kaipu/projects/bic/data/patient_v2.txt"
#patient_sample_type_output_path = "/home/kaipu/projects/bic/data/patient_sample_type.txt"
patient_cancer_sample_type_output_path = "/home/kaipu/projects/bic/data/patient_cancer_sample_type.txt"


#cancer_type_dict = {}
#with open(cancer_type_file_path, 'r') as f:
#    for line in f:
#        line = line.rstrip().split("\t")
#        if line[0] == "cancer_type_id":
#            continue
#        cancer_type_dict[line[1]] = line[0]

#print(cancer_type_dict)


sample_type_dict = {}
with open(sample_type_file_path, 'r') as f:
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] == "sample_type_id":
            continue
        sample_type_dict[line[1]] = line[0]

print(sample_type_dict)


patient_barcode_to_id_dict = {}
patient_list = []
patient_list.append(["patient_id", \
                     "patient_barcode", \
                     "gender", \
                     "age_at_initial_pathologic_diagnosis", \
                     "race_simplified", \
                     "stage", \
                     "stage_early_median_late", \
                     "stage_early_late", \
                     "os", \
                     "os_time", \
                     "cancer_type_id"])


#patient_sample_type_list = []
#patient_sample_type_list.append(["patient_sample_type_id", \
#                                 "sample_barcode", \
#                                 "aliquot_barcode_mirna_seq", \
#                                 "aliquot_barcode_rna_seq", \
#                                 "sub_cluster", \
#                                 "patient_id", \
#                                 "sample_type_id"])

patient_cancer_sample_type_list = []
patient_cancer_sample_type_list.append(["patient_cancer_sample_type_id", \
                                 "sample_barcode", \
                                 "aliquot_barcode_mirna_seq", \
                                 "aliquot_barcode_rna_seq", \
                                 "sub_cluster", \
                                 "patient_id", \
                                 "cancer_type_id", \
                                 "sample_type_id"])

patient_id_count = 1
patient_cancer_sample_type_id_count = 1


with open(aliquot_file_path, 'r') as f:
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] == "Patient_barcode":
            continue
        patient_barcode = line[0]
        gender = line[15]
        age_at_initial_pathologic_diagnosis = line[18]
        race_simplified = line[14].rstrip()
        stage = line[10]
        stage_early_median_late = line[12]
        stage_early_late = line[13]
        os = line[16]
        os_time = line[17]

        cancer_type = line[1]
        #cancer_type_id = cancer_type_dict[cancer_type]

        sample_barcode = line[2]
        aliquot_barcode_mirna_seq = line[3]
        aliquot_barcode_rna_seq = line[4]
        sub_cluster = line[7]
        
        sample_type = line[5]
        sample_type_id = sample_type_dict[sample_type]

        if patient_barcode not in patient_barcode_to_id_dict:
            patient_barcode_to_id_dict[patient_barcode] = patient_id_count
            patient_list.append([str(patient_id_count), \
                                 patient_barcode, \
                                 gender, \
                                 age_at_initial_pathologic_diagnosis, \
                                 race_simplified, \
                                 stage, \
                                 stage_early_median_late, \
                                 stage_early_late, \
                                 os, \
                                 os_time, \
                                 cancer_type])
#                                 cancer_type_id])
            patient_id_count += 1
        
#        patient_sample_type_list.append([str(patient_sample_type_id_count), \
#                                         sample_barcode, \
#                                         aliquot_barcode_mirna_seq, \
#                                         aliquot_barcode_rna_seq, \
#                                         sub_cluster, \
#                                         str(patient_barcode_to_id_dict[patient_barcode]), \
#                                         sample_type_id])
#        
#        patient_sample_type_id_count += 1

        patient_cancer_sample_type_list.append([str(patient_cancer_sample_type_id_count), \
                                         sample_barcode, \
                                         aliquot_barcode_mirna_seq, \
                                         aliquot_barcode_rna_seq, \
                                         sub_cluster, \
                                         str(patient_barcode_to_id_dict[patient_barcode]), \
                                         cancer_type, \
                                         sample_type_id])
        
        patient_cancer_sample_type_id_count += 1


print("patient_id_count: ", patient_id_count)
print("length of patient_list: ", len(patient_list))
#print("patient_sample_type_id_count: ", patient_sample_type_id_count)
#print("length of patient_sample_type_list: ", len(patient_sample_type_list))
print("patient_cancer_sample_type_id_count: ", patient_cancer_sample_type_id_count)
print("length of patient_cancer_sample_type_list: ", len(patient_cancer_sample_type_list))


## check id consistency
#check_patient_sample_type_id = 0
#for i in range(1, len(patient_sample_type_list)):
#    patient_barcode = '-'.join(patient_sample_type_list[i][1].split("-")[:3])
#    if int(patient_sample_type_list[i][5]) == patient_barcode_to_id_dict[patient_barcode]:
#        check_patient_sample_type_id += 1
#
# print("check_patient_sample_type_id: ", check_patient_sample_type_id)

check_cancer_patient_sample_type_id = 0
for i in range(1, len(patient_cancer_sample_type_list)):
    patient_barcode = '-'.join(patient_cancer_sample_type_list[i][1].split("-")[:3])
    if int(patient_cancer_sample_type_list[i][5]) == patient_barcode_to_id_dict[patient_barcode]:
        check_cancer_patient_sample_type_id += 1

print("check_cancer_patient_sample_type_id: ", check_cancer_patient_sample_type_id)

fp = open(patient_output_path, 'w')
for item in patient_list:
    fp.write("\t".join(item)+"\n")
fp.close()


#fp = open(patient_sample_type_output_path, 'w')
#for item in patient_sample_type_list:
#    fp.write("\t".join(item)+"\n")
#fp.close()

fp = open(patient_cancer_sample_type_output_path, 'w')
for item in patient_cancer_sample_type_list:
    fp.write("\t".join(item)+"\n")
fp.close()