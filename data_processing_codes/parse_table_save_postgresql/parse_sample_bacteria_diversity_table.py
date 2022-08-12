#parse sample_bacteria_diversity table


patient_cancer_sample_type_file_path = "/home/kaipu/projects/bic/data/patient_cancer_sample_type.txt"

sample_bacteria_diversity_output_path = "/home/kaipu/projects/bic/data/sample_bacteria_diversity.txt"



taxonomy_level_list = ["Genus", "Family", "Order", "Class", "Phylum"]


patient_cancer_sample_type_dict = {}
with open(patient_cancer_sample_type_file_path) as f:
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] == "patient_cancer_sample_type_id":
            continue
        patient_cancer_sample_type_dict[line[2]] = line[0]

print("len(patient_cancer_sample_type_dict): ", len(patient_cancer_sample_type_dict))
#print(dict(list(patient_cancer_sample_type_dict.items())[0:3]))



sample_bacteria_diversity_list = []
sample_bacteria_diversity_list.append(["sample_bacteria_diversity_id", \
                                       "diversity_method", \
                                       "diversity_value", \
                                       "taxonomy_level_id", \
                                       "patient_cancer_sample_type_id"])

sample_bacteria_diversity_id_count = 1
for i in range(len(taxonomy_level_list)):
    taxonomy_level = taxonomy_level_list[i]
    taxonomy_level_id = taxonomy_level[0]
    print("taxonomy_level_id: ", taxonomy_level_id)

    sample_bacteria_diversity_path = "/home/kaipu/microbiome_TCGA_miRNA/Diversity/TCGA_32_projects_merge_diversity_" + taxonomy_level + ".txt"

    with open(sample_bacteria_diversity_path) as f:
        for line in f:
            line = line.rstrip().split("\t")
            if line[0] == "Patient_barcode":
                continue
            
            aliquot_barcode_mirna_seq = line[3]
            patient_cancer_sample_type_id = patient_cancer_sample_type_dict[aliquot_barcode_mirna_seq]
            evenness = line[7]
            richness = line[8]
            shannon_index = line[9]
            simpson_index_of_diversity = line[10]
            simpson_reciprocal_index = line[11]

            sample_bacteria_diversity_list.append([str(sample_bacteria_diversity_id_count), \
                                                   "Evenness", evenness, \
                                                   taxonomy_level_id, patient_cancer_sample_type_id])
            sample_bacteria_diversity_id_count += 1

            sample_bacteria_diversity_list.append([str(sample_bacteria_diversity_id_count), \
                                                   "Richness", richness, \
                                                   taxonomy_level_id, patient_cancer_sample_type_id])
            sample_bacteria_diversity_id_count += 1
            
            sample_bacteria_diversity_list.append([str(sample_bacteria_diversity_id_count), \
                                                   "Shannon index", shannon_index, \
                                                   taxonomy_level_id, patient_cancer_sample_type_id])
            sample_bacteria_diversity_id_count += 1

            sample_bacteria_diversity_list.append([str(sample_bacteria_diversity_id_count), \
                                                   "Simpson index of diversity", simpson_index_of_diversity, \
                                                   taxonomy_level_id, patient_cancer_sample_type_id])
            sample_bacteria_diversity_id_count += 1

            sample_bacteria_diversity_list.append([str(sample_bacteria_diversity_id_count), \
                                                   "Simpson reciprocal index", simpson_reciprocal_index, \
                                                   taxonomy_level_id, patient_cancer_sample_type_id])
            sample_bacteria_diversity_id_count += 1

            
            

print("sample_bacteria_diversity_id_count: ", sample_bacteria_diversity_id_count)
print("len(sample_bacteria_diversity_list): ", len(sample_bacteria_diversity_list))
print("head of sample_bacteria_diversity_list: ", sample_bacteria_diversity_list[0:4])
print("tail of sample_bacteria_diversity_list: ", sample_bacteria_diversity_list[-2:])



fp = open(sample_bacteria_diversity_output_path, 'w')
for item in sample_bacteria_diversity_list:
    fp.write("\t".join(item)+"\n")
fp.close()