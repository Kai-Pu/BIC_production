#parse bacteria_expression table


taxonomy_level_bacterium_file_path = "/home/kaipu/projects/bic/data/taxonomy_level_bacterium.txt"
#patient_sample_type_file_path = "/home/kaipu/projects/bic/data/patient_sample_type.txt"
patient_cancer_sample_type_file_path = "/home/kaipu/projects/bic/data/patient_cancer_sample_type.txt"

bacteria_expression_output_path = "/home/kaipu/projects/bic/data/bacteria_expression.txt"



taxonomy_level_list = ["Genus", "Family", "Order", "Class", "Phylum"]
taxonomy_level_id_to_level_dict = dict()
#for i in range(len(taxonomy_level_list)):
#    taxonomy_level_id_to_level_dict[i+1] = taxonomy_level_list[i]

taxonomy_level_id_to_level_dict["G"] = "Genus"
taxonomy_level_id_to_level_dict["F"] = "Family"
taxonomy_level_id_to_level_dict["O"] = "Order"
taxonomy_level_id_to_level_dict["C"] = "Class"
taxonomy_level_id_to_level_dict["P"] = "Phylum"

print(taxonomy_level_id_to_level_dict)



#patient_sample_type_dict = {}
#with open(patient_sample_type_file_path) as f:
#    for line in f:
#        line = line.rstrip().split("\t")
#        if line[0] == "patient_sample_type_id":
#            continue
#        patient_sample_type_dict[line[2]] = line[0]
#
#print("len(patient_sample_type_dict): ", len(patient_sample_type_dict))
#print(dict(list(patient_sample_type_dict.items())[0:3]))


patient_cancer_sample_type_dict = {}
with open(patient_cancer_sample_type_file_path) as f:
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] == "patient_cancer_sample_type_id":
            continue
        patient_cancer_sample_type_dict[line[2]] = line[0]

print("len(patient_cancer_sample_type_dict): ", len(patient_cancer_sample_type_dict))
#print(dict(list(patient_cancer_sample_type_dict.items())[0:3]))



bacteria_expression_list = []
#bacteria_expression_list.append(["bacteria_expression_id", \
#                                 "relative_abundance_gmpr_normalized", \
#                                 "count_gmpr_normalized", \
#                                 "taxonomy_level_bacterium_id", \
#                                 "patient_sample_type_id"])

bacteria_expression_list.append(["bacteria_expression_id", \
                                 "relative_abundance_gmpr_normalized", \
                                 "count_gmpr_normalized", \
                                 "taxonomy_level_bacterium_id", \
                                 "patient_cancer_sample_type_id"])

bacteria_expression_id_count = 1
for i in range(len(taxonomy_level_list)):
    taxonomy_level = taxonomy_level_list[i]

    bacteria_relative_abundance_path = "/home/kaipu/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/TCGA_32_projects_merge_relative_abundance_nonFFPE_selected_PT_STN_GMPR_norm_" + taxonomy_level + ".txt"
    
    aliquot_barcode_list = []
    taxonomy_level_bacterium_dict = {}
    with open(taxonomy_level_bacterium_file_path) as f:
        for line in f:
            line = line.rstrip().split("\t")
            if line[0] == "taxonomy_level_bacterium_id":
                continue
#            if (int(line[-1]) == (i+1)):
            if line[-1] == taxonomy_level[0]:
                taxonomy_level_bacterium_dict[line[2]] = line[0]

    print("len(taxonomy_level_bacterium_dict): ", len(taxonomy_level_bacterium_dict))
    #print(dict(list(taxonomy_level_bacterium_dict.items())[0:3]))



    count_gmpr_normalized = "NA"
    with open(bacteria_relative_abundance_path) as f:
        for line in f:
            line = line.rstrip().split("\t")
            if line[0] == taxonomy_level:
                #header = line
                aliquot_barcode_list = line[1:]
                #break
                continue
            
            bacteria_name = line[0]
            taxonomy_level_bacterium_id = taxonomy_level_bacterium_dict[bacteria_name]
            
            baceria_expression = line[1:]     
            relative_abundance_gmpr_normalized = "NA"
            
            for j in range(len(aliquot_barcode_list)):
                bacteria_expression_id = bacteria_expression_id_count
                aliquot_barcode = aliquot_barcode_list[j]
#                patient_sample_type_id = patient_sample_type_dict[aliquot_barcode]
                patient_cancer_sample_type_id = patient_cancer_sample_type_dict[aliquot_barcode]
                relative_abundance_gmpr_normalized = baceria_expression[j]
                bacteria_expression_list.append([str(bacteria_expression_id), \
                                                 relative_abundance_gmpr_normalized, \
                                                 count_gmpr_normalized, \
                                                 taxonomy_level_bacterium_id, \
                                                 patient_cancer_sample_type_id])
#                                                 patient_sample_type_id])
                bacteria_expression_id_count += 1

    #print("relative_abundance", taxonomy_level, len(header))
    #print("head of aliquot_barcode_list: ", aliquot_barcode_list[0:3])

print("bacteria_expression_id_count: ", bacteria_expression_id_count)
print("len(bacteria_expression_list): ", len(bacteria_expression_list))
print("head of bacteria_expression_list: ", bacteria_expression_list[0:4])
print("tail of bacteria_expression_list: ", bacteria_expression_list[-2:])




bacteria_expression_id_count = 1
for i in range(len(taxonomy_level_list)):
    taxonomy_level = taxonomy_level_list[i]

    bacteria_count_gmpr_path = "/home/kaipu/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/TCGA_32_projects_merge_counts_nonFFPE_selected_PT_STN_GMPR_norm_" + taxonomy_level + ".txt"
    #bacteria_count_path = "/home/kaipu/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/TCGA_32_projects_merge_counts_nonFFPE_" + taxonomy_level + ".txt"
    

    #with open(bacteria_count_path) as f:
    #    for line in f:
    #        line = line.rstrip().split("\t")
    #        if line[0] == "Project_name":
    #            project_header = line
    #            break

    #print("count", taxonomy_level, len(project_header))


    with open(bacteria_count_gmpr_path) as f:
        for line in f:
            line = line.rstrip().split("\t")
            if line[0] == taxonomy_level:
                #header = line
                aliquot_barcode_list = line[1:]
                continue
            
            baceria_expression = line[1:]     
            
            for j in range(len(aliquot_barcode_list)):
                count_gmpr_normalized = baceria_expression[j]
                bacteria_expression_list[bacteria_expression_id_count][2] = count_gmpr_normalized
                bacteria_expression_id_count += 1

    #print("relative_abundance", taxonomy_level, len(header))
    #print("head of aliquot_barcode_list: ", aliquot_barcode_list[0:3])

print("bacteria_expression_id_count: ", bacteria_expression_id_count)
print("len(bacteria_expression_list): ", len(bacteria_expression_list))
print("head of bacteria_expression_list: ", bacteria_expression_list[0:4])
print("tail of bacteria_expression_list: ", bacteria_expression_list[-2:])



fp = open(bacteria_expression_output_path, 'w')
for item in bacteria_expression_list:
    fp.write("\t".join(item)+"\n")
fp.close()