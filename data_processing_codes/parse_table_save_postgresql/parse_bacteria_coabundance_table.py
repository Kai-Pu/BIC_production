#parse bacteria_coabundance table


taxonomy_level_bacterium_file_path = "/home/kaipu/projects/bic/data/taxonomy_level_bacterium.txt"

cancer_type_file_path = "/home/kaipu/projects/bic/data/cancer_type.txt"

bacteria_coabundance_output_path = "/home/kaipu/projects/bic/data/bacteria_coabundance.txt"
bacteria_coabundance_num_summary_output_path = "/home/kaipu/projects/bic/data/bacteria_coabundance_num_summary.txt"


taxonomy_level_list = ["Genus", "Family", "Order", "Class", "Phylum"]
taxonomy_level_id_to_level_dict = dict()

taxonomy_level_id_to_level_dict["G"] = "Genus"
taxonomy_level_id_to_level_dict["F"] = "Family"
taxonomy_level_id_to_level_dict["O"] = "Order"
taxonomy_level_id_to_level_dict["C"] = "Class"
taxonomy_level_id_to_level_dict["P"] = "Phylum"

print(taxonomy_level_id_to_level_dict)



cancer_type_list = []
with open(cancer_type_file_path) as f:
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] == "cancer_type_id" or line[0] == "GBM":
            continue
        cancer_type_list.append(line[0])

print("len(cancer_type_list): ", len(cancer_type_list))
print("cancer_type_list: ", cancer_type_list)


bacteria_coabundance_num_summary = []
bacteria_coabundance_num_summary.append(["Taxonomy_level"] + cancer_type_list)





## test 
# cancer_type_list = ["ACC", "BLCA"]
# taxonomy_level_list = ["Phylum", "Order"]

bacteria_coabundance_list = []

bacteria_coabundance_list.append(["bacteria_coabundance_id", \
                                  "cancer_type_id", \
                                  "taxonomy_level_bacterium_id_1_id", \
                                  "taxonomy_level_bacterium_id_2_id", \
                                  "sparscc", \
                                  "p_value"])


bacteria_coabundance_id_count = 1


for taxonomy_level in taxonomy_level_list:
    print("taxonomy_level: ", taxonomy_level)

    taxonomy_level_bacterium_dict = {}
    with open(taxonomy_level_bacterium_file_path) as f:
        for line in f:
            line = line.rstrip().split("\t")
            if line[0] == "taxonomy_level_bacterium_id":
                continue
            if line[-1] == taxonomy_level[0]:
                taxonomy_level_bacterium_dict[line[2]] = line[0]

    print("len(taxonomy_level_bacterium_dict): ", len(taxonomy_level_bacterium_dict))

    tmp_num_summary_list = [taxonomy_level]

    for cancer_type in cancer_type_list:
        print("cancer_type: ", cancer_type)

        bacteria_sparcc_cor_path = "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/" + cancer_type  + "/" + taxonomy_level + "/bacteria_Coabundance/bootstrapTimes10000/" + cancer_type + "_bacteria_coabundance_sparcc_correlation_Percentage_20_" + taxonomy_level + ".txt"

        i = 0
        with open(bacteria_sparcc_cor_path) as f:
            bacteria_list = next(f).rstrip().split("\t")
            #print("bacteria_list: ", bacteria_list)
            print("len(bacteria_list): ", len(bacteria_list))
            tmp_num_summary_list.append(str(len(bacteria_list)))

            for line in f:
                line = line.rstrip().split("\t")
                #print(line)
                bacteria_1 = line[0]
                bacteria_1_id = taxonomy_level_bacterium_dict[bacteria_1]
                i += 1
                for j in range(i, len(line)):
                    bacteria_2 = bacteria_list[j-1]
                    bacteria_2_id = taxonomy_level_bacterium_dict[bacteria_2]
                    sparcc_cor = line[j]
                    if bacteria_1 != bacteria_2:
                        bacteria_coabundance_list.append([ str(bacteria_coabundance_id_count), cancer_type, bacteria_1_id, bacteria_2_id, sparcc_cor ])
                        #print(bacteria_coabundance_list[bacteria_coabundance_id_count])
                        bacteria_coabundance_id_count += 1
    
    bacteria_coabundance_num_summary.append(tmp_num_summary_list)

                    

# print("bacteria_coabundance_id_count: ", bacteria_coabundance_id_count)
# print("len(bacteria_coabundance_list): ", len(bacteria_coabundance_list))
# print("head of bacteria_coabundance_list: ", bacteria_coabundance_list[0:4])
# print("tail of bacteria_coabundance_list: ", bacteria_coabundance_list[-2:])




bacteria_coabundance_id_count = 1

for taxonomy_level in taxonomy_level_list:
    #print("taxonomy_level: ", taxonomy_level)

    taxonomy_level_bacterium_dict = {}
    with open(taxonomy_level_bacterium_file_path) as f:
        for line in f:
            line = line.rstrip().split("\t")
            if line[0] == "taxonomy_level_bacterium_id":
                continue
            if line[-1] == taxonomy_level[0]:
                taxonomy_level_bacterium_dict[line[2]] = line[0]

    #print("len(taxonomy_level_bacterium_dict): ", len(taxonomy_level_bacterium_dict))


    for cancer_type in cancer_type_list:
        #print("cancer_type: ", cancer_type)

        bacteria_sparcc_p_val_path = "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/" + cancer_type  + "/" + taxonomy_level + "/bacteria_Coabundance/bootstrapTimes10000/" + cancer_type + "_bacteria_coabundance_sparcc_p_Val_Percentage_20_" + taxonomy_level + ".txt"


        with open(bacteria_sparcc_p_val_path) as f:
            bacteria_list = next(f).rstrip().split("\t")
            #print("bacteria_list: ", bacteria_list)
            #print("len(bacteria_list): ", len(bacteria_list))
            i = 0

            for line in f:
                line = line.rstrip().split("\t")
                #print(line)
                bacteria_1 = line[0]
                bacteria_1_id = taxonomy_level_bacterium_dict[bacteria_1]
                i += 1
                for j in range(i, len(line)):
                    bacteria_2 = bacteria_list[j-1]
                    bacteria_2_id = taxonomy_level_bacterium_dict[bacteria_2]
                    sparcc_p_val = line[j]
                    if bacteria_1 != bacteria_2:
                        bacteria_coabundance_list[bacteria_coabundance_id_count].append(sparcc_p_val)
                        ##bacteria_coabundance_list[bacteria_coabundance_id_count].append([ str(bacteria_coabundance_id_count), bacteria_1_id, bacteria_2_id, sparcc_p_val ])
                        bacteria_coabundance_id_count += 1



print("bacteria_coabundance_id_count: ", bacteria_coabundance_id_count)
print("len(bacteria_coabundance_list): ", len(bacteria_coabundance_list))
print("head of bacteria_coabundance_list: ", bacteria_coabundance_list[0:4])
print("tail of bacteria_coabundance_list: ", bacteria_coabundance_list[-2:])




# (["bacteria_coabundance_id", \ 0
# "cancer_type_id", \ 1
# "taxonomy_level_bacterium_id_1", \ 2
# "taxonomy_level_bacterium_id_2", \ 3
# "sparscc", \ 4
# "p_value"]) 5

fp = open(bacteria_coabundance_output_path, 'w')
for i in range(len(bacteria_coabundance_list)):
    column_reorder_list = [bacteria_coabundance_list[i][0], \
                           bacteria_coabundance_list[i][4], \
                           bacteria_coabundance_list[i][5], \
                           bacteria_coabundance_list[i][1], \
                           bacteria_coabundance_list[i][2], \
                           bacteria_coabundance_list[i][3]]
    fp.write("\t".join(column_reorder_list)+"\n")
fp.close()


fp = open(bacteria_coabundance_num_summary_output_path, 'w')
for item in bacteria_coabundance_num_summary:
    fp.write("\t".join(item)+"\n")
fp.close()
