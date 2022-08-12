#parse bacteria_human_rna_network table


taxonomy_level_bacterium_file_path = "/home/kaipu/projects/bic/data/taxonomy_level_bacterium.txt"

cancer_type_file_path = "/home/kaipu/projects/bic/data/cancer_type.txt"

bacteria_human_rna_network_output_path = "/home/kaipu/projects/bic/data/bacteria_human_rna_network.txt"
bacteria_human_rna_network_num_summary_output_path = "/home/kaipu/projects/bic/data/bacteria_human_rna_network_num_summary.txt"


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




bacteria_human_rna_network_num_summary = []
bacteria_human_rna_network_num_summary.append(["Taxonomy_level"] + cancer_type_list)



## test 
# cancer_type_list = ["ACC", "BLCA"]
# taxonomy_level_list = ["Phylum", "Order"]

bacteria_human_rna_network_list = []

bacteria_human_rna_network_list.append(["bacteria_human_rna_scc_id", \
                                        "cancer_type_id", \
                                        "taxonomy_level_bacterium_id", \
                                        "human_rna", \
                                        "scc", \
                                        "scc_p_value", \
                                        "fisher_z_transformed_scc"])


bacteria_human_rna_network_id_count = 1


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

        bacteria_scc_cor_path = "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/" + cancer_type  + "/" + taxonomy_level + "/BacteriaGeneCorr/Percentage_20/" + cancer_type + "_Bacteria_Gene_SCC_matrix.txt"

        with open(bacteria_scc_cor_path) as f:
            bacteria_list = next(f).rstrip().split("\t")
            #print("bacteria_list: ", bacteria_list)
            print("len(bacteria_list): ", len(bacteria_list))
            tmp_num_summary_list.append(str(len(bacteria_list)))

            for line in f:
                line = line.rstrip().split("\t")
                #print(line)
                human_rna = line[0]
                for j in range(1, len(line)):
                    bacteria = bacteria_list[j-1]
                    bacteria_id = taxonomy_level_bacterium_dict[bacteria]
                    scc = line[j]
                    bacteria_human_rna_network_list.append([ str(bacteria_human_rna_network_id_count), cancer_type, bacteria_id, human_rna, scc ])
                    bacteria_human_rna_network_id_count += 1
            
    bacteria_human_rna_network_num_summary.append(tmp_num_summary_list)


                    

print("bacteria_human_rna_network_id_count: ", bacteria_human_rna_network_id_count)
print("len(bacteria_human_rna_network_list): ", len(bacteria_human_rna_network_list))
print("head of bacteria_human_rna_network_list: ", bacteria_human_rna_network_list[0:4])
print("tail of bacteria_human_rna_network_list: ", bacteria_human_rna_network_list[-2:])




bacteria_human_rna_network_id_count = 1

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


    for cancer_type in cancer_type_list:

        bacteria_scc_p_value_path = "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/" + cancer_type  + "/" + taxonomy_level + "/BacteriaGeneCorr/Percentage_20/" + cancer_type + "_Bacteria_Gene_SCC_Pvalue_matrix.txt"


        with open(bacteria_scc_p_value_path) as f:
            bacteria_list = next(f).rstrip().split("\t")
            
            for line in f:
                line = line.rstrip().split("\t")
                human_rna = line[0]
                for j in range(1, len(line)):
                    bacteria = bacteria_list[j-1]
                    bacteria_id = taxonomy_level_bacterium_dict[bacteria]
                    scc_p_value = line[j]
                    bacteria_human_rna_network_list[bacteria_human_rna_network_id_count].append(scc_p_value)
                    bacteria_human_rna_network_id_count += 1


print("bacteria_human_rna_network_id_count: ", bacteria_human_rna_network_id_count)
print("len(bacteria_human_rna_network_list): ", len(bacteria_human_rna_network_list))
print("head of bacteria_human_rna_network_list: ", bacteria_human_rna_network_list[0:4])
print("tail of bacteria_human_rna_network_list: ", bacteria_human_rna_network_list[-2:])





bacteria_human_rna_network_id_count = 1

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


    for cancer_type in cancer_type_list:

        bacteria_fisher_z_ransform_scc_path = "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/" + cancer_type  + "/" + taxonomy_level + "/BacteriaGeneCorr/Percentage_20/" + cancer_type + "_Bacteria_Gene_SCC_Fisher_z_transform_matrix.txt"


        with open(bacteria_fisher_z_ransform_scc_path) as f:
            bacteria_list = next(f).rstrip().split("\t")
            
            for line in f:
                line = line.rstrip().split("\t")
                human_rna = line[0]
                for j in range(1, len(line)):
                    bacteria = bacteria_list[j-1]
                    bacteria_id = taxonomy_level_bacterium_dict[bacteria]
                    fisher_z_ransform_scc = line[j]
                    bacteria_human_rna_network_list[bacteria_human_rna_network_id_count].append(fisher_z_ransform_scc)
                    bacteria_human_rna_network_id_count += 1


print("bacteria_human_rna_network_id_count: ", bacteria_human_rna_network_id_count)
print("len(bacteria_human_rna_network_list): ", len(bacteria_human_rna_network_list))
print("head of bacteria_human_rna_network_list: ", bacteria_human_rna_network_list[0:4])
print("tail of bacteria_human_rna_network_list: ", bacteria_human_rna_network_list[-2:])




# ["bacteria_human_rna_scc_id",\ 0
# "human_rna",\ 3
# "scc",\ 4
# "fisher_z_transformed_scc",\  6
# "scc_p_value",\ 5
# "taxonomy_level_bacterium_id",\ 2
# "cancer_type_id"] 1

fp = open(bacteria_human_rna_network_output_path, 'w')
for i in range(len(bacteria_human_rna_network_list)):
    column_reorder_list = [bacteria_human_rna_network_list[i][0], \
                           bacteria_human_rna_network_list[i][3], \
                           bacteria_human_rna_network_list[i][4], \
                           bacteria_human_rna_network_list[i][6], \
                           bacteria_human_rna_network_list[i][5], \
                           bacteria_human_rna_network_list[i][2], \
                           bacteria_human_rna_network_list[i][1]]
    fp.write("\t".join(column_reorder_list)+"\n")
fp.close()





fp = open(bacteria_human_rna_network_num_summary_output_path, 'w')
for item in bacteria_human_rna_network_num_summary:
    fp.write("\t".join(item)+"\n")
fp.close()