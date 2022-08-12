## parse taxonomy_level_bacterium table

taxonomy_level_list = ["Genus", "Family", "Order", "Class", "Phylum"]
taxonomy_level_list_header = ["p__", "c__", "o__", "f__", "g__"]
taxonomy_level_to_id_dict = dict()

#for i in range(len(taxonomy_level_list)):
#    taxonomy_level_to_id_dict[taxonomy_level_list[i]] = i+1

taxonomy_level_to_id_dict["Genus"] = "G"
taxonomy_level_to_id_dict["Family"] = "F"
taxonomy_level_to_id_dict["Order"] = "O"
taxonomy_level_to_id_dict["Class"] = "C"
taxonomy_level_to_id_dict["Phylum"] = "P"


print(taxonomy_level_to_id_dict)

taxonomy_level_list_header_length = len(taxonomy_level_list_header)
print(taxonomy_level_list_header_length)


taxonomy_level_bacterium_id = 1
taxonomy_level_bacterium_list = []
taxonomy_level_bacterium_list.append(["taxonomy_level_bacterium_id", "taxonomy_string_name", "name", "taxonomy_level_id"])

for i in range(len(taxonomy_level_list)):
    taxonomy_level = taxonomy_level_list[i]
    print(taxonomy_level)

    taxaFilePath = "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/Taxonomy/" + \
        "relative_abundance_matched_taxa_" + taxonomy_level + ".txt"
    print(taxaFilePath)
    
    with open(taxaFilePath, 'r') as f:
        for line in f:
            line = line.rstrip().split("\t")
            name = line[0]
            
            taxonomy_name_list = line[1:]
            taxonomy_level_list_header_select = taxonomy_level_list_header[:len(taxonomy_name_list)]
            taxonomy_string_list = []
            for j in range(len(taxonomy_name_list)):
                taxonomy_string_list.append(taxonomy_level_list_header_select[j] + taxonomy_name_list[j])
            
            taxonomy_string_name =  ";".join(taxonomy_string_list)
            taxonomy_level_bacterium_list.append([str(taxonomy_level_bacterium_id), taxonomy_string_name, name, str(taxonomy_level_to_id_dict[taxonomy_level])])
            taxonomy_level_bacterium_id += 1


print("len(taxonomy_level_bacterium_list): ", len(taxonomy_level_bacterium_list))

fp = open("/home/kaipu/projects/bic/data/taxonomy_level_bacterium.txt", 'w')
for item in taxonomy_level_bacterium_list:
    fp.write("\t".join(item)+"\n")
fp.close()