#parse bacteria_associated_gene_ontology table

from os import listdir


cancer_type_file_path = "/home/kaipu/projects/bic/data/cancer_type.txt"

gene_ontology_biological_process_file_path = "/home/kaipu/projects/bic/data/gene_ontology_biological_process.txt"

taxonomy_level_bacterium_file_path = "/home/kaipu/projects/bic/data/taxonomy_level_bacterium.txt"

bacteria_associated_gene_ontology_output_path = "/home/kaipu/projects/bic/data/bacteria_associated_gene_ontology.txt"




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




gene_ontology_biological_process_dict = {}
with open(gene_ontology_biological_process_file_path) as f:
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] == "gene_ontology_biological_process_id":
            continue
        term_name = "GO_" + "_".join(line[1].split(" "))
        gene_ontology_biological_process_dict[term_name] = line[0]

print("len(gene_ontology_biological_process_dict): ", len(gene_ontology_biological_process_dict))
print(dict(list(gene_ontology_biological_process_dict.items())[0:3]))






bacteria_associated_gene_ontology_list = []
bacteria_associated_gene_ontology_list.append([ "bacteria_associated_gene_ontology_id", \
                                                "p_value", \
                                                "adjust_p_value", \
                                                "es", \
                                                "nes", \
                                                "n_more_extreme", \
                                                "size", \
                                                "gene_symbol", \
                                                "cancer_type_id", \
                                                "taxonomy_level_bacterium_id", \
                                                "gene_ontology_biological_process_id" ])




bacteria_associated_gene_ontology_id_count = 1


## subset test
# cancer_type_list = ["ACC", "BRCA"]
# taxonomy_level_list = ["Genus", "Phylum"]


for i in range(len(taxonomy_level_list)):
    taxonomy_level = taxonomy_level_list[i]

    aliquot_barcode_list = []
    taxonomy_level_bacterium_dict = {}
    with open(taxonomy_level_bacterium_file_path) as f:
        for line in f:
            line = line.rstrip().split("\t")
            if line[0] == "taxonomy_level_bacterium_id":
                continue
            if line[-1] == taxonomy_level[0]:
                taxonomy_level_bacterium_dict[line[2]] = line[0]

    print("len(taxonomy_level_bacterium_dict): ", len(taxonomy_level_bacterium_dict))


    for cancer_type in cancer_type_list:
        print("cancer_type: ", cancer_type)

        bacteria_gene_ontology_dir = "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/" + cancer_type + "/" + taxonomy_level + "/BacteriaGeneCorr/Percentage_20/Gene_Ontology/"

        bacteria_list = listdir(bacteria_gene_ontology_dir)
        print("len(bacteria_list): ", len(bacteria_list))

        for bacteria_name in bacteria_list:
            print("bacteria: ", bacteria_name)
            taxonomy_level_bacterium_id = taxonomy_level_bacterium_dict[bacteria_name]

            bacteria_gene_ontology_file_path = bacteria_gene_ontology_dir + "/" + bacteria_name + "/GO_BP.txt"

            with open(bacteria_gene_ontology_file_path) as f:
                next(f) # skip header line
                for line in f:
                    line = line.rstrip().split("\t")
                    # pathway 0
                    # pval 1
                    # padj 2
                    # ES 3
                    # NES 4
                    # nMoreExtreme 5 
                    # size 6 
                    # EntrezID 7
                    # GeneSymbol 8
                    term = line[0]
                    p_value =  line[1]
                    adjust_p_value = line[2]
                    es = line[3]
                    nes = line[4]
                    n_more_extreme = line[5]
                    size = line[6]
                    gene_symbol = line[8]
                    gene_ontology_biological_process_id = gene_ontology_biological_process_dict[term]

                    bacteria_associated_gene_ontology_list.append([ str(bacteria_associated_gene_ontology_id_count), \
                                                                    p_value, \
                                                                    adjust_p_value, \
                                                                    es, \
                                                                    nes, \
                                                                    n_more_extreme, \
                                                                    size, \
                                                                    gene_symbol, \
                                                                    cancer_type, \
                                                                    taxonomy_level_bacterium_id, \
                                                                    gene_ontology_biological_process_id])

                    bacteria_associated_gene_ontology_id_count += 1

    

print("bacteria_associated_gene_ontology_id_count: ", bacteria_associated_gene_ontology_id_count)
print("len(bacteria_associated_gene_ontology_list): ", len(bacteria_associated_gene_ontology_list))
print("head of bacteria_associated_gene_ontology_list: ", bacteria_associated_gene_ontology_list[0:4])
print("tail of bacteria_associated_gene_ontology_list: ", bacteria_associated_gene_ontology_list[-2:])



fp = open(bacteria_associated_gene_ontology_output_path, 'w')
for item in bacteria_associated_gene_ontology_list:
    fp.write("\t".join(item)+"\n")
fp.close()


