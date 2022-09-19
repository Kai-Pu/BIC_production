#parse reactome_pathway table



gsea_msigdb_reactome_file_path = "/home/kaipu/microbiome_TCGA_miRNA/Reactome/c2.cp.reactome.v7.5.1.symbols.gmt"

reactome_pathway_output_path = "/home/kaipu/projects/bic_data/reactome_pathway.txt"




reactome_pathway_list = []

reactome_pathway_list.append([ "reactome_pathway_id", "term", "size", "gene_symbol" ])

reactome_pathway_id = 1
with open(gsea_msigdb_reactome_file_path, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        term = " ".join(line[0][9:].split("_"))
        # gene_symbol = line[2:]
        gene_symbol = sorted(line[2:])
        size = len(gene_symbol)

        reactome_pathway_list.append([ str(reactome_pathway_id), term, str(size), ",".join(gene_symbol) ])
        reactome_pathway_id += 1

print("len(reactome_pathway_list): ", len(reactome_pathway_list))
print("reactome_pathway_id: ", reactome_pathway_id)

print("head or tail of reactome_pathway_list: ")
print(reactome_pathway_list[0])
print(reactome_pathway_list[1])
print(reactome_pathway_list[-1])






fp = open(reactome_pathway_output_path, 'w')
for item in reactome_pathway_list:
    fp.write("\t".join(item)+"\n")
fp.close()