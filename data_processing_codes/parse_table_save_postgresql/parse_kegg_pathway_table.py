#parse kegg_pathway table



gsea_msigdb_kegg_file_path = "/home/kaipu/microbiome_TCGA_miRNA/KEGG/c2.cp.kegg.v7.5.1.symbols.gmt"

kegg_pathway_output_path = "/home/kaipu/projects/bic_data/kegg_pathway.txt"




kegg_pathway_list = []

kegg_pathway_list.append([ "kegg_pathway_id", "term", "size", "gene_symbol" ])

kegg_pathway_id = 1
with open(gsea_msigdb_kegg_file_path, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        term = " ".join(line[0][5:].split("_"))
        # gene_symbol = line[2:]
        gene_symbol = sorted(line[2:])
        size = len(gene_symbol)

        kegg_pathway_list.append([ str(kegg_pathway_id), term, str(size), ",".join(gene_symbol) ])
        kegg_pathway_id += 1

print("len(kegg_pathway_list): ", len(kegg_pathway_list))
print("kegg_pathway_id: ", kegg_pathway_id)

print("head or tail of kegg_pathway_list: ")
print(kegg_pathway_list[0])
print(kegg_pathway_list[1])
print(kegg_pathway_list[-1])






fp = open(kegg_pathway_output_path, 'w')
for item in kegg_pathway_list:
    fp.write("\t".join(item)+"\n")
fp.close()