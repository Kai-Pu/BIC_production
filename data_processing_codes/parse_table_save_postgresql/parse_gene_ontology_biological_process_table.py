#parse gene_ontology_biological_process table



gsea_msigdb_bp_file_path = "/home/kaipu/microbiome_TCGA_miRNA/Gene_Ontology/c5.go.bp.v7.2.symbols.gmt"

gene_ontology_biological_process_output_path = "/home/kaipu/projects/bic/data/gene_ontology_biological_process.txt"




gene_ontology_biological_process_list = []

gene_ontology_biological_process_list.append([ "gene_ontology_biological_process_id", "term", "size", "gene_symbol" ])

gene_ontology_biological_process_id = 1
with open(gsea_msigdb_bp_file_path, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        term = " ".join(line[0][3:].split("_"))
        # gene_symbol = line[2:]
        gene_symbol = sorted(line[2:])
        size = len(gene_symbol)

        gene_ontology_biological_process_list.append([ str(gene_ontology_biological_process_id), term, str(size), ",".join(gene_symbol) ])
        gene_ontology_biological_process_id += 1

print("len(gene_ontology_biological_process_list): ", len(gene_ontology_biological_process_list))
print("gene_ontology_biological_process_id: ", gene_ontology_biological_process_id)

print("head or tail of gene_ontology_biological_process_list: ")
print(gene_ontology_biological_process_list[0])
print(gene_ontology_biological_process_list[1])
print(gene_ontology_biological_process_list[-1])






fp = open(gene_ontology_biological_process_output_path, 'w')
for item in gene_ontology_biological_process_list:
    fp.write("\t".join(item)+"\n")
fp.close()