#!/usr/bin/env python

import pandas as pd
import os, sys

working_dir = sys.argv[1]
selected_sample_path = sys.argv[2]
merge_count_table_path = sys.argv[3]
taxa_level = sys.argv[4]


#working_dir = "/Users/chenkaipu/CancerMicorbiome/microbiome_TCGA_miRNA/Count_unique/Merge_count/Merge_count_nonFFPE/test/"
#selected_sample_path = "/Users/chenkaipu/CancerMicorbiome/microbiome_TCGA_miRNA/Count_unique/Merge_count/Merge_count_nonFFPE/test/selected_sample_table.txt"
#merge_count_table_path = "/Users/chenkaipu/CancerMicorbiome/microbiome_TCGA_miRNA/Count_unique/Merge_count/Merge_count_nonFFPE/test/TCGA_3_projects_merge_counts_nonFFPE_Genus.txt"
#taxa_level = "Genus"


## output header line: aliquot barcode
output = merge_count_table_path.split("/")[-1].replace(taxa_level+".txt", "selected_PT_STN_"+taxa_level+".txt")

print(sys.version)
print( "Working directory: " + working_dir )
print( "Selected sample path: " + selected_sample_path )
print( "Merged count table path: " + merge_count_table_path )
print( "Taxa level: " + taxa_level )
print( "Output count filename: " + output )

os.chdir(working_dir)


## columns in selected_sample_path
## Patient_barcode Project_name    Sample_barcode  Aliquot_barcode miRNAseq_file_uuid      Sample_type     Filtered_Aliquot_barcode
## TCGA-AO-A12G    BRCA    TCGA-AO-A12G-01A        TCGA-AO-A12G-01A-11R-A10I-13    8d28b296-94b2-4a07-bd7b-506e3a9af6c7    Primary Tumor
selected_sample_list = [taxa_level]

with open(selected_sample_path, 'r') as f:
	for line in f:
		line = line.rstrip().split("\t")
		if line[0] == "Patient_barcode":
			continue
		selected_sample_list.append(line[3])

print( "Selected sample number:", len(selected_sample_list) )



## merge_count_table import as pd object
## two header line!
## line 1: project name
## line 2: aliquot barcode
df = pd.read_table(merge_count_table_path, sep="\t", skiprows=1, header=0, index_col=False)
print(df)
#print(df.columns)
print( "Merged count table dimension: {}".format(df.shape) )



## extract columns of selected aliquot barcode and save to output
df_selected = df.loc[:, selected_sample_list]
print(df_selected)
print( "Selected count table dimension: {}".format(df_selected.shape) )
print( df_selected.iloc[0:3, 0:3] )

df_selected.to_csv(output, header=True, index=False, sep="\t", mode="w")


