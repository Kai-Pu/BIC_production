#!/usr/bin/env python

import pandas as pd
import os, sys

working_dir = sys.argv[1]
project_table_path = sys.argv[2]
taxa_level = sys.argv[3]
## columns of project_table_path (without header line):
## project_name, count_table_path, ffpe_table_path, non_ffpe_table_path

#working_dir = "/Users/chenkaipu/CancerMicorbiome/microbiome_TCGA_miRNA/Count_unique/Merge_count/Merge_count_nonFFPE/test"
#project_table_path = "/Users/chenkaipu/CancerMicorbiome/microbiome_TCGA_miRNA/Count_unique/Merge_count/Merge_count_nonFFPE/test/project_table.txt"
#taxa_level = "Genus"


print(sys.version)
print( "Working directory: " + working_dir )
print( "Project path: " + project_table_path )

os.chdir(working_dir)


## import project table and save in a list
## (index 0: project name; index 1: project count path; 
##  index 2: ffpe table path; index 3: non_ffpe table path)
project_list = []

with open(project_table_path, 'r') as f:
	for line in f:
		line = line.rstrip().split("\t")
		#print(line[0])
		#print(line[1])
		project_list.append(line)

print( "Project number:", len(project_list) )

output = "TCGA_" + str(len(project_list)) + "_projects_merge_counts_nonFFPE_" + taxa_level + ".txt"
output_summary = "TCGA_" + str(len(project_list)) + "_projects_merge_counts_nonFFPE_" + taxa_level + "_summary.txt"


## import ffpe table for removal 
## (key: project; value: list containing miRNA-seq_file_uuid)
remove_dict = {}
## No biospecimen information is found in the SKCM TCGA-XV-AB01-01A => remove it
## $ grep TCGA-XV-AB01-01A /TCGA-BRCA/SKCM/bam/gdc_manifest.2019-10-24.txt 
## 331a1cd8-c042-4c31-877f-817df309a1fb	TCGA-XV-AB01-01A-12R-A404-13_mirna_gdc_realn.bam	603e8991d96c6de4020b7dc28756b841	162091402	released
remove_dict["SKCM"] = ["331a1cd8-c042-4c31-877f-817df309a1fb"]
for i in range(len(project_list)):
	if project_list[i][0] not in remove_dict: 
		remove_dict[project_list[i][0]] = []
	with open(project_list[i][2], 'r') as f:
		for line in f:
			line = line.rstrip().split("\t")
			#print(project_list[i][0])
			#print(line)
			if line[0] == "Patient_barcode":
				continue
			## If empty row existed, line is [""]
			#if line != [""]:
				#remove_dict[project_list[i][0]].append(line[3])
			try: 
				remove_dict[project_list[i][0]].append(line[3])
			except:
				print("Empty row existed in:", project_list[i][2])

#print("FFPE dict: ", remove_dict)


#remove_dict["a"] = ["A111", "B0000"]
#remove_dict["b"] = ["9481", "B200"]
#remove_dict["c"] = []

## import non_ffpe for replacement of miRNA-seq_file_uuid to aliquot_barcode 
## (key: project; value: dict with key: miRNA-seq_file_uuid, and value: aliquot_barcode
uuid_2_aliq_dict = {}
for i in range(len(project_list)):
	uuid_2_aliq_dict[project_list[i][0]] = {}
	with open(project_list[i][3], 'r') as f:
		for line in f:
			line = line.rstrip().split("\t")
			if line[0] == "Patient_barcode":
				continue
			#print(line)
			uuid_2_aliq_dict[project_list[i][0]][line[3]] = line[2]

#print("uuid 2 aliquot_barcode dict: ", uuid_2_aliq_dict)


## save count tables into dataframe inside a dict
## summarize the sample numbers, taxa numbers (without unknown) of projects
result_summary = [ ["Project_name", "miRNA-seq_file_num", "taxa_num"] ]
df_dict = {}
for i in range(len(project_list)):
	df_dict[project_list[i][0]] = pd.read_table(project_list[i][1], sep="\t", header=0, index_col=0)
	
	print("Before removal of ffpe: " )
	print( df_dict[project_list[i][0]] )
	print( project_list[i][0] + " dimension: {}".format(df_dict[project_list[i][0]].shape) )

	print("After removal of ffpe: " )
	df_dict[project_list[i][0]].drop(columns=remove_dict[project_list[i][0]], inplace=True)
	## Note: inplace=False: it will return a copy and do not change the original one
	#df_dict[project_list[i][0]] = df_dict[project_list[i][0]].drop(columns=remove_dict[project_list[i][0]])
	print( df_dict[project_list[i][0]] )
	print( project_list[i][0] + " dimension: {}".format(df_dict[project_list[i][0]].shape) )

	print("After removal of ffpe and replace uuid to aliquot_barcode: " )
	df_dict[project_list[i][0]].rename(columns=uuid_2_aliq_dict[project_list[i][0]], inplace=True)
	print( df_dict[project_list[i][0]] )
	print( project_list[i][0] + " dimension: {}".format(df_dict[project_list[i][0]].shape) )

	summary = [ project_list[i][0], len(df_dict[project_list[i][0]].columns), len(df_dict[project_list[i][0]].index)-1 ]
	result_summary.append(summary)




## merge count table accross projects which doesn't contain FFPE samples and unknown taxa
df_merge = pd.concat(df_dict, axis=1, join="outer", sort=True)

print( "Merged table: ")
print( df_merge )
print( "Merged dimension: {}".format(df_merge.shape) )

print("After removal of unknown row: ")
df_merge.drop(index="unknown", inplace=True)
print( df_merge )
print( "Merged dimension: {}".format(df_merge.shape) )

## replace NaN with 0
print("Replace NaN with 0: ")
df_merge.fillna(0, inplace=True)
print( df_merge )
print( "Merged dimension: {}".format(df_merge.shape) )

## turn float points back to integers
df_merge = df_merge.astype(int)
print("Turn float points back to integers: ")
print( df_merge )
print( "Merged dimension: {}".format(df_merge.shape) )


#print("Column names: ")
#print(list(df_merge.columns))

result_summary.append( [ "Sum", len(df_merge.columns), len(df_merge.index) ] )
print(result_summary)


## save to txt file
headers = list(df_merge.columns)
header_project = ["Project_name"]
header_taxa = [taxa_level]
for i in range(len(headers)):
	header_project.append(headers[i][0])
	header_taxa.append(headers[i][1])
with open(output, 'w') as f:
	f.write( "\t".join(header_project) + "\n")
	f.write( "\t".join(header_taxa) + "\n")
df_merge.to_csv(output, header=False, index=True, sep="\t", mode='a')


with open(output_summary, 'w') as f:
	for line in result_summary:
		f.write("\t".join(map(str, line)) + "\n")

