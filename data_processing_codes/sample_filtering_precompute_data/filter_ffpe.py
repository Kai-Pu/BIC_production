#!/usr/bin/env python

## overlap the sample barcode created by TCGAbiolinks and manifest file
## save the patient barocde, sample barcode, aliquot barcode, id, sample type


import os, sys

output_dir = sys.argv[1]
project_manifest_path = sys.argv[2]

print("Python version information\n", sys.version)
print("Working/Output directory: ", output_dir)
print("TCGA project table path: ", project_manifest_path)

os.chdir(output_dir)

project_name_list = []
manifest_path_list = []
biospecimen_path_list = []
output_list_ffpe = []
output_list_non_ffpe = []
output_result_summary = "TCGA_ffpe_non_ffpe_summary.txt"

with open(project_manifest_path, 'r') as f:
	for line in f:
		line = line.rstrip().split("\t")
		project_name = line[0]
		manifest_path = line[1]
		biospecimen_path = line[2]
		output_ffpe = "TCGA-" + project_name + "_patient_sample_aliquot_barcode_uuid_sample_type_ffpe.txt"
		output_non_ffpe = "TCGA-" + project_name + "_patient_sample_aliquot_barcode_uuid_sample_type_non_ffpe.txt"
		project_name_list.append(project_name)
		manifest_path_list.append(manifest_path)
		biospecimen_path_list.append(biospecimen_path)
		output_list_ffpe.append(output_ffpe)
		output_list_non_ffpe.append(output_non_ffpe)



result_summary = []
result_summary.append(["TCGA_project", "Case_number", "FFPE_number", "Non-FFPE_number", "Sample_barcode_not_in_biospecimen", "Duplicate_sample_barcode_in_biospecimen"])

for i in range(len(project_name_list)):
	biospecimen_dict = {}
	ffpe_list = []
	non_ffpe_list = []
	duplicate_sample_barcode_in_biospecimen = []
	not_found_sample_barcode_list = []
	count_biospecimen_line = 0
	count_manifest_line = 0

	
	print("===== TCGA-" + project_name_list[i] + " =====")

	
	with open(biospecimen_path_list[i], 'r') as f:
		for line in f:
			line = line.rstrip().split("\t")
			if line[0] == "submitter_id":
				continue
			if line[0] not in biospecimen_dict:
				biospecimen_dict[line[0]] = line[1:3]
			else:
				duplicate_sample_barcode_in_biospecimen.append(line[0])
			count_biospecimen_line += 1

	print("Count biospecimen line:", count_biospecimen_line)
	print("Number of sample barcode in biospecimen_dict:", len(biospecimen_dict))

	
	with open(manifest_path_list[i], 'r') as f:
		for line in f:
			line = line.rstrip().split("\t")
			if line[0] == "id":
				continue
			file_uuid = line[0]
			aliquot_barcode = line[1].split("_")[0]
			patient_barcode = "-".join(aliquot_barcode.split("-")[:3])
			sample_barcode = "-".join(aliquot_barcode.split("-")[:4])
			count_manifest_line += 1
			try:
				is_ffpe = biospecimen_dict[sample_barcode][0]
				sample_type = biospecimen_dict[sample_barcode][1]
			except:
				is_ffpe = "NA"
				sample_type = "NA"
				not_found_sample_barcode_list.append(sample_barcode)

			if is_ffpe == "FALSE":
				non_ffpe_list.append([patient_barcode, sample_barcode, aliquot_barcode, file_uuid, is_ffpe, sample_type])

			elif is_ffpe == "TRUE":
				ffpe_list.append([patient_barcode, sample_barcode, aliquot_barcode, file_uuid, is_ffpe, sample_type])

	print("Count manifest line:", count_manifest_line)
	print("Numbers of FFPE samples:", len(ffpe_list))
	print("Numbers of non-FFPE samples:", len(non_ffpe_list))
	print("Not found in biospecimen:", ",".join(not_found_sample_barcode_list))
	result_summary.append([project_name_list[i], str(count_manifest_line), str(len(ffpe_list)), str(len(non_ffpe_list)), ",".join(not_found_sample_barcode_list), ",".join(duplicate_sample_barcode_in_biospecimen)])

	fp = open(output_list_ffpe[i], "w")
	fp.write("\t".join(["Patient_barcode", "Sample_barcode", "Aliquot_barcode", "miRNA-seq_file_uuid", "is_ffpe", "Sample_type"]) + "\n")
	for line in ffpe_list:
		fp.write("\t".join(line) + "\n")
	fp.close()


	fp = open(output_list_non_ffpe[i], "w")
	fp.write("\t".join(["Patient_barcode", "Sample_barcode", "Aliquot_barcode", "miRNA-seq_file_uuid", "is_ffpe", "Sample_type"]) + "\n")
	for line in non_ffpe_list:
		fp.write("\t".join(line) + "\n")
	fp.close()


fp = open(output_result_summary, "w")
for line in result_summary:
	fp.write("\t".join(line) + "\n")
fp.close()





