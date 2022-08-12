#!/usr/bin/env python

## remain one seq count for each sample type for each patient

## I. sample type:
## Primary Tumor (PT) and Solid Tissue Normal (STN)
## Metastatic (M), Recurrent Tumor (RT), LAML (188 Primary Blood Derived Cancer - Peripheral Blood) and
## other (1 Additional Metastatic + 11 Additional - New Primary)

## II. Replicate keeping priority:
## - For those with different analyte codes: remain R over T/H
## - For same analyte code: remain the higher value one after ordered

## save the patient barocde, sample barcode, aliquot barcode, id, sample type


import os, sys

non_ffpe_table_dir = sys.argv[1]
project_name_path = sys.argv[2]

print("Python version information\n", sys.version)
print("NonFFPE table directory: " + non_ffpe_table_dir)
print("TCGA project name path: " + project_name_path)

os.chdir(non_ffpe_table_dir)


## create following folders if not existed
create_folder_list = ["nonFFPE_PT", "nonFFPE_STN", "nonFFPE_M", "nonFFPE_RT", "nonFFPE_other", "nonFFPE_all", "nonFFPE_PT_STN"]

for folder in create_folder_list:
	if not os.path.exists(folder):
		os.makedirs(folder)


## import TCGA project and saved as list
project_name_list = []
with open(project_name_path, 'r') as f:
	for line in f:
		line = line.rstrip()
		project_name_list.append(line)


result_selected_summary = []
result_selected_summary.append(["Project_name", "Primary_tumor", "Solid_tissue_normal", "Metastatic", "Recurrent_tumor", "Other"])

all_selected_list = []
all_filtered_list = []

PT_selected_list = []
STN_selected_list = []
M_selected_list = []
RT_selected_list = []
other_selected_list = []
PT_STN_selected_list = []



for project_name in project_name_list:
	## columns in non_ffpe_table: 
	## "Patient_barcode", "Sample_barcode", "Aliquot_barcode", "miRNA-seq_file_uuid", "is_ffpe", "Sample_type"
	non_ffpe_filename = "TCGA-" + project_name + "_patient_sample_aliquot_barcode_uuid_sample_type_non_ffpe.txt"

	## use patient_barcode as key and aliquot_barcode as values to find replicate samples
	PT_dict = {}
	STN_dict = {}
	M_dict = {}
	RT_dict = {}
	other_dict = {}

	## use aliquot_barcode as key and store other informations as values
	all_aliquot_dict = {}

	with open(non_ffpe_filename, 'r') as f:
		for line in f:
			line = line.rstrip().split("\t")

			if line[0] == "Patient_barcode":
				continue

			patient_barcode = line[0]
			sample_barcode = line[1]
			aliquot_barcode = line[2]
			file_uuid = line[3]
			sample_type = line[5]

			information = [patient_barcode, project_name, sample_barcode, aliquot_barcode, file_uuid, sample_type]

			all_aliquot_dict[aliquot_barcode] = information


			if sample_type == "Primary Tumor":
				if patient_barcode not in PT_dict:
					PT_dict[patient_barcode] = [aliquot_barcode]
				else:
					PT_dict[patient_barcode].append(aliquot_barcode)

			elif sample_type == "Solid Tissue Normal":
				if patient_barcode not in STN_dict:
					STN_dict[patient_barcode] = [aliquot_barcode]
				else:
					STN_dict[patient_barcode].append(aliquot_barcode)

			elif sample_type == "Metastatic":
				if patient_barcode not in M_dict:
					M_dict[patient_barcode] = [aliquot_barcode]
				else:
					M_dict[patient_barcode].append(aliquot_barcode)			

			elif sample_type == "Recurrent Tumor":
				if patient_barcode not in RT_dict:
					RT_dict[patient_barcode] = [aliquot_barcode]
				else:
					RT_dict[patient_barcode].append(aliquot_barcode)

			else:
				if patient_barcode not in other_dict:
					other_dict[patient_barcode] = [aliquot_barcode]
				else:
					other_dict[patient_barcode].append(aliquot_barcode)				

		


	for pat_bcr, ali_bcr_list in PT_dict.items():
		if len(ali_bcr_list) >1:
			ali_bcr_selected = "NA"
			ali_bcr_list_sort = sorted(ali_bcr_list, reverse=True)
			analyte_code_list = []
			# TCGA-06-0678-11A-32R-A36C-13
			for ali_bcr in ali_bcr_list_sort:
				analyte_code = ali_bcr.split("-")[4][-1]
				analyte_code_list.append(analyte_code)

			## priority R>H>T
			if "R" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("R")]
			elif "H" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("H")]
			elif "T" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("T")]

			ali_bcr_list_sort.remove(ali_bcr_selected)
			ali_bcr_filtered_list = ali_bcr_list_sort
			for ali_bcr_filtered in ali_bcr_filtered_list:
				all_filtered_list.append(all_aliquot_dict[ali_bcr_filtered])
		else:
			ali_bcr_selected = ali_bcr_list[0]
			ali_bcr_filtered_list = [""]

		PT_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )
		PT_STN_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )
		all_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )


	for pat_bcr, ali_bcr_list in STN_dict.items():
		if len(ali_bcr_list) >1:
			ali_bcr_selected = "NA"
			ali_bcr_list_sort = sorted(ali_bcr_list, reverse=True)
			analyte_code_list = []
			for ali_bcr in ali_bcr_list_sort:
				analyte_code = ali_bcr.split("-")[4][-1]
				analyte_code_list.append(analyte_code)

			## priority R>H>T
			if "R" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("R")]
			elif "H" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("H")]
			elif "T" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("T")]

			ali_bcr_list_sort.remove(ali_bcr_selected)
			ali_bcr_filtered_list = ali_bcr_list_sort
			for ali_bcr_filtered in ali_bcr_filtered_list:
				all_filtered_list.append(all_aliquot_dict[ali_bcr_filtered])
		else:
			ali_bcr_selected = ali_bcr_list[0]
			ali_bcr_filtered_list = [""]

		STN_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )
		PT_STN_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )
		all_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )


	for pat_bcr, ali_bcr_list in M_dict.items():
		if len(ali_bcr_list) >1:
			ali_bcr_selected = "NA"
			ali_bcr_list_sort = sorted(ali_bcr_list, reverse=True)
			analyte_code_list = []
			for ali_bcr in ali_bcr_list_sort:
				analyte_code = ali_bcr.split("-")[4][-1]
				analyte_code_list.append(analyte_code)

			## priority R>H>T
			if "R" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("R")]
			elif "H" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("H")]
			elif "T" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("T")]

			ali_bcr_list_sort.remove(ali_bcr_selected)
			ali_bcr_filtered_list = ali_bcr_list_sort
			for ali_bcr_filtered in ali_bcr_filtered_list:
				all_filtered_list.append(all_aliquot_dict[ali_bcr_filtered])
		else:
			ali_bcr_selected = ali_bcr_list[0]
			ali_bcr_filtered_list = [""]

		M_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )
		all_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )


	for pat_bcr, ali_bcr_list in RT_dict.items():
		if len(ali_bcr_list) >1:
			ali_bcr_selected = "NA"
			ali_bcr_list_sort = sorted(ali_bcr_list, reverse=True)
			analyte_code_list = []
			for ali_bcr in ali_bcr_list_sort:
				analyte_code = ali_bcr.split("-")[4][-1]
				analyte_code_list.append(analyte_code)

			## priority R>H>T
			if "R" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("R")]
			elif "H" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("H")]
			elif "T" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("T")]

			ali_bcr_list_sort.remove(ali_bcr_selected)
			ali_bcr_filtered_list = ali_bcr_list_sort
			for ali_bcr_filtered in ali_bcr_filtered_list:
				all_filtered_list.append(all_aliquot_dict[ali_bcr_filtered])
		else:
			ali_bcr_selected = ali_bcr_list[0]
			ali_bcr_filtered_list = [""]

		RT_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )
		all_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )


	for pat_bcr, ali_bcr_list in other_dict.items():
		if len(ali_bcr_list) >1:
			ali_bcr_selected = "NA"
			ali_bcr_list_sort = sorted(ali_bcr_list, reverse=True)
			analyte_code_list = []
			for ali_bcr in ali_bcr_list_sort:
				analyte_code = ali_bcr.split("-")[4][-1]
				analyte_code_list.append(analyte_code)

			## priority R>H>T
			if "R" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("R")]
			elif "H" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("H")]
			elif "T" in analyte_code_list:
				ali_bcr_selected = ali_bcr_list_sort[analyte_code_list.index("T")]

			ali_bcr_list_sort.remove(ali_bcr_selected)
			ali_bcr_filtered_list = ali_bcr_list_sort
			for ali_bcr_filtered in ali_bcr_filtered_list:
				all_filtered_list.append(all_aliquot_dict[ali_bcr_filtered])
		else:
			ali_bcr_selected = ali_bcr_list[0]
			ali_bcr_filtered_list = [""]

		other_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )
		all_selected_list.append( all_aliquot_dict[ali_bcr_selected] + [",".join(ali_bcr_filtered_list)] )


	result_selected_summary.append( [project_name, len(PT_dict), len(STN_dict), len(M_dict), len(RT_dict), len(other_dict)] )
	print("=== TCGA-" + project_name + " ===")
	print("PT number:", len(PT_dict))
	print("STN number:", len(STN_dict))
	print("M number:", len(M_dict))
	print("RT number:", len(RT_dict))
	print("other number:", len(other_dict))


header_line = "\t".join( ["Patient_barcode", "Project_name", "Sample_barcode", "Aliquot_barcode", "miRNAseq_file_uuid", "Sample_type", "Filtered_Aliquot_barcode"] )
header_line_filetered = "\t".join( ["Patient_barcode", "Project_name", "Sample_barcode", "Aliquot_barcode", "miRNAseq_file_uuid", "Sample_type"] )

fp = open("./nonFFPE_PT/TCGA_miRNAseq_selected_samples_nonFFPE_Primary_Tumor.txt", 'w')
fp.write( header_line + "\n" )
for line in PT_selected_list:
	fp.write( "\t".join(line) + "\n" )
fp.close()


fp = open("./nonFFPE_STN/TCGA_miRNAseq_selected_samples_nonFFPE_Solid_Tissue_Normal.txt", 'w')
fp.write( header_line + "\n" )
for line in STN_selected_list:
	fp.write( "\t".join(line) + "\n" )
fp.close()


fp = open("./nonFFPE_PT_STN/TCGA_miRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal.txt", 'w')
fp.write( header_line + "\n" )
for line in PT_STN_selected_list:
	fp.write( "\t".join(line) + "\n" )
fp.close()


fp = open("./nonFFPE_M/TCGA_miRNAseq_selected_samples_nonFFPE_Metastatic.txt", 'w')
fp.write( header_line + "\n" )
for line in M_selected_list:
	fp.write( "\t".join(line) + "\n" )
fp.close()


fp = open("./nonFFPE_RT/TCGA_miRNAseq_selected_samples_nonFFPE_Recurrent_Tumor.txt", 'w')
fp.write( header_line + "\n" )
for line in RT_selected_list:
	fp.write( "\t".join(line) + "\n" )
fp.close()


fp = open("./nonFFPE_other/TCGA_miRNAseq_selected_samples_nonFFPE_Other.txt", 'w')
fp.write( header_line + "\n" )
for line in other_selected_list:
	fp.write( "\t".join(line) + "\n" )
fp.close()


fp = open("./nonFFPE_all/TCGA_miRNAseq_selected_samples_nonFFPE_all.txt", 'w')
fp.write( header_line + "\n" )
for line in all_selected_list:
	fp.write( "\t".join(line) + "\n" )
fp.close()


fp = open("./nonFFPE_all/TCGA_miRNAseq_filtered_samples_nonFFPE_all.txt", 'w')
fp.write( header_line_filetered + "\n" )
for line in all_filtered_list:
	fp.write( "\t".join(line) + "\n" )
fp.close()


fp = open("./nonFFPE_all/TCGA_miRNAseq_selected_samples_nonFFPE_all_summary.txt", 'w')
for line in result_selected_summary:
	fp.write( "\t".join(map(str, line)) + "\n" )
fp.close()
