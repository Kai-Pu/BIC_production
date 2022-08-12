## remain one seq data for each sample type for each patient (mRNAseq)

import sys
working_dir = sys.argv[1]
mrnaseq_aliquot_bcr_path = sys.argv[2]
mirnaseq_nonFFPE_pt_stn_path = sys.argv[3]
output = "/home/kaipu/microbiome_TCGA_miRNA/FFPE_nonFFPE/nonFFPE_PT_STN/TCGA_miRNAseq_and_mRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal.txt"

## input RNAseq aliquot barcode and store as dict (key: sample_bcr, val: aliquot_bar)
mrnaseq_sample_bcr_dict = dict()
with open(mrnaseq_aliquot_bcr_path, 'r') as f:
	for line in f:
		mrnaseq_aliquot_bcr = line.rstrip()
		mrnaseq_sample_bcr = "-".join(mrnaseq_aliquot_bcr.split("-")[0:4])
		if mrnaseq_sample_bcr not in mrnaseq_sample_bcr_dict:
			mrnaseq_sample_bcr_dict[mrnaseq_sample_bcr] = [mrnaseq_aliquot_bcr]
		else:
			mrnaseq_sample_bcr_dict[mrnaseq_sample_bcr].append(mrnaseq_aliquot_bcr)

print( "Num. of sample bcr: ", len(mrnaseq_sample_bcr_dict) )
print( mrnaseq_sample_bcr )
print( mrnaseq_sample_bcr_dict[mrnaseq_sample_bcr] )


## input miRNAseq nonFFPE PT and STN remaining one table as list
mirnaseq_nonFFPE_pt_stn_list = []
with open(mirnaseq_nonFFPE_pt_stn_path, 'r') as f:
	for line in f:
		line = line.rstrip().split("\t")
		if line[0] == "Patient_barcode":
			continue
		mirnaseq_nonFFPE_pt_stn_list.append(line)

print( "Num. of remaining miRNAseq (PT and STN): ", len(mirnaseq_nonFFPE_pt_stn_list) )
print( '\t'.join(mirnaseq_nonFFPE_pt_stn_list[0]) )


## match mRNAseq sample_bcr with miRNAseq sample_bcr and remain the one with higher sorted value if different plate num existed
result_list = []
result_list.append( ["Patient_barcode", "Project_name", "Sample_barcode", "Aliquot_barcode_miRNAseq", "Aliquot_barcode_RNAseq", "Sample_type", "Filtered_Aliquot_barcode_RNAseq"] )
unmatched_mirnaseq_sample_bcr_list = []
for i in range(len(mirnaseq_nonFFPE_pt_stn_list)):
	mirnaseq_sample_bcr = mirnaseq_nonFFPE_pt_stn_list[i][2]
	matched_mrnaseq_aliquot_bcr_list = []
	sorted_matched_mrnaseq_aliquot_bcr_list = []
	matched_mrnaseq_aliquot_bcr = ''
	filtered_matched_mrnaseq_aliquot_bcr = ''
	try:
		matched_mrnaseq_aliquot_bcr_list = mrnaseq_sample_bcr_dict[mirnaseq_sample_bcr]
		#print(matched_mrnaseq_aliquot_bcr_list)
	except:
		unmatched_mirnaseq_sample_bcr_list.append(mirnaseq_sample_bcr)
		#print(mirnaseq_sample_bcr)

	if len(matched_mrnaseq_aliquot_bcr_list) > 1:
		sorted_matched_mrnaseq_aliquot_bcr_list = sorted(matched_mrnaseq_aliquot_bcr_list)
		matched_mrnaseq_aliquot_bcr = sorted_matched_mrnaseq_aliquot_bcr_list[0]
		filtered_matched_mrnaseq_aliquot_bcr = ",".join(sorted_matched_mrnaseq_aliquot_bcr_list[1:])
		print("sorted_matched_mrnaseq_aliquot_bcr_list: ", sorted_matched_mrnaseq_aliquot_bcr_list)
		print("matched_mrnaseq_aliquot_bcr: ", matched_mrnaseq_aliquot_bcr)
		print("filtered_matched_mrnaseq_aliquot_bcr: ", filtered_matched_mrnaseq_aliquot_bcr)
	elif len(matched_mrnaseq_aliquot_bcr_list) == 1:
		matched_mrnaseq_aliquot_bcr = matched_mrnaseq_aliquot_bcr_list[0]
		#print("matched_mrnaseq_aliquot_bcr: ", matched_mrnaseq_aliquot_bcr)

	result = mirnaseq_nonFFPE_pt_stn_list[i][0:4] + [matched_mrnaseq_aliquot_bcr] + [mirnaseq_nonFFPE_pt_stn_list[i][5]] + [filtered_matched_mrnaseq_aliquot_bcr]

	result_list.append(result)



print( mirnaseq_sample_bcr )
print( "Num. of unmatched miRNAseq sample barcode: ", len(unmatched_mirnaseq_sample_bcr_list) )
#print( "unmatched miRNAseq sample barcode: " + ", ".join(unmatched_mirnaseq_sample_bcr_list ) )








## output the mRNA and miRNA match table
fp = open(output, 'w')
for line in result_list:
	fp.write( "\t".join(line) + "\n" )
fp.close()




