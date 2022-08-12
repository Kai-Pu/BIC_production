## parse statistics table

cancer_type_file_path = "/home/kaipu/projects/bic/data/cancer_type.txt"
patient_file_path = "/home/kaipu/projects/bic/data/patient.txt"
patient_cancer_sample_type_file_path = "/home/kaipu/projects/bic/data/patient_cancer_sample_type.txt"

statistics_output = "/home/kaipu/projects/bic/data/statistics.txt"

cancer_type_list = []
with open(cancer_type_file_path) as f:
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] == "cancer_type_id":
            continue
        cancer_type_list.append(line[0:2])

print("len(cancer_type_list): ", len(cancer_type_list))
print("cancer_type_list: ", cancer_type_list)



patient_cancer_sample_type_list = []
with open(patient_cancer_sample_type_file_path) as f:
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] == "patient_cancer_sample_type_id":
            continue
        patient_cancer_sample_type_list.append(line[-2:])

print("len(patient_cancer_sample_type_list): ", len(patient_cancer_sample_type_list))
print("patient_cancer_sample_type_list[0:5]: ", patient_cancer_sample_type_list[0:5])



patient_list = []
with open(patient_file_path) as f:
    for line in f:
        line = line.rstrip().split("\t")
        if line[0] == "patient_id":
            continue
        patient_list.append(line)

print("len(patient_list): ", len(patient_list))
print("patient_list[0:5]: ", patient_list[0:5])




for i in range(len(cancer_type_list)):
    cancer_type = cancer_type_list[i][0]
    print("cancer_type: ", cancer_type)
    number_tumor = 0
    number_normal = 0
    survival_data_available = "-"
    race_data_available = "-"
    stage_data_available = "-"
    
    for j in range(len(patient_cancer_sample_type_list)):
        if patient_cancer_sample_type_list[j][0] == cancer_type:
            if patient_cancer_sample_type_list[j][1] == "TP":
                number_tumor += 1
            if patient_cancer_sample_type_list[j][1] == "TN":
                number_normal += 1


    for k in range(len(patient_list)):
        if patient_list[k][-1] == cancer_type:
            if patient_list[k][-2] != "NA":
                survival_data_available = "available"
            if patient_list[k][4] != "NA":
                if patient_list[k][4] != "Other":
                    race_data_available = "available"
            if patient_list[k][5] != "NA":
                stage_data_available = "available"
            
    cancer_type_list[i] = cancer_type_list[i] + [str(number_tumor), str(number_normal), survival_data_available, race_data_available, stage_data_available]




print("cancer_type_list: ", cancer_type_list)

fp = open(statistics_output, 'w')
fp.write("\t".join(["cancer_type", "number_of_tumor_sample", "number_of_normal_sample_adjacent_to_tumor", "survival_data_available", "race_data_available", "stage_data_available"])+"\n")
for item in cancer_type_list:
    tmp = item[1] + " (" + item[0] + ")"
    tmp = [tmp] + item[2:]
    fp.write("\t".join(tmp)+"\n")
fp.close()



