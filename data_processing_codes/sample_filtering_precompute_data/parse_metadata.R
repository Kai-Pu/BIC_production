## fix output under following directory
workingDirRoot <- "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/"
setwd(workingDirRoot)
getwd()


## alliquot ids of miRNAseq and RNAseq with best k in Cluster column
aliquotInfoPath <- "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/TCGA_miRNAseq_and_mRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal.txt"
file.exists(aliquotInfoPath)

aliquotInfo <- read.table(aliquotInfoPath, header=T, sep="\t", check.names=F)
dim(aliquotInfo)
aliquotInfo[1:2, ]



## clinical data downloaded from PanCanAtlas
clinicalInfoPath <- "/home/kaipu/microbiome_TCGA_miRNA/PanCanAtlas/TCGA-CDR.txt"
file.exists(clinicalInfoPath)

clinicalInfo <- read.table(clinicalInfoPath, header=T, sep="\t", check.names=F, row.names=1)
dim(clinicalInfo)
clinicalInfo[1:2, ]


## check overlapping between patient barcode of clinical (downloaded from PanCanAtlas) and miRNAseq
## There are 40 patient barcode of aliquotnfo not found in clinicalInfo
miRNAseqPatBar <- unique(as.character(aliquotInfo$Patient_barcode))
length(miRNAseqPatBar)
clinicalPatBar <- unique(as.character(clinicalInfo$bcr_patient_barcode))
length(clinicalPatBar)
length(intersect(miRNAseqPatBar, clinicalPatBar))

length(miRNAseqPatBar) - length(intersect(miRNAseqPatBar, clinicalPatBar))

miRNAseqPatBarOnly <- miRNAseqPatBar[! miRNAseqPatBar %in% clinicalPatBar]
aliquotInfoMiRNAseqOnly <- aliquotInfo[c(aliquotInfo$Patient_barcode %in% miRNAseqPatBarOnly), ]


aggregate(Patient_barcode ~ Project_name + Sample_type,
          data=aliquotInfoMiRNAseqOnly, FUN=length)

table(clinicalInfo$clinical_stage[clinicalInfo$type=="OV"])
table(clinicalInfo$ajcc_pathologic_tumor_stage[clinicalInfo$type=="OV"])
clinicalInfo[1, ]

## match patient_barcode of clinicalInfo with aliquotInfo's
clinicalInfoMatchIndex <- match(aliquotInfo$Patient_barcode, clinicalInfo$bcr_patient_barcode)
clinicalInfoMatch <- clinicalInfo[clinicalInfoMatchIndex, ]
dim(clinicalInfoMatch)
dim(aliquotInfo)
clinicalInfoMatch[1:3, ]
sum(as.vector(is.na(clinicalInfoMatch[which(is.na(clinicalInfoMatchIndex)), ])))/ncol(clinicalInfoMatch)
table(as.character(aliquotInfo$Patient_barcode)==as.character(clinicalInfoMatch$bcr_patient_barcode))
clinicalInfoMatch$bcr_patient_barcode <- aliquotInfo$Patient_barcode

table(clinicalInfo[, c("type", "ajcc_pathologic_tumor_stage")])
table(clinicalInfo[, c("type", "clinical_stage")])



ajcc_pathologic_tumor_stage <- as.character(clinicalInfoMatch$ajcc_pathologic_tumor_stage)
ajcc_pathologic_tumor_stage[which(aliquotInfo$Sample_type == "Solid Tissue Normal")] <- NA
table(ajcc_pathologic_tumor_stage)
sum(table(ajcc_pathologic_tumor_stage))
clinicalInfoMatch[which(ajcc_pathologic_tumor_stage=="I/II NOS"), ]

ajcc_pathologic_tumor_stage[which(ajcc_pathologic_tumor_stage=="[Discrepancy]")] <- NA
ajcc_pathologic_tumor_stage[which(ajcc_pathologic_tumor_stage=="[Not Applicable]")] <- NA
ajcc_pathologic_tumor_stage[which(ajcc_pathologic_tumor_stage=="[Not Available]")] <- NA
ajcc_pathologic_tumor_stage[which(ajcc_pathologic_tumor_stage=="[Unknown]")] <- NA
ajcc_pathologic_tumor_stage[which(ajcc_pathologic_tumor_stage=="I/II NOS")] <- NA
ajcc_pathologic_tumor_stage[which(ajcc_pathologic_tumor_stage=="Stage X")] <- NA

ajcc_pathologic_tumor_stage <- gsub(ajcc_pathologic_tumor_stage, pattern="Stage ", replacement="")
ajcc_pathologic_tumor_stage <- gsub(ajcc_pathologic_tumor_stage, pattern="A", replacement="")
ajcc_pathologic_tumor_stage <- gsub(ajcc_pathologic_tumor_stage, pattern="B", replacement="")
ajcc_pathologic_tumor_stage <- gsub(ajcc_pathologic_tumor_stage, pattern="C", replacement="")
ajcc_pathologic_tumor_stage <- gsub(ajcc_pathologic_tumor_stage, pattern="S", replacement="")

table(ajcc_pathologic_tumor_stage)


clinical_stage <- as.character(clinicalInfoMatch$clinical_stage)
clinical_stage[which(aliquotInfo$Sample_type == "Solid Tissue Normal")] <- NA
table(clinical_stage)
sum(table(clinical_stage))

clinical_stage[which(clinical_stage=="[Discrepancy]")] <- NA
clinical_stage[which(clinical_stage=="[Not Applicable]")] <- NA
clinical_stage[which(clinical_stage=="[Not Available]")] <- NA

clinical_stage <- gsub(clinical_stage, pattern="Stage ", replacement="")
clinical_stage <- toupper(clinical_stage)
clinical_stage <- gsub(clinical_stage, pattern="A", replacement="")
clinical_stage <- gsub(clinical_stage, pattern="B", replacement="")
clinical_stage <- gsub(clinical_stage, pattern="C", replacement="")
clinical_stage <- gsub(clinical_stage, pattern="S", replacement="")
clinical_stage <- gsub(clinical_stage, pattern="1", replacement="")
clinical_stage <- gsub(clinical_stage, pattern="2", replacement="")

table(clinical_stage)



aliquotInfo[1, ]
aliquotInfo$ajcc_pathologic_tumor_stage <- ajcc_pathologic_tumor_stage
table(aliquotInfo$ajcc_pathologic_tumor_stage)

aliquotInfo$clinical_stage <- clinical_stage
table(aliquotInfo$clinical_stage)

aliquotInfo$Stage <- NA
aliquotInfo$Stage_source <- NA


i <- 1
for (i in 1:nrow(colorCode31)) {
  projectName <- as.character(colorCode31[i, 1])
  projectRowIndex <- which(aliquotInfo$Project_name %in% projectName)
  if (length(table(aliquotInfo$ajcc_pathologic_tumor_stage[projectRowIndex]))>1) {
    aliquotInfo$Stage[projectRowIndex] <- aliquotInfo$ajcc_pathologic_tumor_stage[projectRowIndex]
    aliquotInfo$Stage_source[projectRowIndex] <- rep("ajcc_pathologic_tumor_stage", length(projectRowIndex))
  } else if (length(table(aliquotInfo$clinical_stage[projectRowIndex]))>1) {
    aliquotInfo$Stage[projectRowIndex] <- aliquotInfo$clinical_stage[projectRowIndex]
    aliquotInfo$Stage_source[projectRowIndex] <- rep("clinical_stage", length(projectRowIndex))
  } 
}


table(aliquotInfo$Project_name, aliquotInfo$Stage)
table(aliquotInfo$Project_name, aliquotInfo$Stage_source)


aliquotInfo$Stage_EarMedLat <- NA
aliquotInfo$Stage_EarMedLat[which(aliquotInfo$Stage=="I")] <- "I"
aliquotInfo$Stage_EarMedLat[which(aliquotInfo$Stage=="II" | aliquotInfo$Stage=="III")] <- "II&III"
aliquotInfo$Stage_EarMedLat[which(aliquotInfo$Stage=="IV")] <- "IV"
table(aliquotInfo$Project_name, aliquotInfo$Stage_EarMedLat)

aliquotInfo$Stage_EarLat <- NA
aliquotInfo$Stage_EarLat[which(aliquotInfo$Stage=="I" | aliquotInfo$Stage=="II")] <- "I&II"
aliquotInfo$Stage_EarLat[which(aliquotInfo$Stage=="III" | aliquotInfo$Stage=="IV")] <- "III&IV"
table(aliquotInfo$Project_name, aliquotInfo$Stage_EarLat)


## add simplified race column to aliquotInfo
aliquotInfo$Race_simplify <- "Other"
aliquotInfo$Race_simplify[which(clinicalInfoMatch$race=="ASIAN")] <- "Asian"
aliquotInfo$Race_simplify[which(clinicalInfoMatch$race=="BLACK OR AFRICAN AMERICAN")] <- "African American "
aliquotInfo$Race_simplify[which(clinicalInfoMatch$race=="WHITE")] <- "White"
table(aliquotInfo$Race_simplify)

head(aliquotInfo)



## add gender
aliquotInfo$Gender <- clinicalInfoMatch$gender




## add column from clinicalInfoMatch to aliquotinfo
aliquotInfo <- cbind(aliquotInfo, clinicalInfoMatch[, c("OS", "OS.time")])
head(aliquotInfo)
clinicalInfoMatch[1, ]
table(clinicalInfoMatch$gender)


## add age_at_initial_pathologic_diagnosis
aliquotInfo$age_at_initial_pathologic_diagnosis <- clinicalInfoMatch$age_at_initial_pathologic_diagnosis




## save additional columns (Race_simplify, ) for aliquotInfo
write.table(aliquotInfo, file="TCGA_miRNAseq_and_mRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal_v2.txt",
           sep="\t", col.names=T, row.names=F, quote=F)