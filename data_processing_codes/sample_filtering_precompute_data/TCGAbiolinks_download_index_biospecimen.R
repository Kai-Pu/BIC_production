## TCGAbiolinks: find FFPE from biospecimen 20200817

library("TCGAbiolinks")
setwd("/home/kaipu/microbiome_TCGA_miRNA/biospecimen")

## load TCGA project table
TCGA_project_df <- read.table("../TCGA_project_name_case_sample_mirna_num.txt", 
                              sep="\t", header=T)
dim(TCGA_project_df)
head(TCGA_project_df, n=3)


## download indexed biospecimen
?GDCquery_clinic

#i <- 1
for (i in 1:nrow(TCGA_project_df)) {
  projectName <- paste("TCGA-", TCGA_project_df[i, 1], sep="")
  biospecimen <- GDCquery_clinic(project=projectName, type="biospecimen", save.csv=F)
  outputFile <- paste("TCGA-", TCGA_project_df[i, 1], 
                      "_biospecimen_indexed_TCGAbiolinks.txt", sep="")
  biospecimen_output <- cbind(biospecimen$submitter_id, 
                              biospecimen$is_ffpe,
                              biospecimen$sample_type,
                              biospecimen$updated_datetime,
                              biospecimen$created_datetime)
  colnames(biospecimen_output) <- c("submitter_id", 
                                    "is_ffpe",
                                    "sample_type",
                                    "updated_datetime",
                                    "created_datetime")
  write.table(biospecimen_output, file=outputFile, sep="\t",
              col.names=T, row.names=F, quote=F)
}

#head(biospecimen_output)

## test TCGA-BLCA
#biospecimen <- GDCquery_clinic(project="TCGA-BLCA", type="biospecimen")
#dim(biospecimen)
#summary(biospecimen)
#head(biospecimen, n=2)
#class(biospecimen)



## download xml and parse (skipped this part)
#query <- GDCquery(project="TCGA-BLCA", data.category="Biospecimen", file.type="xml")
#GDCdownload(query)
#biospecimen_xml_sample <- GDCprepare_clinic(query, clinical.info = "sample")
#biospecimen_xml_bio_patient <- GDCprepare_clinic(query, clinical.info = "bio_patient")
#biospecimen_xml_analyte <- GDCprepare_clinic(query, clinical.info = "analyte")
#biospecimen_xml_aliquot <- GDCprepare_clinic(query, clinical.info = "aliquot")
#dim(biospecimen_xml_aliquot)
#biospecimen_xml_protocol <- GDCprepare_clinic(query, clinical.info = "protocol")
