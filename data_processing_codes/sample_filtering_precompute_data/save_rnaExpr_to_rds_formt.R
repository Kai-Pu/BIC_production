## fix output under following directory
workingDirRoot <- "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/"
setwd(workingDirRoot)
getwd()


## human RNA expression downloaded from PanCanAtlas
rnaExprPath <- "/home/kaipu/microbiome_TCGA_miRNA/PanCanAtlas/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"
file.exists(rnaExprPath)

rnaExpr <- read.table(rnaExprPath, header=T, sep="\t", check.names=F, row.names=1) 
dim(rnaExpr)
rnaExpr[1:3, 1:3]

saveRDS(rnaExpr, file="rnaExpr.rds")## save it into rnaExpr.rds for faster input next time 


rnaExpr2 <- readRDS("rnaExpr.rds")
identical(rnaExpr, rnaExpr2)
dim(rnaExpr2)
rnaExpr2[1:3, 1:3]


