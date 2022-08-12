### use Phylum level for skmeans clustering
### remember to modify the following variables

taxaLevel <-  "Phylum" 



## fix output under following directory
workingDirRoot <- "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/"
setwd(workingDirRoot)
getwd()



## load library
library("skmeans")
library('ggplot2')
library("ggpubr")
library("gridExtra")
library("pheatmap")
library("cluster")
require("dplyr")

## bacteria relative abundance
bacteriaExprPath <- paste0("/home/kaipu/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/TCGA_32_projects_merge_relative_abundance_nonFFPE_selected_PT_STN_GMPR_norm_", taxaLevel, ".txt")
file.exists(bacteriaExprPath)

bacteriaExpr <- read.table(bacteriaExprPath, header=T, sep="\t", check.names=F, row.names=1)
dim(bacteriaExpr)
bacteriaExpr[1:3, 1:3]
rownames(bacteriaExpr)

## alliquot ids of miRNAseq and RNAseq
aliquotInfoPath <- "/home/kaipu/microbiome_TCGA_miRNA/FFPE_nonFFPE/nonFFPE_PT_STN/TCGA_miRNAseq_and_mRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal.txt"
file.exists(aliquotInfoPath)

aliquotInfo <- read.table(aliquotInfoPath, header=T, sep="\t", check.names=F)
dim(aliquotInfo)
aliquotInfo[1:2, ]



## color codes for projects (without LAML and GBM)
colorCode31Path <- "/home/kaipu/microbiome_TCGA_miRNA/PanCanAtlas/color_codes_without_LAML_GBM.txt"
file.exists(colorCode31Path)

colorCode31 <- read.table(colorCode31Path, header=T, check.names=F, comment.char="")
dim(colorCode31)


## check if the order of bacteriaExpr and aliquotInfo are the same
table(colnames(bacteriaExpr) == aliquotInfo$Aliquot_barcode_miRNAseq)
table(aliquotInfo$Sample_type)
table(aliquotInfo$Project_name)


## do spherical kmeans and find best k by highest mean silhouette value
projectBestK <- data.frame(ProjectName=as.character(colorCode31[, 1]), k=NA)
projectBestKSampleNum <- list()
colorCode31[, 1]
aliquotInfo$Cluster <- NA
#i <- 1
for (i in 1:nrow(colorCode31)) {
  
  projectName <- as.character(colorCode31[i, 1])
  projectNameIndex <- which(aliquotInfo$Project_name %in% projectName)
  #length(projectNameIndex)
  bacteriaExprProject <- bacteriaExpr[, projectNameIndex]
  dim(bacteriaExprProject)
  #table(colnames(bacteriaExprProject) == aliquotInfo$Aliquot_barcode_miRNAseq[projectNameIndex])
  
  projectSkmeans2to10 <- matrix(NA, ncol=9, nrow=ncol(bacteriaExprProject))
  rownames(projectSkmeans2to10) <- colnames(bacteriaExprProject)
  colnames(projectSkmeans2to10) <- 2:10
  dim(projectSkmeans2to10)
  projectSkmeans2to10[1:3, ]
  
  projectSkmeans2to10SilhouetteMean <- as.data.frame(matrix(NA, ncol=2, nrow=9))
  projectSkmeans2to10SilhouetteMean[, 1] <- 2:10
  colnames(projectSkmeans2to10SilhouetteMean) <- c("k", "MeanSilhouette")
  
  dir.create(file.path(getwd(), projectName, taxaLevel, "Skmeans", "Silhouette"), recursive=T, showWarnings=F)
  
  #j <- 2
  for (j in 2:10) {
    projectSkmeans <- skmeans(t(bacteriaExprProject), k=j)
    projectSkmeans2to10[, j-1] <- projectSkmeans$cluster
    projectSilhouette <- silhouette(projectSkmeans$cluster, dist(t(bacteriaExprProject)))
    #table(projectSkmeans$cluster == projectSilhouette[, 1])
    #table(names(projectSkmeans$cluster)==colnames(bacteriaExprProject))
    projectSkmeans2to10SilhouetteMean[j-1, 2] <- mean(projectSilhouette[, 3])
    
    silhouetteK <- projectSilhouette[, 1:3]
    rownames(silhouetteK) <- colnames(bacteriaExprProject)
    silhouetteKoutput <- file.path(getwd(), projectName, taxaLevel, "Skmeans", "Silhouette", paste0(projectName, "_Skmeans_Silhouette_K_", j, ".txt"))
    write.table(silhouetteK, file=silhouetteKoutput, sep="\t", col.names=T, row.names=T, quote=F)
  }
  
  silhouetteMeanOutput <- file.path(getwd(), projectName, taxaLevel, "Skmeans", "Silhouette", paste0(projectName, "_Skmeans_Silhouette_K_Mean_summary.txt"))
  write.table(projectSkmeans2to10SilhouetteMean, file=silhouetteMeanOutput, sep="\t", col.names=T, row.names=F, quote=F)
  
  silhouetteMeanPlotOutput <- file.path(getwd(), projectName, taxaLevel, "Skmeans", "Silhouette", paste0(projectName,"_Skmeans_Silhouette_K_Mean_plot.svg"))
  svg(filename=silhouetteMeanPlotOutput, width=5, height=5)
  plot(c(1, projectSkmeans2to10SilhouetteMean[, 1]), c(0, projectSkmeans2to10SilhouetteMean[, 2]),
       type="o", pch=20, frame=T, 
       xlab="", ylab="",
       cex.axis=0.5, cex.lab=0.7,
       yaxt="n", xaxt="n", las=1, tck=-0.05, mgp=c(3,0.5,0), col="#3498db",
       main=paste0(projectName, "\n(spherical k-means clustering)"))
  axis(side=1, at=c(1, projectSkmeans2to10SilhouetteMean[, 1]), cex.axis=0.7, tck=-0.01, mgp=c(3,0.2,0))
  axis(side=2, at=seq(0, max(round(projectSkmeans2to10SilhouetteMean[, 2], 2))+0.10, by=0.05), 
       cex.axis=0.7, tck=-0.01, mgp=c(3,0.5,0), las=1)
  mtext(side=1, text="Number of clusters (k)", line=1.5, cex=1)
  mtext(side=2, text="Average Silhouette width", line=2.5, cex=1)
  dev.off()
  
  
  projectSkmeans2to10SilhouetteMean <- projectSkmeans2to10SilhouetteMean[order(projectSkmeans2to10SilhouetteMean$MeanSilhouette, decreasing=T), ]
  bestK <- projectSkmeans2to10SilhouetteMean[1, 1]
  
  projectBestK[i, 2] <- bestK
  print(projectName)
  print(table(projectSkmeans2to10[, bestK-1]))
  projectBestKSampleNum[[i]] <- table(projectSkmeans2to10[, bestK-1])
  
  projectSkmeans2to10Output <- file.path(getwd(), projectName, taxaLevel, "Skmeans", "Skmeans_2to10.txt")
  write.table(projectSkmeans2to10, file=projectSkmeans2to10Output, sep="\t", col.names=T, row.names=T, quote=F)
  
  aliquotInfo$Cluster[projectNameIndex] <- projectSkmeans2to10[, bestK-1]
  
}

table(aliquotInfo$Project_name)

## save best k found by mean silhoutte of skmeans in all project
projectBestKoutput <- file.path(getwd(), paste0(taxaLevel, "_Skmeans_best_K_by_SilhouetteMean.txt"))
write.table(projectBestK, file=projectBestKoutput, sep="\t", col.names=T, row.names=T, quote=F)

## save sample number of conditions
colnames(aliquotInfo)
?aggregate
resultSummary <- aggregate(Aliquot_barcode_miRNAseq ~ Project_name + Cluster + Sample_type,
                           data=aliquotInfo, FUN=length)
resultSummary <- resultSummary[order(resultSummary$Project_name, resultSummary$Cluster, resultSummary$Sample_type), ]
colnames(resultSummary) <- c("Project_name", "Cluster_in_each_project", "Sample_type", "Frequency")

resultSummaryOutput <- file.path(getwd(), paste0(taxaLevel, "_Skmeans_best_K_by_SilhouetteMean_SampleNum.txt"))
write.table(resultSummary, file=resultSummaryOutput, sep="\t", col.names=T, row.names=F, quote=F)

aliquotInfo[1:3, ]
aliquotInfoOutput <- file.path(getwd(), "TCGA_miRNAseq_and_mRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal.txt")
write.table(aliquotInfo, file=aliquotInfoOutput, sep="\t", col.names=T, row.names=F, quote=F)
