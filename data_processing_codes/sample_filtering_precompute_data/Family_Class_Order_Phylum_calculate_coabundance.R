## calculate correlation between Family/ Class/ Order/ Phylum

## spherical k-means for clustering by Phylum level
skmeansClustBy <- "Phylum" 

## fix output under following directory
workingDirRoot <- "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/"
setwd(workingDirRoot)
getwd()

## some cutoff settings
samplePercentCutoff <- 0.2 ## for bacteria co-abundance and bacteria-human rna SCC 
samplePercentLabel <- paste0("Percentage_", round(samplePercentCutoff*100, 2), "_")

# bootstrapTimes <- 100
# bootstrapTimes <- 1000
bootstrapTimes <- 10000

## load library
library("SpiecEasi")
library("pheatmap")

## aliquot ids of miRNAseq and RNAseq with best k in Cluster column
aliquotInfoPath <- "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/TCGA_miRNAseq_and_mRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal_v2.txt"
file.exists(aliquotInfoPath)

aliquotInfo <- read.table(aliquotInfoPath, header=T, sep="\t", check.names=F)
dim(aliquotInfo)
aliquotInfo[1:2, ]



## color codes for projects (without LAML and GBM)
colorCode31Path <- "/home/kaipu/microbiome_TCGA_miRNA/PanCanAtlas/color_codes_without_LAML_GBM.txt"
file.exists(colorCode31Path)

colorCode31 <- read.table(colorCode31Path, header=T, check.names=F, comment.char="")
dim(colorCode31)



taxaLevelList <- c("Family", "Class", "Order", "Phylum")


for (i in 1:length(taxaLevelList)) {
  taxaLevel <- taxaLevelList[i]
  print(taxaLevel)
  
  
  ## bacteria relative abundance
  bacteriaExprPath <- paste0("/home/kaipu/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/TCGA_32_projects_merge_relative_abundance_nonFFPE_selected_PT_STN_GMPR_norm_", taxaLevel, ".txt")
  file.exists(bacteriaExprPath)
  
  bacteriaExpr <- read.table(bacteriaExprPath, header=T, sep="\t", check.names=F, row.names=1)
  dim(bacteriaExpr)
  bacteriaExpr[1:3, 1:3]
  rownames(bacteriaExpr)[1:2]
  
  
  ## bacteria count
  bacteriaExprCountPath <- paste0("/home/kaipu/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/TCGA_32_projects_merge_counts_nonFFPE_selected_PT_STN_GMPR_norm_", taxaLevel, ".txt")
  file.exists(bacteriaExprCountPath)
  
  bacteriaExprCount <- read.table(bacteriaExprCountPath, header=T, sep="\t", check.names=F, row.names=1)
  dim(bacteriaExprCount)
  bacteriaExprCount[1:3, 1:3]
  rownames(bacteriaExprCount)[1:2]
  
  
  
  ## check if the order of bacteriaExpr and aliqoutInfo are the same
  print(table(colnames(bacteriaExpr) == aliquotInfo$Aliquot_barcode_miRNAseq))
  print(table(colnames(bacteriaExprCount) == aliquotInfo$Aliquot_barcode_miRNAseq))
  
  
  
  
  resultSummary <- matrix(NA, nrow=nrow(colorCode31), ncol=3)
  colnames(resultSummary) <- c("bacteriaNumRemain", "Sparsity (%)", "Neff")
  rownames(resultSummary) <- colorCode31[, 1]
  
  #i <- 9
  for (i in 1:nrow(colorCode31)) {
    projectName <- as.character(colorCode31[i, 1])
    print(projectName)
    
    dir.create(file.path(getwd(), projectName, taxaLevel, "bacteria_Coabundance", 
                         paste0("bootstrapTimes", bootstrapTimes)), recursive=T, showWarnings=F)
    
    
    projectIndex <- which(aliquotInfo$Project_name == projectName)
    length(projectIndex)
    
    bacteriaExprProject <- bacteriaExpr[, projectIndex]
    print(dim(bacteriaExprProject))
    bacteriaExprProject[1:3, 1:3]
    
    bacteriaExprCountProject <- bacteriaExprCount[, projectIndex]
    print(dim(bacteriaExprCountProject))
    bacteriaExprCountProject[1:3, 1:3]
    
    
    aliquotInfoProject <- aliquotInfo[projectIndex, ]
    print(dim(aliquotInfoProject))
    aliquotInfoProject[1:3, ]
    
    ## check again the order of bacteriaExprProject and aliquotInfoProject are the same
    table(colnames(bacteriaExprProject) ==  aliquotInfoProject$Aliquot_barcode_miRNAseq)
    table(colnames(bacteriaExprCountProject) ==  aliquotInfoProject$Aliquot_barcode_miRNAseq)
    
    
    
    
    
    
    ## calculate sparscc for genus expressed >= samplePercentCutoff
    bacteriaExistProject <- rowSums(bacteriaExprProject>0)
    bacteriaExistPercentProject <- bacteriaExistProject/ncol(bacteriaExprProject)
    length(bacteriaExistPercentProject)
    #hist(bacteriaExistPercentProject)
    
    table(bacteriaExistPercentProject >= samplePercentCutoff)
    
    bacteriaExprProjectPass <- bacteriaExprProject[which(bacteriaExistPercentProject>=samplePercentCutoff), ]
    print(dim(bacteriaExprProjectPass))
    
    
    bacteriaExprCountProjectPass <- bacteriaExprCountProject[which(bacteriaExistPercentProject>=samplePercentCutoff), ]
    print(dim(bacteriaExprCountProjectPass))
    
    
    ## check entropy effective number (Neff: average of entropy effective number)
    # tmp <- as.matrix(bacteriaExprProject*log(bacteriaExprProject))
    # dim(tmp)
    # tmp[is.nan(tmp)] <- 0
    # H <- -colSums(tmp)
    # length(H)
    # entEffNum <- exp(H)
    # entEffNumMean <- mean(entEffNum)
    
    tmp2 <- as.matrix(bacteriaExprProjectPass*log(bacteriaExprProjectPass))
    dim(tmp2)
    tmp2[1:5, 1:5]
    tmp2[is.nan(tmp2)] <- 0
    H2 <- -colSums(tmp2)
    entEffNum2 <- exp(H2)
    entEffNumMean2 <- mean(entEffNum2)
    
    
    sparsity <- sum(unlist(bacteriaExprProject==0))/(ncol(bacteriaExprProject)*nrow(bacteriaExprProject))*100
    passSparsity <- sum(unlist(bacteriaExprProjectPass==0))/(ncol(bacteriaExprProjectPass)*nrow(bacteriaExprProjectPass))*100
    
    
    resultSummary[i, 1] <- nrow(bacteriaExprCountProjectPass)
    resultSummary[i, 2] <- passSparsity
    resultSummary[i, 3] <- entEffNumMean2
    
    # ?sparcc
    # sparccResult <- sparcc(t(bacteriaExprProjectPass))
    # dim(sparccResult$Cor)
    # dim(sparccResult$Cov)
    # 
    # sparccResult2 <- sparcc(t(bacteriaExprCountProjectPass))
    # dim(sparccResult2$Cor)
    # dim(sparccResult2$Cov)
    # hist(sparccResult2$Cor, breaks=100, xlim=c(-0.5, 0.5))
    # heatmap(sparccResult2$Cor)
    
    #sparccbootResult <- sparccboot(t(bacteriaExprProjectPass), R=100, ncpus=6)
    #summary(sparccbootResult)
    #sparccbootResultpVal <- matrix(NA, nrow=nrow(bacteriaExprProjectPass), ncol=nrow(bacteriaExprProjectPass))
    #sparccbootResultCor <- matrix(NA, nrow=nrow(bacteriaExprProjectPass), ncol=nrow(bacteriaExprProjectPass))
    #sparccbootResultpVal[lower.tri(sparccbootResultpVal, diag=F)] <- pval.sparccboot(sparccbootResult)$pvals
    #sparccbootResultCor[lower.tri(sparccbootResultCor, diag=F)] <- pval.sparccboot(sparccbootResult)$cors
    #table(sparccbootResultpVal<0.05)
    
    #hist(sparccbootResultCor, breaks=100, xlab="sparCC correlation coefficient", main="Relative abundance as input")
    #length(sparccbootResultCor)
    #(nrow(bacteriaExprProjectPass)*nrow(bacteriaExprProjectPass)-nrow(bacteriaExprProjectPass))/2
    
    #heatmap(sparccbootResultCor, na.rm=T)
    #table(is.nan(sparccbootResultCor))
    #table(is.na(sparccbootResultCor))
    
    
    ?sparccboot
    sparccbootResult2 <- sparccboot(t(bacteriaExprCountProjectPass), R=bootstrapTimes, ncpus=6)
    #summary(sparccbootResult2)
    length(pval.sparccboot(sparccbootResult2)$pvals)
    length(pval.sparccboot(sparccbootResult2)$cors)
    
    sparccbootResult2pVal <- matrix(NA, nrow=nrow(bacteriaExprCountProjectPass), ncol=nrow(bacteriaExprCountProjectPass))
    sparccbootResult2Cor <- matrix(NA, nrow=nrow(bacteriaExprCountProjectPass), ncol=nrow(bacteriaExprCountProjectPass))
    sparccbootResult2pVal[lower.tri(sparccbootResult2pVal, diag=F)] <- pval.sparccboot(sparccbootResult2)$pvals
    sparccbootResult2Cor[lower.tri(sparccbootResult2Cor, diag=F)] <- pval.sparccboot(sparccbootResult2)$cors
    print(table(sparccbootResult2pVal<0.05))
    #  print(table(sparccbootResult2pVal<0.01))
    #  print(table(sparccbootResult2pVal<0.001))
    print(table(sparccbootResult2pVal<0.0001))
    colnames(sparccbootResult2Cor) <- rownames(bacteriaExprCountProjectPass)
    rownames(sparccbootResult2Cor) <- rownames(bacteriaExprCountProjectPass)
    colnames(sparccbootResult2pVal) <- rownames(bacteriaExprCountProjectPass)
    rownames(sparccbootResult2pVal) <- rownames(bacteriaExprCountProjectPass)
    #table(sparccbootResult2$t0==pval.sparccboot(sparccbootResult2)$cors)
    
    sparccbootResult2Cor[1:5, 1:4]
    
    
    sparccbootResult2pValT <- t(sparccbootResult2pVal)
    
    sparccbootResult2pValFull <- sparccbootResult2pVal
    sparccbootResult2pValFull[!is.na(sparccbootResult2pValT)] <- sparccbootResult2pValT[!is.na(sparccbootResult2pValT)]
    sparccbootResult2pValFull[1:3, 1:3]
    table(sparccbootResult2pValFull==t(sparccbootResult2pValFull))
    #sparccbootResult2pValFull[is.na(sparccbootResult2pValFull)] <- 1
    
    
    sparccbootResult2CorT <- t(sparccbootResult2Cor)
    
    sparccbootResult2CorFull <- sparccbootResult2Cor
    sparccbootResult2CorFull[!is.na(sparccbootResult2CorT)] <- sparccbootResult2CorT[!is.na(sparccbootResult2CorT)]
    sparccbootResult2CorFull[1:3, 1:3]
    table(sparccbootResult2CorFull==t(sparccbootResult2CorFull))
    #sparccbootResult2CorFull[is.na(sparccbootResult2CorFull)] <- 1
    
    
    
    #min(sparccbootResult2CorFull)
    #range(as.vector(sparccbootResult2CorFull))
    
    corMax <- max(abs(range(sparccbootResult2CorFull, na.rm=T)))
    col_breaks <- seq(-corMax, corMax, by=0.01)
    col_heatmap <- colorRampPalette(c("#30A9DE", "white", "#E53A40"))(length(col_breaks))
    p <- pheatmap(sparccbootResult2CorFull,
                  color=col_heatmap, fontface="bold", fontsize=5, breaks=col_breaks)
    
    heatmapOutput <- file.path(getwd(), projectName, taxaLevel, "bacteria_Coabundance", paste0("bootstrapTimes", bootstrapTimes), 
                               paste0(projectName, "_heatmap_bacteria_coabundance_sparcc_", samplePercentLabel, taxaLevel, ".svg"))
    svg(heatmapOutput, width=10, height=10)
    print(p)
    dev.off()
    #dev.new()
    
    histOutput <- file.path(getwd(), projectName, taxaLevel, "bacteria_Coabundance", paste0("bootstrapTimes", bootstrapTimes),
                            paste0(projectName, "_hist_bacteria_coabundance_sparcc_correlation_", samplePercentLabel, taxaLevel, ".svg"))
    #svg(histOutput, width=8, height=4)
    svg(histOutput, width=4, height=4)
    #par(mfrow=c(1,2))
    hist(sparccbootResult2Cor, breaks=100, xlab="sparCC correlation coefficient", main=projectName)
    #hist(pval.sparccboot(sparccbootResult2)$cors, breaks=100)
    #hist(sparccbootResult2pVal[which(sparccbootResult2pVal<0.05)], breaks=100)
    #hist(sparccbootResult2Cor[which(sparccbootResult2pVal<0.05)], breaks=100,
    #     xlab="sparCC correlation coefficient", main=paste0(projectName, " (bootstrap, p< 0.05)"))
    dev.off()
    #plot(sparccbootResult2Cor, -log10(sparccbootResult2pVal))
    #abline(h=-log10(0.05), col="green")
    #hist(pval.sparccboot(sparccbootResult2)$cors[which(pval.sparccboot(sparccbootResult2)$pvals<0.05)], breaks=100)
    #length(sparccbootResult2Cor)
    #(nrow(bacteriaExprCountProjectPass)*nrow(bacteriaExprCountProjectPass)-nrow(bacteriaExprCountProjectPass))/2
    #(nrow(bacteriaExprCountProjectPass)*nrow(bacteriaExprCountProjectPass)-nrow(bacteriaExprCountProjectPass))/2 + nrow(bacteriaExprCountProjectPass)
    
    #summary(as.vector(sparccbootResult2Cor))
    #heatmap(sparccbootResult2Cor, na.rm=T)
    #table(is.nan(sparccbootResult2Cor))
    #table(is.na(sparccbootResult2Cor))
    #table(is.na(pval.sparccboot(sparccbootResult2)$cors))
    
    # pargs=list(seed=666,ncores=6)
    # spiec <- spiec.easi(as.matrix(t(bacteriaExprCountProjectPass)), 
    #                     method='glasso', nlambda=30, pulsar.params=pargs)
    # summary(spiec)
    # cor <- cov2cor(apply(getOptCov(spiec), 2, as.numeric))
    # dim(cor)
    # heatmap(cor)
    # hist(cor, breaks=100, xlim=c(-0.5, 0.5))
    # row.names(cor) <- rownames(bacteriaExprCountProjectPass)
    # colnames(cor) <- rownames(bacteriaExprCountProjectPass)
    
    ## test convert a vector into a lower triangle of the  matrix
    # a <- 1:6
    # b <- matrix(NA, nrow=4, ncol=4)
    # b[lower.tri(b, diag=FALSE)] <- a
    # b <- t(b)
    # b
    
    
    #?pval.sparccboot
    ## save corr and p
    
    corOutput <- file.path(getwd(), projectName, taxaLevel, "bacteria_Coabundance", paste0("bootstrapTimes", bootstrapTimes), 
                           paste0(projectName, "_bacteria_coabundance_sparcc_correlation_", samplePercentLabel, taxaLevel, ".txt"))
    write.table(sparccbootResult2CorFull, file=corOutput, sep="\t", col.names=T, row.names=T, quote=F)
    
    pValOutput <- file.path(getwd(), projectName, taxaLevel, "bacteria_Coabundance", paste0("bootstrapTimes", bootstrapTimes), 
                            paste0(projectName, "_bacteria_coabundance_sparcc_p_Val_", samplePercentLabel, taxaLevel, ".txt"))
    write.table(sparccbootResult2pValFull, file=pValOutput, sep="\t", col.names=T, row.names=T, quote=F)
    
  }
  
  
  summaryOutput <- file.path(getwd(), 
                             paste0(taxaLevel, "_bacteria_coabundance_summary.txt"))
  write.table(resultSummary, file=summaryOutput, sep="\t", col.names=T, row.names=T, quote=F)
  
  resultSummary
}








