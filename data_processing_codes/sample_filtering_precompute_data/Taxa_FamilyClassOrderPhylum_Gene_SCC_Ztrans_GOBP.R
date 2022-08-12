### do bacteria gene SCC calculation, SCC z-transformation and bacteria FGSEA (GO BP) 
### remember to modify the following variables

## fix output under following directory
workingDirRoot <- "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/"
setwd(workingDirRoot)
getwd()

## some cutoff settings
samplePercentCutoff <- 0.2 ## for bacteria co-abundance and bacteria-human rna SCC 
samplePercentLabel <- paste0("Percentage_", round(samplePercentCutoff*100, 2), "")

## load library
library("Hmisc")
library("fgsea")
library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")
#library("ComplexHeatmap")
#library("circlize")



## alliquot ids of miRNAseq and RNAseq with best k in Cluster column
aliquotInfoPath <- "/home/kaipu/microbiome_TCGA_miRNA/TCGA_projects/TCGA_miRNAseq_and_mRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal.txt"
file.exists(aliquotInfoPath)

aliquotInfo <- read.table(aliquotInfoPath, header=T, sep="\t", check.names=F)
dim(aliquotInfo)
aliquotInfo[1:2, ]



## color codes for projects (without LAML and GBM)
colorCode31Path <- "/home/kaipu/microbiome_TCGA_miRNA/PanCanAtlas/color_codes_without_LAML_GBM.txt"
file.exists(colorCode31Path)

colorCode31 <- read.table(colorCode31Path, header=T, check.names=F, comment.char="")
dim(colorCode31)


## human RNA expression downloaded from PanCanAtlas
rnaExpr <- readRDS("rnaExpr.rds")
dim(rnaExpr)
rnaExpr[1:40, 1:3]


## GO BP downloaded from Molecular Signatures Database
GOBPpath <- "/home/kaipu/microbiome_TCGA_miRNA/Gene_Ontology/c5.go.bp.v7.2.entrez.gmt"

GOBP <- gmtPathways(GOBPpath)
head(GOBP)
class(GOBP)
length(GOBP)



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
  rownames(bacteriaExpr)
  
  
  
  
  ## check if the order of bacteriaExpr and aliquotInfo are the same
  print(table(colnames(bacteriaExpr) == aliquotInfo$Aliquot_barcode_miRNAseq))
  #table(aliquotInfo$Sample_type)
  #sort(table(aliquotInfo$Project_name[which(aliquotInfo$Aliquot_barcode_RNAseq!="")]))
  
  
  ## remain samples both existed in Aliquot_barcode_miRNAseq and  Aliquot_barcode_RNAseq
  aliquotInfo[1, ]
  table(aliquotInfo$Aliquot_barcode_RNAseq=="")
  
  
  
  
  #cutoffList <- c(0.1, 0.15, 0.2, 0.25)
  #bacteriaSamPerCutoffSummary <- matrix(NA, nrow=nrow(colorCode31), ncol=4)
  #rownames(bacteriaSamPerCutoffSummary) <- colorCode31[, 1]
  #colnames(bacteriaSamPerCutoffSummary) <- c(paste0(">=", cutoffList))
  #pie(bacteriaSamPerCutoffSummary[, 3], col=as.character(colorCode31[, 5]))
  
  #bacteriaSamPerCutoffSummaryOutput <- file.path(getwd(), paste0(taxaLevel, "_bacteria_sample_percentage_cutoff_remained_taxa_num.txt"))
  #write.table(bacteriaSamPerCutoffSummary, file=bacteriaSamPerCutoffSummaryOutput, sep="\t", col.names=T, row.names=T, quote=F)
  
  
  
  remainedNumSummary <- matrix(NA, nrow=nrow(colorCode31), ncol=4)
  rownames(remainedNumSummary) <- colorCode31[, 1]
  colnames(remainedNumSummary) <- c("Sample_Num", "Sample_SCC_Num", paste0("Remain_", taxaLevel, "_Num"),"Remain_Gene_Num" )
  remainedNumSummary[1, ]
  
  
  i <- 1
  for (i in 1:nrow(colorCode31)) {
    projectName <- as.character(colorCode31[i, 1])
    print(projectName)
    remainedNumSummary[i, 1] <- length(which(c(aliquotInfo$Project_name %in% projectName)))
    
    projectNameBothIndex <- which(c(aliquotInfo$Project_name %in% projectName) & aliquotInfo$Aliquot_barcode_RNAseq!="")
    sampleSize <- length(projectNameBothIndex)
    remainedNumSummary[i, 2] <- sampleSize
    
    dir.create(file.path(getwd(), projectName, taxaLevel, "BacteriaGeneCorr", samplePercentLabel, "Gene_Ontology"), recursive=T, showWarnings=F)
    
    bacteriaExprBothProject <- bacteriaExpr[, projectNameBothIndex]
    dim(bacteriaExprBothProject)
    sum(colnames(bacteriaExprBothProject) == aliquotInfo$Aliquot_barcode_miRNAseq[projectNameBothIndex])
    
    rnaExprBothProjectIndex <- match(aliquotInfo$Aliquot_barcode_RNAseq[projectNameBothIndex],
                                     colnames(rnaExpr))
    length(rnaExprBothProjectIndex)
    rnaExprBothProject <- rnaExpr[, rnaExprBothProjectIndex]
    dim(rnaExprBothProject)
    sum(colnames(rnaExprBothProject) == aliquotInfo$Aliquot_barcode_RNAseq[projectNameBothIndex])
    
    
    #  ## check cutoff and remained taxa num
    #  k <- 1
    #  for (k in 1:length(cutoffList)) {
    #    bacteriaExprBothProjectSamplePercent <- rowSums(bacteriaExprBothProject>0)/ncol(bacteriaExprBothProject)
    #    bacteriaExprBothProjectSamplePercentPass <- bacteriaExprBothProject[which(bacteriaExprBothProjectSamplePercent>=cutoffList[k]), ]
    #    bacteriaSamPerCutoffSummary[i, k] <- nrow(bacteriaExprBothProjectSamplePercentPass)
    #  }
    #}
    #  {   
    
    bacteriaExprBothProjectSamplePercent <- rowSums(bacteriaExprBothProject>0)/ncol(bacteriaExprBothProject)
    plot(bacteriaExprBothProjectSamplePercent*100, xlab=taxaLevel, ylab="Percentage of sample",
         main=projectName, col=as.character(colorCode31[i, 5]), pch=19, cex=0.7)
    abline(h=20, lty="dashed", lwd=1)
    bacteriaExprBothProjectSamplePercentPass <- as.matrix(bacteriaExprBothProject[which(bacteriaExprBothProjectSamplePercent>=samplePercentCutoff), ])
    dim(bacteriaExprBothProjectSamplePercentPass)
    remainedNumSummary[i, 3] <- nrow(bacteriaExprBothProjectSamplePercentPass)
    
    
    rnaExprBothProjectSamplePercent <- rowSums(rnaExprBothProject>0)/ncol(rnaExprBothProject)
    rnaExprBothProjectSamplePercentPass <- as.matrix(rnaExprBothProject[which(rnaExprBothProjectSamplePercent>=samplePercentCutoff), ])
    dim(rnaExprBothProjectSamplePercentPass)
    remainedNumSummary[i, 4] <- nrow(rnaExprBothProjectSamplePercentPass)
    
    #?cor
    print("Calculate SCC...")
    bacteriaGeneSCC <- cor(t(rnaExprBothProjectSamplePercentPass), 
                           t(bacteriaExprBothProjectSamplePercentPass), 
                           method="spearman")
    dim(bacteriaGeneSCC)
    bacteriaGeneSCC[1:3, 1:3]
    range(unlist(bacteriaGeneSCC))
    
    
    ## handle gene id and gene symbol (some gene symbol are rename, update it by org.Hs.eg.db)
    geneInfo <- matrix(unlist(strsplit(rownames(bacteriaGeneSCC), split="|", fixed=T)), ncol=2, byrow=T)
    colnames(geneInfo) <- c("Gene_symbol", "Entrez_id")
    rownames(geneInfo) <- rownames(bacteriaGeneSCC)
    head(geneInfo)
    dim(geneInfo)
    
    #length(unique(geneInfo[, 1]))
    #nrow(geneInfo)
    #nrow(geneInfo[which(geneInfo[, 1]=="?"), ])
    #geneInfoMerge[which(geneInfoMerge[, 1]=="SLC35E2"), ]
    
    geneEntrez2symbol <- bitr(geneID=geneInfo[, 2], fromType="ENTREZID", 
                              toType="SYMBOL", OrgDb="org.Hs.eg.db", drop=F)
    geneEntrez2symbol[1:5, ]
    str(geneEntrez2symbol)
    dim(geneEntrez2symbol)
    
    table(geneInfo[, 2] == geneEntrez2symbol[, 1])
    geneInfoMerge <- cbind(geneInfo, geneEntrez2symbol)
    dim(geneInfoMerge)
    geneInfoMerge[1:3,]
    
    nrow(geneInfoMerge[which(geneInfo[, 1] != geneEntrez2symbol[, 2]), ])
    geneInfoMerge[which(is.na(geneEntrez2symbol$SYMBOL)), ]
    
    rownames(bacteriaGeneSCC) <- paste(geneInfoMerge[, 4], geneInfoMerge[, 3], sep="|")
    
    bacteriaGeneSCCoutput <- file.path(getwd(), projectName, taxaLevel, "BacteriaGeneCorr", samplePercentLabel,
                                       paste0(projectName, "_Bacteria_Gene_SCC_matrix.txt"))
    write.table(bacteriaGeneSCC, file=bacteriaGeneSCCoutput, sep="\t", col.names=T, row.names=T, quote=F)
    #bacteriaGeneSCC[1,1]
    
    #?rcorr
    print("Calculate P value of SCC...")
    bacteriaGeneSCCpVal <- t(sapply(1:nrow(rnaExprBothProjectSamplePercentPass), function(x) {
      sapply(1:nrow(bacteriaExprBothProjectSamplePercentPass), function(y) {
        rcorr(rnaExprBothProjectSamplePercentPass[x,], bacteriaExprBothProjectSamplePercentPass[y,],type="spearman")[[3]][1,2]
        #rcorr(rnaExprBothProjectSamplePercentPass[1,], bacteriaExprBothProjectSamplePercentPass[1,],type="spearman")[[1]][1,2]
      })
    }))
    dim(bacteriaGeneSCCpVal)
    bacteriaGeneSCCpVal[1:3, 1:3]
    colnames(bacteriaGeneSCCpVal) <- colnames(bacteriaGeneSCC)
    rownames(bacteriaGeneSCCpVal) <- rownames(bacteriaGeneSCC)
    
    bacteriaGeneSCCpValOutput <- file.path(getwd(), projectName, taxaLevel, "BacteriaGeneCorr", samplePercentLabel,
                                           paste0(projectName, "_Bacteria_Gene_SCC_Pvalue_matrix.txt"))
    write.table(bacteriaGeneSCCpVal, file=bacteriaGeneSCCpValOutput, sep="\t", col.names=T, row.names=T, quote=F)
    
    
    
    print("Apply Fisher's Z-transformation for SCC...")
    bacteriaGeneSCCzTransform <- (sqrt(sampleSize-3)/2)*log((1+bacteriaGeneSCC)/(1-bacteriaGeneSCC))
    dim(bacteriaGeneSCCzTransform)
    bacteriaGeneSCCzTransform[1:3, 1:3]
    range(unlist(bacteriaGeneSCCzTransform))
    #hist(unlist(bacteriaGeneSCC), breaks=100)
    #hist(unlist(bacteriaGeneSCCzTransform), breaks=100)
    #colSums(abs(bacteriaGeneSCCzTransform)>1.96)
    
    bacteriaGeneSCCzTransformOutput <- file.path(getwd(), projectName, taxaLevel, "BacteriaGeneCorr", samplePercentLabel,
                                                 paste0(projectName, "_Bacteria_Gene_SCC_Fisher_z_transform_matrix.txt"))
    write.table(bacteriaGeneSCCzTransform, file=bacteriaGeneSCCzTransformOutput, sep="\t", col.names=T, row.names=T, quote=F)
    
    
    
    
    print("Start GOBP analysis by FGSEA...")
    j <- 1
    ## GO BP analysis with FGSEA
    for (j in 1:ncol(bacteriaGeneSCCzTransform)) {
      bacteriaName <- colnames(bacteriaGeneSCCzTransform)[j]
      print(bacteriaName)
      orderDecreasing <- order(bacteriaGeneSCCzTransform[, j], decreasing=T)
      #head(orderDecreasing)
      length(orderDecreasing)
      
      preranks <- bacteriaGeneSCCzTransform[orderDecreasing, j]
      names(preranks) <- geneInfoMerge[orderDecreasing, 2]
      #length(preranks)
      head(preranks)
      tail(preranks)
      #length(which(table(preranks)==1)) #check tie numbers: few
      #sum(preranks>1.96)
      #sum(preranks< -1.96)
      
      
      GOBPresult <- fgsea(pathways=GOBP,
                          stats=preranks, 
                          nperm=1000, minSize=15, maxSize=500)
      
      GOBPresultPadjSort <- as.data.frame(GOBPresult[order(GOBPresult$padj, decreasing=F), ])
      head(GOBPresultPadjSort)
      dim(GOBPresultPadjSort)
      
      ## leadingEdge is list => convert to string and save in Gene column, and save it
      GOBPresultPadjSort$EntrezID <- NA
      GOBPresultPadjSort$GeneSymbol <- NA
      k <- 1
      for (k in 1:nrow(GOBPresultPadjSort)) {
        leadingEdgeGene <- unlist(GOBPresultPadjSort[k, 8])
        GOBPresultPadjSort$EntrezID[k] <- paste(leadingEdgeGene, collapse=",")
        
        entrez2symbolIndex <- match(leadingEdgeGene, geneInfoMerge[, 2])
        length(entrez2symbolIndex)
        
        GOBPresultPadjSort$GeneSymbol[k] <- paste(geneInfoMerge[entrez2symbolIndex, 4], collapse=",")
      }
      
      dir.create(file.path(getwd(), projectName, taxaLevel, "BacteriaGeneCorr", samplePercentLabel, "Gene_Ontology", bacteriaName), recursive=T, showWarnings=F)
      GOBPresultPadjSortOutput <- file.path(getwd(), projectName, taxaLevel, "BacteriaGeneCorr", samplePercentLabel, "Gene_Ontology", bacteriaName, "GO_BP.txt")
      write.table(GOBPresultPadjSort[, -8], file=GOBPresultPadjSortOutput, sep="\t", col.names=T, row.names=F, quote=F)
      
      collapsedPathways <- collapsePathways(GOBPresult[order(pval)][padj < 0.05], 
                                            GOBP, preranks)
      length(collapsedPathways$mainPathways)
      
      GOBPresultPadjSigSortOutput <- file.path(getwd(), projectName, taxaLevel, "BacteriaGeneCorr", samplePercentLabel, "Gene_Ontology", bacteriaName, "GO_BP_main_adjp_sig.txt")
      write.table(GOBPresultPadjSort[which(GOBPresultPadjSort$pathway %in% collapsedPathways$mainPathways), -8],
                  file=GOBPresultPadjSigSortOutput, sep="\t", col.names=T, row.names=F, quote=F)
      
    }
  }
  
  remainedNumSummaryOutput <- file.path(file.path(getwd(), 
                                                  paste0(taxaLevel, "_Gene_SCC_remained_sample_num_summary.txt")))
  write.table(remainedNumSummary, file=remainedNumSummaryOutput, sep="\t", col.names=T, row.names=T, quote=F)

  
}








