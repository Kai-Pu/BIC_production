## calculate diversities and evenness

setwd("~/microbiome_TCGA_miRNA/Diversity/")

## taxalevel
taxaLevel <- "Genus"
relativeAbun <- paste("~/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/TCGA_32_projects_merge_relative_abundance_nonFFPE_selected_PT_STN_GMPR_norm_", taxaLevel, ".txt", sep="")
diversityOutput <- paste("TCGA_32_projects_merge_diversity_", taxaLevel, ".txt", sep="")

## load library
library("vegan")
library('ggplot2')
library("ggpubr")
library("gridExtra")


## relative abundance
data <- read.table(relativeAbun, header=T, sep="\t", check.names=F)
dim(data)
str(data)
data[1:3, 1:3]


## load aliquot barcode table
aliquotInfo <- read.table("~/microbiome_TCGA_miRNA/FFPE_nonFFPE/nonFFPE_PT_STN/TCGA_miRNAseq_selected_samples_nonFFPE_Primary_Tumor_and_Solid_Tissue_Normal.txt",
                          header=T, sep="\t", check.names=F)
dim(aliquotInfo)
head(aliquotInfo, n=3)
str(aliquotInfo)

## project color code
colorCode <- read.table("~/microbiome_TCGA_miRNA/PanCanAtlas/color_codes_without_LAML.txt",
                        header=T, sep="\t", check.names=F, comment.char="")
dim(colorCode)
head(colorCode, n=2)

## check sample order
table(colnames(data)[2:ncol(data)] == aliquotInfo$Aliquot_barcode)


?diversity
## calculate shannon, simpson and simpson_reciprocal index, richness, evenness
#margin => 1: calculate by row; 2: calculate by column

## test 
#a <- matrix(c(1:9), ncol=3)
#a
#a_relatAbun <- t(t(a)/colSums(a))
#colSums(a_relatAbun)
#diversity(a, index="shannon", MARGIN=2)
#diversity(a_relatAbun, index="shannon", MARGIN=2)
#diversity(a, index="simpson", MARGIN=2)
#diversity(a_relatAbun, index="simpson", MARGIN=2)
#diversity(a, index="invsimpson", MARGIN=2)
#diversity(a_relatAbun, index="invsimpson", MARGIN=2)

#b <- a
#b[1:2,1] <- 0
#b[3, 3] <- 0
#colSums(b>0)

Richness <- colSums(data[, 2:ncol(data)]>0) #numbers of taxa in each samples
length(Richness)

ShannonIndex <- diversity(as.matrix(data[, 2:ncol(data)]), index="shannon", MARGIN=2)
length(ShannonIndex)

Evenness <- ShannonIndex/log(Richness)
length(Evenness)

SimpsonIndexDiv <- diversity(as.matrix(data[, 2:ncol(data)]), index="simpson", MARGIN=2)
length(SimpsonIndexDiv)

SimpsonRecipIndex <- diversity(as.matrix(data[, 2:ncol(data)]), index="invsimpson", MARGIN=2)
length(SimpsonRecipIndex)


## combine into data.frame
aliquotInfo_diversity <- cbind(aliquotInfo,
                               ShannonIndex=ShannonIndex,
                               SimpsonIndexDiv=SimpsonIndexDiv,
                               SimpsonRecipIndex=SimpsonRecipIndex,
                               Richness=Richness,
                               Evenness=Evenness)
dim(aliquotInfo_diversity)
head(aliquotInfo_diversity, n=2)
#table(rownames(aliquotInfo_diversity)==aliquotInfo_diversity$Aliquot_barcode)


## save al diversities into a single file
write.table(aliquotInfo_diversity, file=diversityOutput,
            sep="\t", col.names=T, row.names=F, quote=F)




## extract primay tumor samples for visualization
aliquotInfo_diversity_PT <- aliquotInfo_diversity[which(aliquotInfo_diversity$Sample_type=="Primary Tumor"), ]
dim(aliquotInfo_diversity_PT)
table(aliquotInfo_diversity_PT$Project_name)

project_PT_remain <- sort(as.character(unique(aliquotInfo_diversity_PT$Project_name)))



## filter rows of colorCode if there is 0 sample containing primary tumor
colorCodeRemain <- colorCode[-which(!(colorCode$Project_name %in% project_PT_remain)), ]


## find median of upper-quartile of each porjects for all diversities 
UpperQuartMatrix <- matrix(NA, nrow=length(project_PT_remain), ncol=5)
rownames(UpperQuartMatrix) <- project_PT_remain
colnames(UpperQuartMatrix) <- c("ShannonIndex", "SimpsonIndexDiv", "SimpsonRecipIndex", "Richness", "Evenness")
dim(UpperQuartMatrix)
head(UpperQuartMatrix, n=2)

i <- 1
aliquotInfo_diversity_PT[1, ]
for (i in 1:length(project_PT_remain)) {
    projectRowIndex <- which(aliquotInfo_diversity_PT$Project_name==project_PT_remain[i])
    length(projectRowIndex)
    ShannonIndexUpperQuart <- quantile(aliquotInfo_diversity_PT$ShannonIndex[projectRowIndex])[4]
    SimpsonIndexDivUpperQuart <- quantile(aliquotInfo_diversity_PT$SimpsonIndexDiv[projectRowIndex])[4]
    SimpsonRecipIndexUpperQuart <- quantile(aliquotInfo_diversity_PT$SimpsonRecipIndex[projectRowIndex])[4]
    RichnessUpperQuart <- quantile(aliquotInfo_diversity_PT$Richness[projectRowIndex])[4]
    EvennessUpperQuart <- quantile(aliquotInfo_diversity_PT$Evenness[projectRowIndex])[4]
    UpperQuartMatrix[i, ] <- c(ShannonIndexUpperQuart, SimpsonIndexDivUpperQuart, SimpsonRecipIndexUpperQuart, RichnessUpperQuart, EvennessUpperQuart)
}
?apply
UpperQuartMatrix <- as.data.frame(UpperQuartMatrix)
apply(UpperQuartMatrix, 2, median)



LowerQuartMatrix <- matrix(NA, nrow=length(project_PT_remain), ncol=5)
rownames(LowerQuartMatrix) <- project_PT_remain
colnames(LowerQuartMatrix) <- c("ShannonIndex", "SimpsonIndexDiv", "SimpsonRecipIndex", "Richness", "Evenness")
dim(LowerQuartMatrix)
head(LowerQuartMatrix, n=2)

i <- 1
aliquotInfo_diversity_PT[1, ]
for (i in 1:length(project_PT_remain)) {
  projectRowIndex <- which(aliquotInfo_diversity_PT$Project_name==project_PT_remain[i])
  length(projectRowIndex)
  ShannonIndexLowerQuart <- quantile(aliquotInfo_diversity_PT$ShannonIndex[projectRowIndex])[2]
  SimpsonIndexDivLowerQuart <- quantile(aliquotInfo_diversity_PT$SimpsonIndexDiv[projectRowIndex])[2]
  SimpsonRecipIndexLowerQuart <- quantile(aliquotInfo_diversity_PT$SimpsonRecipIndex[projectRowIndex])[2]
  RichnessLowerQuart <- quantile(aliquotInfo_diversity_PT$Richness[projectRowIndex])[2]
  EvennessLowerQuart <- quantile(aliquotInfo_diversity_PT$Evenness[projectRowIndex])[2]
  LowerQuartMatrix[i, ] <- c(ShannonIndexLowerQuart, SimpsonIndexDivLowerQuart, SimpsonRecipIndexLowerQuart, RichnessLowerQuart, EvennessLowerQuart)
}
?apply
LowerQuartMatrix <- as.data.frame(LowerQuartMatrix)
apply(LowerQuartMatrix, 2, median)

## boxplot
?svg
svg(filename=paste("Shannon_index_PT_", taxaLevel, ".svg", sep=""), 
    width=6, height=5)
ggplot(aliquotInfo_diversity_PT, 
       aes(x=Project_name, y=ShannonIndex, fill=Project_name)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
  xlab("") + ylab("Shannon Index") + 
  scale_fill_manual(values=as.vector(colorCodeRemain$`HexColor#`)) + 
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"),
        axis.title=element_text(size=10)) +
  geom_hline(aes(yintercept=median(UpperQuartMatrix$ShannonIndex)),
             colour="#FF0000", linetype="dashed") #+
#  geom_hline(aes(yintercept=median(LowerQuartMatrix$ShannonIndex)),
#             colour="#00B050", linetype="dashed")
dev.off()


svg(filename=paste("Simpson_index_Diversity_PT_", taxaLevel, ".svg", sep=""), 
    width=6, height=5)
ggplot(aliquotInfo_diversity_PT, 
       aes(x=Project_name, y=SimpsonIndexDiv, fill=Project_name)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
  xlab("") + ylab("Simpson's Index of Diversity") + 
  scale_fill_manual(values=as.vector(colorCodeRemain$`HexColor#`)) + 
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"),
        axis.title=element_text(size=10)) +
  geom_hline(aes(yintercept=median(UpperQuartMatrix$SimpsonIndexDiv)),
             colour="#FF0000", linetype="dashed") #+
#  geom_hline(aes(yintercept=median(LowerQuartMatrix$SimpsonIndexDiv)),
#             colour="#00B050", linetype="dashed")
dev.off()


svg(filename=paste("Simpson_Reciprocal_index_PT_", taxaLevel, ".svg", sep=""), 
    width=6, height=5)
ggplot(aliquotInfo_diversity_PT, 
       aes(x=Project_name, y=SimpsonRecipIndex, fill=Project_name)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
  xlab("") + ylab("Simpson's Reciprocal Index") + 
  scale_fill_manual(values=as.vector(colorCodeRemain$`HexColor#`)) + 
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"),
        axis.title=element_text(size=10)) +
  geom_hline(aes(yintercept=median(UpperQuartMatrix$SimpsonRecipIndex)),
             colour="#FF0000", linetype="dashed") #+
#  geom_hline(aes(yintercept=median(LowerQuartMatrix$SimpsonRecipIndex)),
#             colour="#00B050", linetype="dashed")
dev.off()


svg(filename=paste("Richness_PT_", taxaLevel, ".svg", sep=""), 
    width=6, height=5)
ggplot(aliquotInfo_diversity_PT, 
       aes(x=Project_name, y=Richness, fill=Project_name)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
  xlab("") + ylab("Richness") + 
  scale_fill_manual(values=as.vector(colorCodeRemain$`HexColor#`)) + 
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"),
        axis.title=element_text(size=10)) +
  geom_hline(aes(yintercept=median(UpperQuartMatrix$Richness)),
             colour="#FF0000", linetype="dashed") #+
#  geom_hline(aes(yintercept=median(LowerQuartMatrix$Richness)),
#             colour="#00B050", linetype="dashed")
dev.off()


svg(filename=paste("Richness_PT_ylim_0_400_", taxaLevel, ".svg", sep=""), 
    width=6, height=5)
ggplot(aliquotInfo_diversity_PT, 
       aes(x=Project_name, y=Richness, fill=Project_name)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
  xlab("") + ylab("Richness") + ylim(0, 400) + 
  scale_fill_manual(values=as.vector(colorCodeRemain$`HexColor#`)) + 
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"),
        axis.title=element_text(size=10)) +
  geom_hline(aes(yintercept=median(UpperQuartMatrix$Richness)),
             colour="#FF0000", linetype="dashed") #+
#  geom_hline(aes(yintercept=median(LowerQuartMatrix$Richness)),
#             colour="#00B050", linetype="dashed")
dev.off()



svg(filename=paste("Evenness_PT_", taxaLevel, ".svg", sep=""), 
    width=6, height=5)
ggplot(aliquotInfo_diversity_PT, 
       aes(x=Project_name, y=Evenness, fill=Project_name)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
  xlab("") + ylab("Evenness") + 
  scale_fill_manual(values=as.vector(colorCodeRemain$`HexColor#`)) + 
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=8, colour="black"),
        axis.text.y=element_text(size=8, colour="black"),
        axis.title=element_text(size=10)) +
#  geom_hline(aes(yintercept=median(UpperQuartMatrix$Evenness)),
#             colour="#FF0000", linetype="dashed") #+
  geom_hline(aes(yintercept=median(LowerQuartMatrix$Evenness)),
             colour="#00B050", linetype="dashed") 
dev.off()

aliquotInfo_diversity[1, ]





## extract projects whose sample size of both PT and STN are >= 10
table(aliquotInfo$Project_name)
project_PT_num <- table(aliquotInfo$Project_name[which(aliquotInfo$Sample_type=="Primary Tumor")])
project_STN_num <- table(aliquotInfo$Project_name[which(aliquotInfo$Sample_type=="Solid Tissue Normal")])
sessionInfo()


project_PTandSTN <- intersect(names(project_PT_num[which(project_PT_num>=10)]), names(project_STN_num[which(project_STN_num>=10)]))


## reorder the boxplot (STN, then PT)
aliquotInfo_diversity$Sample_type <- factor(aliquotInfo_diversity$Sample_type, 
                                  levels=c("Solid Tissue Normal", "Primary Tumor"))

str(aliquotInfo_diversity_PT$Sample_type)

pie(c(1,1), col=c("#30B7E7", "#F08080"))
colorPTandSTN <- c("#30B7E7", "#F08080")






## save all plots into a single figure
plot_list <- list()
for (i in 1:length(project_PTandSTN)) {
  project_df <- aliquotInfo_diversity[which(aliquotInfo_diversity$Project_name==project_PTandSTN[i]), ]
  plot_list[[i]] <- ggplot(project_df,
              aes(x=Sample_type, y=ShannonIndex, colour=Sample_type)) +
    geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
    xlab("") + ylab("Shannon index") + 
    scale_colour_manual(values=colorPTandSTN) +
    theme(legend.position="none", panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), panel.background=element_blank(), 
          axis.line=element_line(colour="black", size=0.25),
          axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5, size=8, colour="black"),
          axis.text.y=element_text(size=8, colour="black"),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12, hjust=0.5)) +
    ggtitle(project_PTandSTN[i]) +
    stat_compare_means(method="wilcox.test", size=3, 
                       aes(label=paste0("p = ", ..p.format..)))
}
svg(filename=paste("Shannon_index_STNvsPT_N10.svg", sep=""),
    width=15, height=9)
grid.arrange(grobs=plot_list,ncol=5)
dev.off()




plot_list <- list()
i <- 1
for (i in 1:length(project_PTandSTN)) {
  project_df <- aliquotInfo_diversity[which(aliquotInfo_diversity$Project_name==project_PTandSTN[i]), ]
  plot_list[[i]] <- ggplot(project_df,
                           aes(x=Sample_type, y=SimpsonIndexDiv, colour=Sample_type)) +
    geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
    xlab("") + ylab("Simpson's Index of Diversity") + xlim(fixYRange[1], fixYRange[2]) + 
    scale_colour_manual(values=colorPTandSTN) +
    theme(legend.position="none", panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), panel.background=element_blank(), 
          axis.line=element_line(colour="black", size=0.25),
          axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5, size=8, colour="black"),
          axis.text.y=element_text(size=8, colour="black"),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12, hjust=0.5)) +
    ggtitle(project_PTandSTN[i]) +
    stat_compare_means(method="wilcox.test", size=3, 
                       aes(label=paste0("p = ", ..p.format..)))
}
svg(filename=paste("Simpson_index_Diversity_STNvsPT_N10.svg", sep=""),
    width=15, height=9)
grid.arrange(grobs=plot_list,ncol=5)
dev.off()


plot_list <- list()
for (i in 1:length(project_PTandSTN)) {
  project_df <- aliquotInfo_diversity[which(aliquotInfo_diversity$Project_name==project_PTandSTN[i]), ]
  plot_list[[i]] <- ggplot(project_df,
                           aes(x=Sample_type, y=SimpsonRecipIndex, colour=Sample_type)) +
    geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
    xlab("") + ylab("Simpson's Reciprocal Index") + 
    scale_colour_manual(values=colorPTandSTN) +
    theme(legend.position="none", panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), panel.background=element_blank(), 
          axis.line=element_line(colour="black", size=0.25),
          axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5, size=8, colour="black"),
          axis.text.y=element_text(size=8, colour="black"),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12, hjust=0.5)) +
    ggtitle(project_PTandSTN[i]) +
    stat_compare_means(method="wilcox.test", size=3, 
                       aes(label=paste0("p = ", ..p.format..)))
}
svg(filename=paste("Simpson_Reciprocal_index_STNvsPT_N10.svg", sep=""),
    width=15, height=9)
grid.arrange(grobs=plot_list,ncol=5)
dev.off()


plot_list <- list()
for (i in 1:length(project_PTandSTN)) {
  project_df <- aliquotInfo_diversity[which(aliquotInfo_diversity$Project_name==project_PTandSTN[i]), ]
  plot_list[[i]] <- ggplot(project_df,
                           aes(x=Sample_type, y=Richness, colour=Sample_type)) +
    geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
    xlab("") + ylab("Richness") + 
    scale_colour_manual(values=colorPTandSTN) +
    theme(legend.position="none", panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), panel.background=element_blank(), 
          axis.line=element_line(colour="black", size=0.25),
          axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5, size=8, colour="black"),
          axis.text.y=element_text(size=8, colour="black"),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12, hjust=0.5)) +
    ggtitle(project_PTandSTN[i]) +
    stat_compare_means(method="wilcox.test", size=3, 
                       aes(label=paste0("p = ", ..p.format..)))
}
svg(filename=paste("Richness_STNvsPT_N10.svg", sep=""),
    width=15, height=9)
grid.arrange(grobs=plot_list,ncol=5)
dev.off()


plot_list <- list()
for (i in 1:length(project_PTandSTN)) {
  project_df <- aliquotInfo_diversity[which(aliquotInfo_diversity$Project_name==project_PTandSTN[i]), ]
  plot_list[[i]] <- ggplot(project_df,
                           aes(x=Sample_type, y=Evenness, colour=Sample_type)) +
    geom_boxplot(outlier.size=0.1, position=position_dodge(1), lwd=0.5)+ 
    xlab("") + ylab("Evenness") + 
    scale_colour_manual(values=colorPTandSTN) +
    theme(legend.position="none", panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), panel.background=element_blank(), 
          axis.line=element_line(colour="black", size=0.25),
          axis.text.x=element_text(angle=0, vjust=0.5, hjust=0.5, size=8, colour="black"),
          axis.text.y=element_text(size=8, colour="black"),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12, hjust=0.5)) +
    ggtitle(project_PTandSTN[i]) +
    stat_compare_means(method="wilcox.test", size=3, 
                       aes(label=paste0("p = ", ..p.format..)))
}
svg(filename=paste("Evenness_STNvsPT_N10.svg", sep=""),
    width=15, height=9)
grid.arrange(grobs=plot_list,ncol=5)
dev.off()




## combine 14 projects (PT and STN) into single boxplot 
aliquotInfo_diversity_PTandSTN <- aliquotInfo_diversity[which(aliquotInfo_diversity$Project_name %in% project_PTandSTN), ]
table(aliquotInfo_diversity_PTandSTN$Project_name)

aliquotInfo_diversity_PTandSTN$Sample_type <- factor(aliquotInfo_diversity_PTandSTN$Sample_type,
                                                     levels=c("Solid Tissue Normal",
                                                            "Primary Tumor"))
table(aliquotInfo_diversity_PTandSTN$Sample_type)


svg(filename=paste("Shannon_index_STNvsPT_N10_merge.svg", sep=""),
    width=12, height=8)
ggplot(aliquotInfo_diversity_PTandSTN,
       aes(x=Project_name, y=ShannonIndex, colour=Sample_type)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(0.9), lwd=0.5)+ 
  xlab("") + ylab("Shannon index") + 
  scale_colour_manual(values=colorPTandSTN) +
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=10, colour="black"),
        axis.text.y=element_text(size=10, colour="black"),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12, hjust=0.5)) +
  stat_compare_means(method="wilcox.test", size=4, label="p.signif",
                     aes(group=Sample_type))
dev.off()



svg(filename=paste("Simpson_index_Diversity_STNvsPT_N10_merge.svg", sep=""),
    width=12, height=8)
ggplot(aliquotInfo_diversity_PTandSTN,
       aes(x=Project_name, y=SimpsonIndexDiv, colour=Sample_type)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(0.9), lwd=0.5)+ 
  xlab("") + ylab("Simpson's Index of Diversity") + 
  scale_colour_manual(values=colorPTandSTN) +
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=10, colour="black"),
        axis.text.y=element_text(size=10, colour="black"),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12, hjust=0.5)) +
  stat_compare_means(method="wilcox.test", size=4, label="p.signif",
                     aes(group=Sample_type))
dev.off()



svg(filename=paste("Simpson_Reciprocal_index_STNvsPT_N10_merge.svg", sep=""),
    width=12, height=8)
ggplot(aliquotInfo_diversity_PTandSTN,
       aes(x=Project_name, y=SimpsonRecipIndex, colour=Sample_type)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(0.9), lwd=0.5)+ 
  xlab("") + ylab("Simpson's Reciprocal Index") + 
  scale_colour_manual(values=colorPTandSTN) +
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=10, colour="black"),
        axis.text.y=element_text(size=10, colour="black"),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12, hjust=0.5)) +
  stat_compare_means(method="wilcox.test", size=4, label="p.signif",
                     aes(group=Sample_type))
dev.off()



svg(filename=paste("Richness_STNvsPT_N10_merge.svg", sep=""),
    width=12, height=8)
ggplot(aliquotInfo_diversity_PTandSTN,
       aes(x=Project_name, y=Richness, colour=Sample_type)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(0.9), lwd=0.5)+ 
  xlab("") + ylab("Richness") + 
  scale_colour_manual(values=colorPTandSTN) +
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=10, colour="black"),
        axis.text.y=element_text(size=10, colour="black"),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12, hjust=0.5)) +
  stat_compare_means(method="wilcox.test", size=4, label="p.signif",
                     aes(group=Sample_type))
dev.off()



svg(filename=paste("Evenness_STNvsPT_N10_merge.svg", sep=""),
    width=12, height=8)
ggplot(aliquotInfo_diversity_PTandSTN,
       aes(x=Project_name, y=Evenness, colour=Sample_type)) +
  geom_boxplot(outlier.size=0.1, position=position_dodge(0.9), lwd=0.5)+ 
  xlab("") + ylab("Evenness") + 
  scale_colour_manual(values=colorPTandSTN) +
  theme(legend.position="none", panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), panel.background=element_blank(), 
        axis.line=element_line(colour="black", size=0.25),
        axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5, size=10, colour="black"),
        axis.text.y=element_text(size=10, colour="black"),
        axis.title=element_text(size=12),
        plot.title = element_text(size=12, hjust=0.5)) +
  stat_compare_means(method="wilcox.test", size=4, label="p.signif",
                     aes(group=Sample_type))
dev.off()



# save entire list of environment objects
save.image(file="Diversity2.RData")
quit(save='no')
