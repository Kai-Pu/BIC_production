setwd("/home/kaipu/microbiome_TCGA_miRNA/Count_unique/Merge_count_nonFFPE/")

taxaLevel <- "Order"
countfileName <- paste("TCGA_32_projects_merge_counts_nonFFPE_selected_PT_STN_", taxaLevel, ".txt", sep="")

output_gmpr_norm_count <- gsub("_STN_", "_STN_GMPR_norm_", countfileName)
output_gmpr_norm_relative_abundance <- gsub("counts", "relative_abundance", output_gmpr_norm_count)
output_gmpr_norm_NAs_samples <- paste("GMPR_size_factor_NA_aliquot_barcodes_", taxaLevel, ".txt", sep="")
countfileName
output_gmpr_norm_count
output_gmpr_norm_relative_abundance
output_gmpr_norm_NAs_samples

require(matrixStats)


GMPR <- function (comm, intersect.no = 10, ct.min = 1, trace = TRUE) {
  # Computes the GMPR size factor
  #
  # Args:
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
  #   ct.min: the minimum number of counts required to calculate ratios
  #
  # Returns:
  #   a vector of the size factors with attribute 'NSS'. Samples with distinct sets of features will be output as NA.
  #         NSS:   number of samples with significant sharing (> intersect.no) including itself
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  if (trace) cat('Begin GMPR size factor calculation ...\n')
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {		
    if (i %% 50 == 0) {
      cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))		
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n', 
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  if (trace) cat('Completed!\n')
  if (trace) cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  names(gmpr) <- names(comm.no) <- colnames(comm)
  attr(gmpr, 'NSS') <- comm.no
  return(gmpr)
}

# load data
data <- read.table(countfileName, header=TRUE, sep='\t', quote="\"", check.names=F)
dim(data)
data[1:3, 1:5]

# do GMPR normalization
gmpr.size.factor <- GMPR(as.matrix(data[, 2:ncol(data)]), intersect.no=10) 

summary(is.na(gmpr.size.factor))
length(gmpr.size.factor)
head(gmpr.size.factor)


gmpr.size.factor2 <- GMPR(as.matrix(data[, 2:ncol(data)]), intersect.no=9) 
summary(is.na(gmpr.size.factor2))
gmpr.size.factor3 <- GMPR(as.matrix(data[, 2:ncol(data)]), intersect.no=8) 
summary(is.na(gmpr.size.factor3))
gmpr.size.factor4 <- GMPR(as.matrix(data[, 2:ncol(data)]), intersect.no=7) 
summary(is.na(gmpr.size.factor4))
gmpr.size.factor5 <- GMPR(as.matrix(data[, 2:ncol(data)]), intersect.no=6) 
summary(is.na(gmpr.size.factor5))
gmpr.size.factor6 <- GMPR(as.matrix(data[, 2:ncol(data)]), intersect.no=5) 
summary(is.na(gmpr.size.factor6))
gmpr.size.factor7 <- GMPR(as.matrix(data[, 2:ncol(data)]), intersect.no=4) 
summary(is.na(gmpr.size.factor7))
gmpr.size.factor8 <- GMPR(as.matrix(data[, 2:ncol(data)]), intersect.no=3) 
summary(is.na(gmpr.size.factor8))

summary(is.na(gmpr.size.factor))
summary(is.na(gmpr.size.factor2))
summary(is.na(gmpr.size.factor3))
summary(is.na(gmpr.size.factor4))
summary(is.na(gmpr.size.factor5))
summary(is.na(gmpr.size.factor6))
#summary(is.na(gmpr.size.factor7))
#summary(is.na(gmpr.size.factor8))



## Remember to modify these!!!
## use gmpr.size.factor for normalization (intersect.no=5) 
length(gmpr.size.factor6)
head(gmpr.size.factor6)




## select gmpr.size.factors which are all not NAs
data_gmpr <- t(t(as.matrix(data[, 2:ncol(data)]))/gmpr.size.factor6)
dim(data_gmpr)
gmpr.size.factor[1:4]
data_gmpr[1:3, 1:5]
data[1:3, 1:5]
rownames(data_gmpr) <- data[, 1]

## cehck NA numbers 
sum(colSums(is.na(data_gmpr)))
nrow(data_gmpr) * sum(colSums(is.na(data_gmpr)))


# turn normalized count into relative abundance
col_sum <- colSums(data_gmpr)
length(col_sum)
table(is.na(col_sum))

data_gmpr_relative_abundance <- t(t(data_gmpr)/col_sum)
head(colSums(data_gmpr_relative_abundance))
tail(colSums(data_gmpr_relative_abundance))

table(round(colSums(data_gmpr_relative_abundance), 5)==1)


## save normalized count and relative abundance to txt files
table(colnames(data)[2:ncol(data)]==colnames(data_gmpr))

data_gmpr_df <- cbind(data[, 1], as.data.frame(data_gmpr))
colnames(data_gmpr_df) <- colnames(data)
data_gmpr_df[1:3, 1:3]
dim(data_gmpr_df)

row_sum <- rowSums(data_gmpr>0)
hist(row_sum)
row_sum[which(row_sum==0)]

write.table(data_gmpr_df, file=output_gmpr_norm_count, sep="\t",
            row.names=F, col.names=T, quote=F)


data_gmpr_relative_abundance_df <- cbind(data[, 1], as.data.frame(data_gmpr_relative_abundance))
colnames(data_gmpr_relative_abundance_df) <- colnames(data)
data_gmpr_relative_abundance_df[1:5, 1:5]
dim(data_gmpr_relative_abundance_df)

write.table(data_gmpr_relative_abundance_df, file=output_gmpr_norm_relative_abundance, sep="\t",
            row.names=F, col.names=T, quote=F)


## output aliqout barcodes whose normalized factors were fond to be NA
aliquot_barcode_NAs <- names(gmpr.size.factor6)[which(is.na(gmpr.size.factor6))]
write.table(aliquot_barcode_NAs, file=output_gmpr_norm_NAs_samples, 
            sep="\n", col.names=F, row.names=F, quote=F)


# save entire list of environment objects
save.image(file=paste("GMPR_norm_", taxaLevel, ".RData", sep=""))
quit(save='no')
