# Normalization function
normalization <- function(counts){
  # Filter genes with all zero counts
  counts <- as.matrix(counts)
  counts <- ceiling(counts)
  sampleNum <- ncol(counts)
  if(any(rowSums(counts) == 0))
    message("Removing ", sum(rowSums(counts) == 0), " rows of genes with all zero counts")
  counts_NAZ <- counts[rowSums(counts) != 0,]
  geneNum_NAZ <- nrow(counts_NAZ)
  
  # Normalization
  GEOmean <- rep(NA,geneNum_NAZ)
  for (i in 1:geneNum_NAZ)
  {
    gene_NZ <- counts_NAZ[i,counts_NAZ[i,] > 0]
    GEOmean[i] <- exp(sum(log(gene_NZ), na.rm=TRUE) / length(gene_NZ))
  }
  S <- rep(NA, sampleNum)
  counts_norm <- counts_NAZ
  for (j in 1:sampleNum)
  {
    sample_j <- counts_NAZ[,j]/GEOmean
    S[j] <- median(sample_j[which(sample_j != 0)])
    counts_norm[,j] <- counts_NAZ[,j]/S[j]
  }
  counts_norm <- ceiling(counts_norm)
  return(counts_norm)
}




