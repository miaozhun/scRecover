# Estimate dropout gene number in a cell
estDropoutNum <- function(sample, depth){
  sample <- as.numeric(sample)
  sample <- ceiling(sample)
  histCounts <- cbind(as.numeric(names(table(sample))), as.numeric(unname(table(sample))))
  histCounts <- histCounts[-1,]
  colnames(histCounts) <- c("j", "n_j")
  preseqR <- ztnb.rSAC(histCounts)
  geneNumPredict <- round(preseqR(depth))
  dropoutNum <- geneNumPredict - sum(sample != 0)
  return(dropoutNum)
}




