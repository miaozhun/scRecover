# Estimate dropout gene number in a cell
estDropoutNum <- function(sample = NULL, histCounts = NULL, depth = 100, return = "dropoutNum"){
  library(preseqR)
  if(!is.null(sample)){
    sample <- as.numeric(sample)
    sample <- ceiling(sample)
    histCounts <- cbind(as.numeric(names(table(sample))), as.numeric(unname(table(sample))))
    histCounts <- histCounts[-1,]
    colnames(histCounts) <- c("j", "n_j")
  }
  preseqR <- ztnb.rSAC(histCounts)
  geneNumPredict <- round(preseqR(depth))
  dropoutNum <- geneNumPredict - sum(histCounts[,2])
  if(return == "dropoutNum")
    return(dropoutNum)
  else
    return(geneNumPredict)
}




