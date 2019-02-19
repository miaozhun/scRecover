#' estDropoutNum: Estimate dropout gene number in a cell
#'
#' This function is used to estimate dropout gene number in a cell for single-cell RNA-seq (scRNA-seq) data. It takes a non-negative vector of scRNA-seq raw read counts of a cell as input.
#'
#' @param sample A vector of a cell's raw read counts for each gene.
#' @param depth Relative sequencing depth to be predicted compared with initial sample depth, should between 2-100, default is 20.
#' @param histCounts Optional. Only needed when \code{sample} is blank or \code{sample = NULL}. A histogram table of raw read counts for the cell.
#' @param return A character for choosing the return value type of the function. Either be "dropoutNum" (default) for dropout gene number or "geneNumPredict" for all expressed gene number predicted.
#' @return
#' The dropout gene number (or all expressed gene number) predicted in a cell.
#'
#' @author Zhun Miao.
#' @seealso
#' \code{\link{scRecover}}, for the imputation of single-cell RNA-seq data.
#'
#' \code{\link{scRecoverTest}}, a test dataset for scRecover.
#'
#' @examples
#' # Load test data for scRecover
#' data(scRecoverTest)
#'
#' # Estimate dropout gene number in a cell
#' estDropoutNum(sample = counts[,1])
#'
#' # Estimate all expressed gene number in a cell
#' estDropoutNum(sample = counts[,1], return = "geneNumPredict")
#'
#'
#' @import stats
#' @importFrom utils read.csv write.csv
#' @importFrom parallel detectCores
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb fitdistr
#' @importFrom pscl zeroinfl
#' @importFrom bbmle mle2
#' @importFrom gamlss gamlssML
#' @importFrom preseqR ztnb.rSAC
#' @importFrom scImpute scimpute
#' @importFrom SAVER saver
#' @importFrom Rmagic magic
#' @importFrom BiocParallel bpparam bplapply
#' @export



# Estimate dropout gene number in a cell
estDropoutNum <- function(sample = NULL, depth = 20, histCounts = NULL, return = "dropoutNum"){
  # Invalid input control
  if(is.null(sample) & is.null(histCounts))
    stop("One of 'sample' and 'histCounts' must be specified")
  if(!is.null(sample) & !is.null(histCounts))
    stop("Only one of 'sample' and 'histCounts' should be specified")

  if(!is.vector(sample) & !is.data.frame(sample) & !is.integer(sample) & !is.numeric(sample) & !is.double(sample))
    stop("Wrong data type of 'sample'")
  if(sum(is.na(sample)) > 0)
    stop("NA detected in 'sample'")
  if(sum(sample < 0) > 0)
    stop("Negative value detected in 'sample'")
  if(all(sample == 0))
    stop("All elements of 'sample' are zero")

  if(!is.numeric(depth) & !is.integer(depth))
    stop("Data type of 'depth' is not numeric or integer")
  if(length(depth) != 1)
    stop("Length of 'depth' is not one")

  if(!is.character(return))
    stop("Data type of 'return' is not character")
  if(length(return) != 1)
    stop("Length of 'return' is not one")

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
  else if(return == "geneNumPredict")
    return(geneNumPredict)
  else if(return == "transcriptNum")
    return(geneNumPredict)
  else
    return(dropoutNum)
}




