#' estDropoutNum: Estimate dropout gene number in a cell
#'
#' This function is used to estimate dropout gene number in a cell for single-cell RNA-seq (scRNA-seq) data. It takes a non-negative vector of scRNA-seq raw read counts of a cell as input.
#'
#' @param sample A cell's raw read counts for each gene, could be a vector or a SingleCellExperiment object.
#' @param depth Relative sequencing depth to be predicted compared with initial sample depth, should between 0-100, default is 20.
#' @param histCounts Optional. Only needed when \code{sample} is blank or \code{sample = NULL}. A histogram table of raw read counts for the cell.
#' @param return A character for choosing the return value type of the function. "dropoutNum" (default) for dropout gene number, "geneNumPredict" for all expressed gene number predicted, "transcriptNum" for all transcript number predicted.
#' @return
#' The dropout gene number (or all expressed gene number) predicted in a cell.
#'
#' @author Zhun Miao.
#' @seealso
#' \code{\link{scRecover}}, for imputation of single-cell RNA-seq data.
#'
#' \code{\link{countsSampling}}, for downsampling the read counts in a cell.
#'
#' \code{\link{normalization}}, for normalization of single-cell RNA-seq data.
#'
#' \code{\link{scRecoverTest}}, a test dataset for scRecover.
#'
#' @examples
#' # Load test data
#' data(scRecoverTest)
#'
#' # Estimate dropout gene number in a cell
#' estDropoutNum(sample = counts[,1], return = "dropoutNum")
#'
#' # Estimate all expressed gene number in a cell
#' estDropoutNum(sample = counts[,1], return = "geneNumPredict")
#'
#'
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import penalized
#' @importFrom methods is
#' @importFrom kernlab specc
#' @importFrom rsvd rpca
#' @importFrom graphics hist
#' @importFrom stats complete.cases dgamma dnbinom dnorm median plogis prcomp quantile rgamma rnorm sd uniroot
#' @importFrom utils read.csv read.table write.csv write.table
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb fitdistr
#' @importFrom pscl zeroinfl
#' @importFrom bbmle mle2
#' @importFrom gamlss gamlssML
#' @importFrom preseqR ztnb.rSAC
#' @importFrom SAVER saver
#' @importFrom BiocParallel bpparam bplapply
#' @export



# Estimate dropout gene number in a cell
estDropoutNum <- function(sample = NULL, depth = 20, histCounts = NULL, return = "dropoutNum"){
  # Handle SingleCellExperiment
  if(is(sample, "SingleCellExperiment")){
    if(!requireNamespace("SingleCellExperiment"))
      stop("To use SingleCellExperiment as input, you should install the package firstly")
    if(dim(sample)[2] != 1)
      stop("'sample' should only contain one cell")
    else
      sample <- counts(sample)
  }
  
  # Invalid input control
  if(is.null(sample) & is.null(histCounts))
    stop("One of 'sample' and 'histCounts' must be specified")
  if(!is.null(sample) & !is.null(histCounts))
    stop("Only one of 'sample' and 'histCounts' should be specified")

  if(!is.vector(sample) & !is.data.frame(sample) & !is.matrix(sample) & class(sample)[1] != "dgCMatrix" & !is.integer(sample) & !is.numeric(sample) & !is.double(sample))
    stop("Wrong data type of 'sample'")
  if(sum(is.na(sample)) > 0)
    stop("NA detected in 'sample'")
  if(sum(sample < 0) > 0)
    stop("Negative value detected in 'sample'")
  if(all(sample == 0))
    stop("All elements of 'sample' are zero")

  if(!is.numeric(depth) & !is.integer(depth))
    stop("Data type of 'depth' is not numeric or integer")

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
  suppressWarnings(preseqR <- ztnb.rSAC(histCounts))
  geneNumPredict <- round(preseqR(depth))
  dropoutNum <- geneNumPredict - sum(histCounts[,2])

  if(return %in% c("geneNumPredict", "transcriptNum")){
    return(geneNumPredict)
  }else{
    return(dropoutNum)
  }

}




