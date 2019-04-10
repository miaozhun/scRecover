#' countsSampling: Downsampling the read counts in a cell
#'
#' This function is used to downsample the read counts in a cell for single-cell RNA-seq (scRNA-seq) data. It takes a non-negative vector of scRNA-seq raw read counts of a cell as input.
#'
#' @param counts A cell's raw read counts for each gene, could be a vector or a SingleCellExperiment object.
#' @param fraction Fraction of reads to be downsampled, should be between 0-1, default is 0.1.
#' @return
#' A vector of the downsampled read counts of each gene in the cell.
#'
#' @author Zhun Miao.
#' @seealso
#' \code{\link{scRecover}}, for imputation of single-cell RNA-seq data.
#'
#' \code{\link{estDropoutNum}}, for estimating dropout gene number in a cell.
#'
#' \code{\link{normalization}}, for normalization of single-cell RNA-seq data.
#'
#' \code{\link{scRecoverTest}}, a test dataset for scRecover.
#'
#' @examples
#' # Load test data
#' data(scRecoverTest)
#'
#' # Downsample the read counts in oneCell
#' oneCell.down <- countsSampling(counts = oneCell, fraction = 0.1)
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
#' @importFrom Rmagic magic
#' @importFrom BiocParallel bpparam bplapply
#' @export



# Downsampling the read counts in a cell
countsSampling <- function(counts, fraction = 0.1){
  # Handle SingleCellExperiment
  if(is(counts, "SingleCellExperiment")){
    if(!requireNamespace("SingleCellExperiment"))
      stop("To use SingleCellExperiment as input, you should install the package firstly")
    if(dim(counts)[2] != 1)
      stop("'counts' should only contain one cell")
    else
      counts <- counts(counts)
  }
  
  # Invalid input control
  if(!is.vector(counts) & !is.data.frame(counts) & !is.matrix(counts) & class(counts)[1] != "dgCMatrix" & !is.integer(counts) & !is.numeric(counts) & !is.double(counts))
    stop("Wrong data type of 'counts'")
  if(sum(is.na(counts)) > 0)
    stop("NA detected in 'counts'")
  if(sum(counts < 0) > 0)
    stop("Negative value detected in 'counts'")
  if(all(counts == 0))
    stop("All elements of 'counts' are zero")

  if(!is.numeric(fraction) & !is.double(fraction))
    stop("Data type of 'fraction' is not numeric or double")
  if(length(fraction) != 1)
    stop("Length of 'fraction' is not one")
  if(fraction < 0 | fraction > 1)
    stop("'fraction' should be between 0 and 1")

  # Downsample
  n <- floor(sum(counts) * fraction)
  readsGet <- sort(sample(seq_len(sum(counts)), n))
  cumCounts <-  c(0, cumsum(counts))
  counts_New <- counts
  counts_New[] <- hist(readsGet, breaks = cumCounts, plot=FALSE)$count
  return(counts_New)
}




