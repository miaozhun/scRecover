#' normalization: Normalization for single-cell RNA-seq data
#'
#' This function is used to normalize single-cell RNA-seq (scRNA-seq) data. It takes a non-negative matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object as input.
#'
#' @param counts A non-negative integer matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object which contains the read counts matrix. The rows of the matrix are genes and columns are samples/cells.
#' @return
#' A normalized scRNA-seq read counts matrix.
#'
#' @author Zhun Miao.
#' @seealso
#' \code{\link{scRecover}}, for imputation of single-cell RNA-seq data.
#'
#' \code{\link{estDropoutNum}}, for estimating dropout gene number in a cell.
#'
#' \code{\link{countsSampling}}, for downsampling the read counts in a cell.
#'
#' \code{\link{scRecoverTest}}, a test dataset for scRecover.
#'
#' @examples
#' # Load test data
#' data(scRecoverTest)
#'
#' # Normalization of counts
#' counts.norm <- normalization(counts = counts)
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



# Normalization function
normalization <- function(counts){
  # Handle SingleCellExperiment
  if(is(counts, "SingleCellExperiment")){
    if(!requireNamespace("SingleCellExperiment"))
      stop("To use SingleCellExperiment as input, you should install the package firstly")
    counts <- counts(counts)
  }

  # Invalid input control
  if(!is.matrix(counts) & !is.data.frame(counts) & class(counts)[1] != "dgCMatrix")
    stop("Wrong data type of 'counts'")
  if(sum(is.na(counts)) > 0)
    stop("NA detected in 'counts'");gc();
  if(sum(counts < 0) > 0)
    stop("Negative value detected in 'counts'");gc();
  if(all(counts == 0))
    stop("All elements of 'counts' are zero");gc();
  if(any(colSums(counts) == 0))
    warning("Library size of zero detected in 'counts'");gc();

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
  for (i in seq_len(geneNum_NAZ))
  {
    gene_NZ <- counts_NAZ[i,counts_NAZ[i,] > 0]
    GEOmean[i] <- exp(sum(log(gene_NZ), na.rm=TRUE) / length(gene_NZ))
  }
  S <- rep(NA, sampleNum)
  counts_norm <- counts_NAZ
  for (j in seq_len(sampleNum))
  {
    sample_j <- counts_NAZ[,j]/GEOmean
    S[j] <- median(sample_j[which(sample_j != 0)])
    counts_norm[,j] <- counts_NAZ[,j]/S[j]
  }
  counts_norm <- ceiling(counts_norm)
  return(counts_norm)
}




