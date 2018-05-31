#' ImputeSingle: Imputation for single-cell RNA-seq data
#'
#' This function is used to impute missing values in single-cell RNA-seq (scRNA-seq) data. It takes a non-negative matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object as input. So users should map the reads (obtained from sequencing libraries of the samples) to the corresponding genome and count the reads mapped to each gene according to the gene annotation to get the raw read counts matrix in advance.
#'
#' @param counts A non-negative integer matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object which contains the read counts matrix. The rows of the matrix are genes and columns are samples/cells.
#' @param Kcluster An integer specifying the number of cell subpopulations. This parameter can be determined based on prior knowledge or clustering of raw data. \code{Kcluster} is used to determine the candidate neighbors of each cell.
#' @param labels A character/integer vector specifying the cell type of each column in the raw count matrix. Only needed when \code{Kcluster = NULL}. Each cell type should have at least two cells for imputation.
#' @param parallel If FALSE (default), no parallel computation is used; if TRUE, parallel computation using \code{BiocParallel}, with argument \code{BPPARAM}.
#' @param BPPARAM An optional parameter object passed internally to \code{\link{bplapply}} when \code{parallel=TRUE}. If not specified, \code{\link{bpparam}()} (default) will be used.
#' @return
#' Imputed counts matrices will be saved in a folder named with prefix 'outputFile_ImputeSingle_' under the current working directory.
#'
#' @author Zhun Miao.
#' @seealso
#' \code{\link{TestImputeSingleData}}, a test dataset for ImputeSingle.
#'
#' @examples
#' # Load test data for ImputeSingle
#' data(TestImputeSingleData)
#'
#' # Run ImputeSingle with Kcluster specified
#' ImputeSingle(counts = counts, Kcluster = 2, parallel = TRUE)
#'
#' # Run ImputeSingle with labels specified
#' ImputeSingle(counts = counts, labels = labels, parallel = TRUE)
#'
#' @import stats
#' @import BiocParallel
#' @import SingleCellExperiment
#' @importFrom utils read.csv write.csv
#' @importFrom MASS glm.nb fitdistr
#' @importFrom pscl zeroinfl
#' @importFrom bbmle mle2
#' @importFrom gamlss gamlssML
#' @importFrom preseqR ztnb.mincount
#' @importFrom scImpute scimpute
#' @importFrom SAVER saver
#' @importFrom Rmagic run_magic
#' @export



ImputeSingle <- function(counts, Kcluster = NULL, labels = NULL, parallel = FALSE, BPPARAM = bpparam()){

  # Handle for SingleCellExperiment
  if(class(counts)[1] == "SingleCellExperiment")
    counts <- counts(counts)

  # Invalid input control
  if(!is.matrix(counts) & !is.data.frame(counts))
    stop("Wrong data type of 'counts'")
  if(sum(is.na(counts)) > 0)
    stop("NA detected in 'counts'")
  if(sum(counts < 0) > 0)
    stop("Negative value detected in 'counts'")
  if(all(counts == 0))
    stop("All elements of 'counts' are zero")
  if(any(colSums(counts) == 0))
    warning("Library size of zero detected in 'counts'")

  if(is.null(Kcluster) & is.null(labels))
    stop("One of 'Kcluster' and 'labels' must be specified")
  if(!is.null(Kcluster) & !is.null(labels))
    stop("Only one of 'Kcluster' and 'labels' should be specified")

  if(!is.null(Kcluster)){
    if(!is.numeric(Kcluster) & !is.integer(Kcluster))
      stop("'Kcluster' must be a numeric or an integer")
    if(length(Kcluster) != 1)
      stop("Length of 'Kcluster' is not one")
  }

  if(!is.null(labels)){
    if(!is.numeric(labels) & !is.character(labels))
      stop("'labels' must be numeric or character")
    if(ncol(counts) != length(labels))
      stop("Length of 'labels' must equal to column number of 'counts'")
    if(min(table(labels)) < 2)
      stop("Too few cells (< 2) in a cluster of 'labels'")
  }

  if(!is.logical(parallel))
    stop("Data type of 'parallel' is not logical")
  if(length(parallel) != 1)
    stop("Length of 'parallel' is not one")

  # File path
  outputFile_path <- paste0("./outputFile_ImputeSingle_", gsub(" ", "-", gsub(":", "-", Sys.time())), "/")
  tempFile_path <- paste0(outputFile_path, "tempFile/")
  dir.create(outputFile_path, showWarnings = FALSE)
  dir.create(tempFile_path, showWarnings = FALSE)
  write.csv(counts, file = paste0(outputFile_path, "raw_data.csv"))
  count_path <- paste0(outputFile_path, "raw_data.csv")

  # Run scImpute
  scimpute(# full path to raw count matrix
    count_path = count_path,
    infile = "csv",              # format of input file
    outfile = "csv",             # format of output file
    out_dir = tempFile_path,     # full path to output directory
    labeled = is.null(Kcluster), # cell type labels not available
    drop_thre = 0.5,             # threshold set on dropout probability
    Kcluster = Kcluster,         # 2 cell subpopulations
    labels = labels,             # Each cell type should have at least two cells for imputation
    ncores = 1)                  # number of cores used in parallel computation
  counts_scImpute <- read.csv(file = paste0(tempFile_path, "scimpute_count.csv"), header = TRUE, row.names = 1)

  # # Run SAVER and MAGIC
  # library(doParallel)
  # cl <- makeCluster(4, outfile = "")
  # registerDoParallel(cl)
  # counts_SAVER <- saver(counts_raw)
  # stopCluster(cl)
  # counts_SAVER <- counts_SAVER$estimate
  # counts_MAGIC <- run_magic(counts_raw, t = 6)
  # write.csv(counts_SAVER, file = paste0(tempFile_path, "SAVER_count.csv"))
  # write.csv(counts_MAGIC, file = paste0(tempFile_path, "MAGIC_count.csv"))

  # Get cluster information
  if(is.null(labels))
    clust <- readRDS(paste0(tempFile_path, "clust.rds"))
  if(!is.null(labels))
    clust <- labels
  names(clust) <- colnames(counts)
  nclust <- sum(!is.na(unique(clust)))

  # Eliminate outliers and normalization
  counts_used <- counts[,!is.na(clust)]
  counts_norm <- normalization(counts_used)
  whether_impute <- counts
  whether_impute[,] <- FALSE

  # # Estimate dropout gene number in each cell
  dropoutNum <- NULL
  if(!parallel){
    for(i in 1:ncol(counts_used)){
      cat("\r",paste0("ImputeSingle is estimating dropout gene number in ", i, " of ", ncol(counts_used), " non-outlier cells"))
      dropoutNum <- c(dropoutNum, estDropoutNum(counts_used[,i]))
    }
    names(dropoutNum) <- colnames(counts_used)
  }else{
    message("ImputeSingle is estimating dropout gene number in ", ncol(counts_used), " non-outlier cells")
    dropoutNum <- do.call(c, bplapply(as.data.frame(counts_used), FUN = estDropoutNum, BPPARAM = BPPARAM))
  }

  # ZINB MLE of each cluster gene by gene
  ZINB_parameters_list <- list()
  P_dropout_cc_list <- list()
  P_dropout_mat <- NULL
  for(cc in 1:nclust){
    message("Processing ", cc, " of ", nclust, " cell clusters")
    cells <- names(clust[clust %in% unique(clust[!is.na(clust)])[cc]])
    counts_norm_cc <- counts_norm[, cells, drop = FALSE]

    # ZINB MLE of each cluster gene by gene
    ZINB_parameters <- NULL
    if(!parallel){
      for(i in 1:nrow(counts_norm_cc)){
        cat("\r",paste0("ImputeSingle is analyzing ", i, " of ", nrow(counts_norm), " genes in cluster ", cc))
        ZINB_parameters <- rbind(ZINB_parameters, mleZINB(counts_norm_cc[i,]))
      }
      row.names(ZINB_parameters) <- row.names(counts_norm_cc)
      message("\r")
    }else{
      message("ImputeSingle is analyzing ", nrow(counts_norm), " genes in cluster ", cc)
      ZINB_parameters <- do.call(rbind, bplapply(as.data.frame(t(counts_norm_cc)), FUN = mleZINB, BPPARAM = BPPARAM))
    }
    ZINB_parameters_list[[cc]] <- ZINB_parameters

    # Estimate dropout probability of each gene
    P_dropout <- NULL
    for(i in 1:nrow(counts_norm_cc)){
      if(any(is.na(ZINB_parameters[i,])))
        P_dropout <- c(P_dropout, 0)
      else
        P_dropout <- c(P_dropout, (1 - ZINB_parameters[i, "theta"]) * dnbinom(x=0, size=ZINB_parameters[i, "size"], prob=ZINB_parameters[i, "prob"]) / (ZINB_parameters[i, "theta"] + (1 - ZINB_parameters[i, "theta"]) * dnbinom(x=0, size=ZINB_parameters[i, "size"], prob=ZINB_parameters[i, "prob"])))
    }
    names(P_dropout) <- row.names(counts_norm_cc)
    P_dropout_mat <- cbind(P_dropout_mat, P_dropout)

    # Get dropout probability for each gene in each cell
    P_dropout_cc <- counts_norm_cc
    P_dropout_cc[,1:ncol(counts_norm_cc)] <- P_dropout
    P_dropout_cc[counts_norm_cc != 0] <- 0
    P_dropout_cc_list[[cc]] <- P_dropout_cc

    # Determine the position need to be imputed
    P_dropout_rank <- apply(-P_dropout_cc, 2, rank)
    whether_impute_cc <- sweep(P_dropout_rank, MARGIN = 2, dropoutNum[cells], FUN = "<=")
    whether_impute[row.names(whether_impute_cc), cells] <- whether_impute_cc
  }
  colnames(P_dropout_mat) <- paste0("CellCluster_", 1:nclust)

  # Get whether_impute matrix for counts
  whether_impute_iz <- whether_impute
  whether_impute_inz <- whether_impute
  whether_impute_inz[counts != 0] <- TRUE

  # Imputation with whether_impute to results of scImpute, SAVER and MAGIC
  counts_scImpute_iz <- counts + counts_scImpute * whether_impute_iz
  counts_scImpute_inz <- counts_scImpute * whether_impute_inz

  # counts_SAVER_iz <- counts + counts_SAVER * whether_impute_iz
  # counts_SAVER_inz <- counts_SAVER * whether_impute_inz
  #
  # counts_MAGIC_iz <- counts + counts_MAGIC * whether_impute_iz
  # counts_MAGIC_inz <- counts_MAGIC * whether_impute_inz

  # Output files
  write.csv(counts_scImpute_iz, file = paste0(outputFile_path, "counts_scImpute_iz.csv"))
  write.csv(counts_scImpute_inz, file = paste0(outputFile_path, "counts_scImpute_inz.csv"))
  # write.csv(counts_SAVER_iz, file = paste0(outputFile_path, "counts_SAVER_iz.csv"))
  # write.csv(counts_SAVER_inz, file = paste0(outputFile_path, "counts_SAVER_inz.csv"))
  # write.csv(counts_MAGIC_iz, file = paste0(outputFile_path, "counts_MAGIC_iz.csv"))
  # write.csv(counts_MAGIC_inz, file = paste0(outputFile_path, "counts_MAGIC_inz.csv"))


}




