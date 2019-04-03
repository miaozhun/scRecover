#' scRecover: Imputation for single-cell RNA-seq data
#'
#' This function is used to impute missing values in single-cell RNA-seq (scRNA-seq) data. It takes a non-negative matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object as input. So users should map the reads (obtained from sequencing libraries of the samples) to the corresponding genome and count the reads mapped to each gene according to the gene annotation to get the raw read counts matrix in advance.
#'
#' @param counts A non-negative integer matrix of scRNA-seq raw read counts or a \code{SingleCellExperiment} object which contains the read counts matrix. The rows of the matrix are genes and columns are samples/cells.
#' @param Kcluster An integer specifying the number of cell subpopulations. This parameter can be determined based on prior knowledge or clustering of raw data. \code{Kcluster} is used to determine the candidate neighbors of each cell.
#' @param labels Optional. Only needed when \code{Kcluster} is blank or \code{Kcluster = NULL}. A character/integer vector specifying the cell type of each column in the raw count matrix. Each cell type should have at least two cells.
#' @param outputDir The path of the output directory. If not specified, a folder named with prefix 'outDir_scRecover_' under the current working directory will be used.
#' @param depth Relative sequencing depth to be predicted compared with initial sample depth, should between 2-100, default is 20.
#' @param SAVER Whether use and improve SAVER in imputation, default is FALSE.
#' @param MAGIC Whether use and improve MAGIC in imputation, default is FALSE.
#' @param UMI Whether use full UMI data, default is FALSE. If TRUE, \code{hist_raw_counts} and \code{hist_RUG_counts} should be specified.
#' @param hist_raw_counts A list contains the histogram table of raw read counts for each cell in \code{counts}.
#' @param hist_RUG_counts A list contains the histogram table of raw UMI-gene counts for each cell in \code{counts}.
#' @param parallel If FALSE (default), no parallel computation is used; if TRUE, parallel computation using \code{BiocParallel}, with argument \code{BPPARAM}.
#' @param BPPARAM An optional parameter object passed internally to \code{\link{bplapply}} when \code{parallel=TRUE}. If not specified, \code{\link{bpparam}()} (default) will be used.
#' @param verbose Whether to show specific calculation progress, default is TRUE.
#' @return
#' Imputed counts matrices will be saved in the output directory specified by \code{outputDir}.
#'
#' @author Zhun Miao.
#' @seealso
#' \code{\link{estDropoutNum}}, for estimating dropout gene number in a cell.
#'
#' \code{\link{countsSampling}}, for downsampling the read counts in a cell.
#'
#' \code{\link{normalization}}, for normalization of single-cell RNA-seq data.
#'
#' \code{\link{scRecoverTest}}, a test dataset for scRecover.
#'
#' @examples
#' # Load test data for scRecover
#' data(scRecoverTest)
#'
#' # Run scRecover with Kcluster specified
#' scRecover(counts = counts, Kcluster = 2)
#'
#' # Or run scRecover with labels specified
#' # scRecover(counts = counts, labels = labels)
#'
#'
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import penalized
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



scRecover <- function(counts, Kcluster = NULL, labels = NULL, outputDir = NULL, depth = 20, SAVER = FALSE, MAGIC = FALSE, UMI = FALSE, hist_raw_counts = NULL, hist_RUG_counts = NULL, parallel = FALSE, BPPARAM = bpparam(), verbose = TRUE){

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
    if(!is.numeric(labels) & !is.integer(labels) & !is.character(labels))
      stop("'labels' must be numeric or integer or character")
    if(ncol(counts) != length(labels))
      stop("Length of 'labels' must equal to column number of 'counts'")
    if(min(table(labels)) < 2)
      stop("Too few cells (< 2) in a cluster of 'labels'")
  }

  if(!is.numeric(depth) & !is.integer(depth))
    stop("Data type of 'depth' is not numeric or integer")
  if(length(depth) != 1)
    stop("Length of 'depth' is not one")

  if(!is.logical(SAVER))
    stop("Data type of 'SAVER' is not logical")
  if(length(SAVER) != 1)
    stop("Length of 'SAVER' is not one")

  if(!is.logical(MAGIC))
    stop("Data type of 'MAGIC' is not logical")
  if(length(MAGIC) != 1)
    stop("Length of 'MAGIC' is not one")

  if(!is.logical(UMI))
    stop("Data type of 'UMI' is not logical")
  if(length(UMI) != 1)
    stop("Length of 'UMI' is not one")

  if(UMI){
    if(is.null(hist_raw_counts) | is.null(hist_RUG_counts))
      stop("'hist_raw_counts' and 'hist_RUG_counts' must be specified if UMI == TRUE")
    if(!is.list(hist_raw_counts) | !is.list(hist_RUG_counts))
      stop("'hist_raw_counts' and 'hist_RUG_counts' must be lists")
    if(length(hist_raw_counts) != ncol(counts) | length(hist_RUG_counts) != ncol(counts))
      stop("'hist_raw_counts' and 'hist_RUG_counts' must have length equal to column number of 'counts'")
  }

  if(!is.logical(parallel))
    stop("Data type of 'parallel' is not logical")
  if(length(parallel) != 1)
    stop("Length of 'parallel' is not one")

  # Preprocessing
  if(is.null(outputDir))
    outputDir <- paste0(tempdir(), "/outDir_scRecover_", gsub(" ", "-", gsub(":", "-", Sys.time())), "/")
  tempFileDir <- paste0(outputDir, "tempFile/")
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tempFileDir, showWarnings = FALSE, recursive = TRUE)
  counts <- as.matrix(counts)
  write.csv(counts, file = paste0(outputDir, "raw_data.csv"))
  count_path <- paste0(outputDir, "raw_data.csv")

  # Run SAVER
  if(SAVER){
    print("========================== Running SAVER ==========================")
    counts_SAVER <- saver(counts, ncores = if(!parallel) 1 else detectCores() - 2)$estimate
    write.csv(counts_SAVER, file = paste0(tempFileDir, "SAVER_count.csv"))
    print("========================== SAVER finished =========================")
  }

  # Run MAGIC
  if(MAGIC){
    print("========================== Running MAGIC ==========================")
    counts_MAGIC <- t(magic(t(counts), n.jobs = if(!parallel) 1 else -3)[[1]])
    write.csv(counts_MAGIC, file = paste0(tempFileDir, "MAGIC_count.csv"))
    print("========================== MAGIC finished =========================")
  }

  # Run scImpute
  print("========================= Running scImpute =========================")
  scimpute(
    count_path = count_path,      # full path to raw count matrix
    infile = "csv",               # format of input file
    outfile = "csv",              # format of output file
    out_dir = tempFileDir,        # full path to output directory
    labeled = is.null(Kcluster),  # cell type labels not available
    drop_thre = 0.5,              # threshold set on dropout probability
    Kcluster = Kcluster,          # 2 cell subpopulations
    labels = labels,              # Each cell type should have at least two cells for imputation
    ncores = if(parallel & .Platform$OS.type != "windows") detectCores() - 2 else 1)    # number of cores used
  # if(parallel & .Platform$OS.type != "windows"){
  #   env <- foreach:::.foreachGlobals
  #   rm(list=ls(name=env), pos=env)
  # }
  counts_scImpute <- read.csv(file = paste0(tempFileDir, "scimpute_count.csv"), header = TRUE, row.names = 1)
  print("========================= scImpute finished ========================")

  # Get cluster information
  if(is.null(labels))
    clust <- readRDS(paste0(tempFileDir, "clust.rds"))
  if(!is.null(labels))
    clust <- labels
  names(clust) <- colnames(counts)
  nclust <- sum(!is.na(unique(clust)))

  # Eliminate outliers and normalization
  counts_used <- counts[,!is.na(clust)]
  counts_norm <- normalization(counts_used)
  hist_raw_counts <- hist_raw_counts[colnames(counts_used)]
  hist_RUG_counts <- hist_RUG_counts[colnames(counts_used)]
  whether_impute <- counts
  whether_impute[,] <- FALSE

  # Estimate dropout number in each cell
  if(!UMI){
    dropoutNum <- rep(NA, ncol(counts_used))
    if(!parallel){
      for(i in 1:ncol(counts_used)){
        if(verbose){
          cat("\r",paste0("scRecover is estimating dropout gene number in ", i, " of ", ncol(counts_used), " non-outlier cells"))
        }
        dropoutNum[i] <- estDropoutNum(sample = counts_used[,i], depth = depth, histCounts = NULL, return = "dropoutNum")
      }
      names(dropoutNum) <- colnames(counts_used)
      message("\r")
    }else{
      message("scRecover is estimating dropout gene number in ", ncol(counts_used), " non-outlier cells")
      dropoutNum <- do.call(c, bplapply(as.data.frame(counts_used), depth = depth, histCounts = NULL, return = "dropoutNum", FUN = estDropoutNum, BPPARAM = BPPARAM))
    }
  }

  # Estimate dropout number and transcript number in each cell (UMI)
  if(UMI){
    dropoutNum <- rep(NA, ncol(counts_used))
    transcriptNum <- rep(NA, ncol(counts_used))
    if(!parallel){
      for(i in 1:ncol(counts_used)){
        if(verbose){
          cat("\r",paste0("scRecover is estimating dropout gene number in ", i, " of ", ncol(counts_used), " non-outlier cells (UMI)"))
        }
        dropoutNum[i] <- estDropoutNum(sample = NULL, depth = depth, histCounts = hist_raw_counts[[i]], return = "dropoutNum")
      }
      names(dropoutNum) <- colnames(counts_used)
      message("\r")

      for(i in 1:ncol(counts_used)){
        if(verbose){
          cat("\r",paste0("scRecover is estimating transcript number in ", i, " of ", ncol(counts_used), " non-outlier cells (UMI)"))
        }
        transcriptNum[i] <- estDropoutNum(sample = NULL, depth = depth, histCounts = hist_RUG_counts[[i]], return = "transcriptNum")
      }
      names(transcriptNum) <- colnames(counts_used)
      message("\r")
    }else{
      message("scRecover is estimating dropout gene number in ", ncol(counts_used), " non-outlier cells (UMI)")
      dropoutNum <- do.call(c, bplapply(hist_raw_counts, sample = NULL, depth = depth, return = "dropoutNum", FUN = estDropoutNum, BPPARAM = BPPARAM))

      message("scRecover is estimating transcript number in ", ncol(counts_used), " non-outlier cells (UMI)")
      transcriptNum <- do.call(c, bplapply(hist_RUG_counts, sample = NULL, depth = depth, return = "transcriptNum", FUN = estDropoutNum, BPPARAM = BPPARAM))
    }
    saveRDS(transcriptNum, file = paste0(tempFileDir, "transcriptNum.rds"))
  }

  # Save dropoutNum
  saveRDS(dropoutNum, file = paste0(tempFileDir, "dropoutNum.rds"))

  # Processing data by cell clusters
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
        if(verbose){
          cat("\r",paste0("scRecover is analyzing ", i, " of ", nrow(counts_norm), " genes in cluster ", cc))
        }
        ZINB_parameters <- rbind(ZINB_parameters, mleZINB(counts_norm_cc[i,]))
      }
      row.names(ZINB_parameters) <- row.names(counts_norm_cc)
      message("\r")
    }else{
      message("scRecover is analyzing ", nrow(counts_norm), " genes in cluster ", cc)
      ZINB_parameters <- do.call(rbind, bplapply(as.data.frame(t(counts_norm_cc)), FUN = mleZINB, BPPARAM = BPPARAM))
    }
    ZINB_parameters_list[[cc]] <- ZINB_parameters

    # Estimate dropout probability of each gene
    P_dropout <- rep(NA, nrow(counts_norm_cc))
    for(i in 1:nrow(counts_norm_cc)){
      if(any(is.na(ZINB_parameters[i,])))
        P_dropout[i] <- 0
      else
        P_dropout[i] <- (1 - ZINB_parameters[i, "theta"]) * dnbinom(x=0, size=ZINB_parameters[i, "size"], prob=ZINB_parameters[i, "prob"]) / (ZINB_parameters[i, "theta"] + (1 - ZINB_parameters[i, "theta"]) * dnbinom(x=0, size=ZINB_parameters[i, "size"], prob=ZINB_parameters[i, "prob"]))
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
  whether_impute_iz[counts != 0] <- FALSE
  whether_impute_inz <- whether_impute
  whether_impute_inz[counts != 0] <- TRUE
  saveRDS(whether_impute_iz, file = paste0(tempFileDir, "whether_impute_iz.rds"))
  saveRDS(whether_impute_inz, file = paste0(tempFileDir, "whether_impute_inz.rds"))

  # Imputation with whether_impute to results of scImpute, SAVER and MAGIC
  counts_scImpute_inz <- as.matrix(counts_scImpute * whether_impute_inz)
  if(SAVER)
    counts_SAVER_inz <- as.matrix(counts_SAVER * whether_impute_inz)
  if(MAGIC)
    counts_MAGIC_inz <- as.matrix(counts_MAGIC * whether_impute_inz)

  if(UMI){
    transcriptNum_all <- colSums(counts)
    transcriptNum_all[names(transcriptNum)] <- transcriptNum
    counts_scImpute_inz <- counts_scImpute_inz %*% diag(transcriptNum_all/colSums(counts_scImpute_inz))
    if(SAVER)
      counts_SAVER_inz <- counts_SAVER_inz %*% diag(transcriptNum_all/colSums(counts_SAVER_inz))
    if(MAGIC)
      counts_MAGIC_inz <- counts_MAGIC_inz %*% diag(transcriptNum_all/colSums(counts_MAGIC_inz))
  }

  # Output files
  write.csv(counts_scImpute_inz, file = paste0(outputDir, "scRecover+scImpute.csv"))
  if(SAVER)
    write.csv(counts_SAVER_inz, file = paste0(outputDir, "scRecover+SAVER.csv"))
  if(MAGIC)
    write.csv(counts_MAGIC_inz, file = paste0(outputDir, "scRecover+MAGIC.csv"))
  print("======================== scRecover finished ========================")
  print(paste0("The output files of scRecover are in ", outputDir))
  print("========================= Congratulations! =========================")

}




