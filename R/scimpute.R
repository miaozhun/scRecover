# scImpute is import from https://github.com/Vivianstats/scImpute
# scImpute original authors: Wei Vivian Li & Jingyi Jessica Li
# scImpute reference: Li, Wei Vivian, and Jingyi Jessica Li. "An accurate and robust imputation method scImpute for single-cell RNA-seq data." Nature communications 9.1 (2018): 997.
scimpute <-
function (count_path, infile = "csv", outfile = "csv", type = "count", out_dir, labeled = FALSE,
          drop_thre = 0.5, Kcluster = NULL, labels = NULL, genelen = NULL, ncores = 5)
{
    if(labeled == TRUE & is.null(labels)){
      stop("'labels' must be specified when 'labeled = TRUE'!")
    }
    if(labeled == FALSE & is.null(Kcluster)){
      stop("'Kcluster' must be specified when 'labeled = FALSE'!")
    }
    if(!(type %in% c("count", "TPM"))){ stop("expression values can be either 'count' or 'TPM'!") }
    if(type == "TPM" & is.null(genelen)){ stop("'genelen' must be specified when type = 'TPM'!") }

    # print(drop_thre)
    print("reading in raw count matrix ...")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    count_lnorm = read_count(filetype = infile, path = count_path, out_dir = out_dir,
                             type = type, genelen = genelen)
    print("reading finished!")

    if(labeled){
      if(length(labels) != ncol(count_lnorm)){
        stop("number of cells does not match number of labels !")
      }
    }
    genenames = rownames(count_lnorm)
    cellnames = colnames(count_lnorm)

    print("imputation starts ...")
    if (!labeled){
      res_imp = imputation_model8(count = count_lnorm, labeled = FALSE,
                                  point = log10(1.01), drop_thre = drop_thre,
                                  Kcluster = Kcluster,
                                  out_dir = out_dir, ncores = ncores)
    }else{
      res_imp = imputation_wlabel_model8(count = count_lnorm, labeled = TRUE,
                                         cell_labels = labels, point = log10(1.01),
                                         drop_thre = drop_thre,
                                         Kcluster = NULL, out_dir = out_dir,
                                         ncores = ncores)
    }
    count_imp = res_imp$count_imp
    outliers = res_imp$outlier
    count_imp = 10^count_imp - 1.01
    rownames(count_imp) = genenames
    colnames(count_imp) = cellnames
    print("writing imputed count matrix ...")
    write_count(count_imp, filetype = outfile, out_dir, type = type, genelen = genelen)
    return(outliers)
}



