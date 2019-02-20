# scRecover

### *Zhun Miao*
### *2018-05-30*


## Introduction

Imputation for single-cell RNA-seq data


## Citation

If you use **`scRecover`** in published research, please cite:

> 


## Installation

To install the *developmental version* from [**GitHub**](https://github.com/miaozhun/scRecover/):

```{r Installation from GitHub, eval = FALSE}
devtools::install_github("miaozhun/scRecover", build_vignettes = TRUE)
```

To load the installed **`scRecover`** in R:

```{r Load scRecover, eval = FALSE}
library(scRecover)
```


## Input

**`scRecover`** takes two inputs: `counts` and one of `Kcluster` or `labels`.

The input `counts` is a scRNA-seq **read counts matrix** or a **`SingleCellExperiment`** object which contains the read counts matrix. The rows of the matrix are genes and columns are cells.

`Kcluster` is an integer specifying the number of cell subpopulations. This parameter can be determined based on prior knowledge or clustering of raw data. `Kcluster` is used to determine the candidate neighbors of each cell.

`labels` is a character/integer vector specifying the cell type of each column in the raw count matrix. Only needed when `Kcluster = NULL`. Each cell type should have at least two cells for imputation.


## Test data

Users can load the test data in **`scRecover`** by

```{r Load scRecoverTest}
library(scRecover)
data(scRecoverTest)
```

The test data `counts` in `scRecoverTest` is a scRNA-seq read counts matrix which has 200 genes (rows) and 150 cells (columns).

```{r counts}
dim(counts)
counts[1:6, 1:6]
```

The object `labels` in `scRecoverTest` is a vector of integer specifying the cell types in the read counts matrix, corresponding to the columns of `counts`.

```{r labels}
length(labels)
table(labels)
```


## Usage

### With read counts matrix input

Here is an example to run **`scRecover`** with read counts matrix input:

```{r demo1, eval = FALSE}
# Load test data for scRecover
data(scRecoverTest)

# Run scRecover with Kcluster specified
scRecover(counts = counts, Kcluster = 2, parallel = TRUE)

# Run scRecover with labels specified
scRecover(counts = counts, labels = labels, parallel = TRUE)
```

### With SingleCellExperiment input

The [`SingleCellExperiment`](http://bioconductor.org/packages/SingleCellExperiment/) class is a widely used S4 class for storing single-cell genomics data. **`scRecover`** also could take the `SingleCellExperiment` data representation as input.

Here is an example to run **`scRecover`** with `SingleCellExperiment` input:

```{r demo2, eval = FALSE}
# Load library and the test data for scRecover
library(scRecover)
library(SingleCellExperiment)
data(scRecoverTest)

# Convert the test data in scRecover to SingleCellExperiment data representation
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)))

# Run scRecover with SingleCellExperiment input sce (Kcluster specified)
scRecover(counts = sce, Kcluster = 2, parallel = TRUE)

# Run scRecover with SingleCellExperiment input sce (labels specified)
scRecover(counts = sce, labels = labels, parallel = TRUE)
```


## Output

Imputed counts matrices will be saved in a folder named with prefix 'outputFile_scRecover_' under the current working directory.


## Parallelization

**`scRecover`** integrates parallel computing function with [`BiocParallel`](http://bioconductor.org/packages/BiocParallel/) package. Users could just set `parallel = TRUE` in function `scRecover` to enable parallelization and leave the `BPPARAM` parameter alone.

```{r demo3, eval = FALSE}
# Load library
library(scRecover)

# Run scRecover with Kcluster specified
scRecover(counts = counts, Kcluster = 2, parallel = TRUE)

# Run scRecover with labels specified
scRecover(counts = counts, labels = labels, parallel = TRUE)
```

Advanced users could use a `BiocParallelParam` object from package `BiocParallel` to fill in the `BPPARAM` parameter to specify the parallel back-end to be used and its configuration parameters.

### For Unix and Mac users

The best choice for Unix and Mac users is to use `MulticoreParam` to configure a multicore parallel back-end:

```{r demo4, eval = FALSE}
# Load library
library(scRecover)
library(BiocParallel)

# Set the parameters and register the back-end to be used
param <- MulticoreParam(workers = 18, progressbar = TRUE)
register(param)

# Run scRecover with 18 cores (Kcluster specified)
scRecover(counts = counts, Kcluster = 2, parallel = TRUE, BPPARAM = param)

# Run scRecover with 18 cores (labels specified)
scRecover(counts = counts, labels = labels, parallel = TRUE, BPPARAM = param)
```

### For Windows users

For Windows users, use `SnowParam` to configure a Snow back-end is a good choice:

```{r demo5, eval = FALSE}
# Load library
library(scRecover)
library(BiocParallel)

# Set the parameters and register the back-end to be used
param <- SnowParam(workers = 8, type = "SOCK", progressbar = TRUE)
register(param)

# Run scRecover with 8 cores (Kcluster specified)
scRecover(counts = counts, Kcluster = 2, parallel = TRUE, BPPARAM = param)

# Run scRecover with 8 cores (labels specified)
scRecover(counts = counts, labels = labels, parallel = TRUE, BPPARAM = param)
```

See the [*Reference Manual*](https://bioconductor.org/packages/release/bioc/manuals/BiocParallel/man/BiocParallel.pdf) of [`BiocParallel`](http://bioconductor.org/packages/BiocParallel/) package for more details of the `BiocParallelParam` class.


## Visualization of results



## Interpretation of results



## Help

Use `browseVignettes("scRecover")` to see the vignettes of **`scRecover`** in R after installation.

Use the following code in R to get access to the help documentation for **`scRecover`**:

```{r help1, eval = FALSE}
# Documentation for scRecover
?scRecover
```

```{r help2, eval = FALSE}
# Documentation for test data
?scRecoverTest
?counts
?labels
```

You are also welcome to contact the author by email for help.


## Author

*Zhun Miao* <<miaoz13@mails.tsinghua.edu.cn>>

MOE Key Laboratory of Bioinformatics; Bioinformatics Division and Center for Synthetic & Systems Biology, TNLIST; Department of Automation, Tsinghua University, Beijing 100084, China.

