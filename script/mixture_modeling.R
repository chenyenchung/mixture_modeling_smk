#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(loo))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(txtplot))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(R.utils))
source("script/utils.R")

# Since we are running out model multiple times, recompiling Stan code is
# not efficient.
# Setting auto_write = TRUE asks rstan to compile the Stan model
# to C++ and save as RDS files
rstan_options(auto_write = TRUE, iter = 5e4)


#####################################################################
##                                                                 ##
##     Run main analysis if command line options were provided     ##
##                                                                 ##
#####################################################################

## Checking Stan options and use default values if not given
#
# ## nIter
# ## How many iterations to run to fit the model
# # Set default value if NULL
# nIter <- as.integer(snakemake@params[["nIter"]])
#
# ## k-fold cross-validation
# numK <- as.integer(snakemake@params[["numK"]])
#
# upm <- read.csv(
#   snakemake@input[['mtx']],
#     row.names = 1
#   )
#
# ### The expression matrix is expected to consist only of numbers
# ### Check if the expression matrix is numeric and show error message
# ### if it is not.
# num_col <- apply(upm, 2, is.numeric)
# if (!all(num_col)) {
#   message(
#     paste0(
#       "The following columns of the expression matrix are not numbers:\n",
#       colnames(upm)[!num_col]
#     )
#   )
#   stop("All columns of the expression matrix must be numbers.")
# }
#
# ## Load metadata csv
# metadata <- read.csv(snakemake@input[['meta']], row.names = 1)
#
# ## TODO: Compare assigning as individual driver or one driver
# ### Set pseudobulk if the column is not already present in the metadata table
# ### This corresponds to the driver used in the original code from Davis et al.
# ### Internally, driver identity represents another factor other than cell type
# ### in the mixture modeling.
# if (!"pseudobulk" %in% colnames(metadata)) {
#   metadata$pseudobulk <- 1
# }
#
#
# ## Ensure the metadata table is the same as column names
# ## as column names are modified at cryptic times even if check.names = FALSE
# metadata$sample <- make.names(metadata$sample)
#
# ### The metadata should describe in order each column of the expression matrix
# ### If not, show error message and stop.
# if (!identical(colnames(upm), metadata$sample)) {
#   stop(
#     "The 'sample' column must be the same as the column names of the expression matrix."
#   )
# }
#
# ### Make sure all columns in metadata are characters
# if (!is.character(metadata$sample)) {
#   metadata$sample <- as.character(metadata$sample)
# }
# if (!is.character(metadata$cell)) {
#   metadata$cell <- as.character(metadata$cell)
# }
# if (!is.character(metadata$pseudobulk)) {
#   metadata$pseudobulk <- as.character(metadata$pseudobulk)
# }

#### Define input object to fit original code
# inDat <- list(
#   nGenes = nrow(upm),
#   nCells = length(unique(metadata$cell)),
#   nSamples = length(unique(metadata$sample)),
#   nDrivers = length(unique(metadata$pseudobulk)),
#   O = upm,
#   logC = upm,
#   cell = metadata$cell,
#   driver = metadata$pseudobulk,
#   sample = colnames(upm)
# )
nIter <- 500
numK <- 10
t1 <- 500
t2 <- 500
upm <- matrix(
  c(
    rnbinom(t1, mu = 1, size = 3), rnbinom(t2, mu = 10, size = 3),
    rnbinom(t1, mu = 0.5, size = 3), rnbinom(t2, mu = 2.5, size = 3),
    rnbinom(t1, mu = 0.01, size = 3) * 100, rnbinom(t2, mu = 1, size = 3) * 100
  ), nrow = 3, byrow = TRUE
)
row.names(upm) <- paste0("gene", seq(nrow(upm)))
metadata <- data.frame(
  sample = seq(t1 + t2),
  cell = c(rep(1L, t1), rep(2L, t2)),
  pseudobulk = 1
)

colnames(upm) <- metadata$sample

inDat <- list(
  nGenes = 1,
  nCells = ncol(upm),
  nSamples = length(unique(metadata$sample)),
  nDrivers = length(unique(metadata$pseudobulk)),
  O = upm,
  logC = upm,
  cell = metadata$cell,
  driver = metadata$pseudobulk,
  sample = colnames(upm)
)

## Setup individual STAN temporary file path for each run
# temp_path <- paste(strsplit(snakemake@output[[1]], "/")[[1]][1:2], collapse = "/")

temp_path <- "script/STAN"
# dir.create(temp_path)

# stan_scripts <- list.files("script/STAN/", pattern = "\\.stan$")
# file.copy(
#   paste0("script/STAN/", stan_scripts),
#   temp_path
# )


# Define paths that point to the Stan scripts
level1.bim.stanFn <- paste0(temp_path, "/bimodal.stan")
level1.bim.cv.stanFn <- paste0(temp_path, "/bimodal_CV.stan")
level1.uni.stanFn <- paste0(temp_path, "/unimodal.stan")
level1.uni.cv.stanFn <- paste0(temp_path, "/unimodal_CV.stan")


## Extract values to be modeled from the input object
## This is more meaningful as Davis et al. used the object for other
## purposes.
## In my case, keeping this is just to avoid massively changing the data
## processing logic of the modeling code.
geneList <- row.names(upm)
genes <- geneList
nGenes <- length(geneList)
nSamples <- inDat$nSamples
samples <- inDat$sample
cells <- inDat$cell
nCells <- length(unique(inDat$cell))
drivers <- inDat$driver
nDrivers <- length(unique(inDat$driver))
logC <- inDat$logC


####################################################################
##                                                                ##
##                      Model fitting                             ##
##                                                                ##
####################################################################

####### Level 1: Fit Unimodal and Bimodal mixture models with Stan

### Run on one gene to compile Stan models
### This only useful when auto_write is set to TRUE.

## Bimodal
stanOptions <- list()
stanOptions$bim <- list(
  file = level1.bim.stanFn,
  data = list(
    nSamples = nSamples,
    nCells = nCells,
    nDrivers = nDrivers,
    logE = as.numeric(inDat$logC[genes[1],]),
    cell = as.integer(factor(inDat$cell, levels=unique(inDat$cell))),
    driver = as.integer(factor(inDat$driver, levels=unique(inDat$driver)))
  ),
  iter    = nIter,
  chains  = 4,
  cores   = parallel::detectCores(),
  verbose = FALSE,
  ## Defalt for RStan is 0.8
  control = list(adapt_delta = 0.98),
  ## Suppress printing progress to stdout if run as a script
  refresh = ifelse(
    interactive(),
    max(nIter/10, 1),
    -1
  )
)

stanOptions$bim$fit <- do.call(stan, stanOptions$bim)

## Unimodal for compilation
stanOptions$uni <- stanOptions$bim
stanOptions$uni$file <- level1.uni.stanFn
stanOptions$uni$fit <- NULL
stanOptions$uni$fit <- do.call(stan, stanOptions$uni)

## Partitioning for k-fold CV
kpart1 <- inferState.partitionK(
  n = nSamples,
  k = numK,
  labels = as.integer(as.factor(inDat$driver))
)

## k-fold CV bimodal for compilation
stanOptions$bim.cv <- stanOptions$bim
stanOptions$bim.cv$file <- level1.bim.cv.stanFn
stanOptions$bim.cv$data <- list(
  nSamples_t = length(which(kpart1 != 1)),
  nSamples_h = length(which(kpart1 == 1)),
  logE_t = as.numeric(inDat$logC[genes[1], which(kpart1 != 1)]),
  logE_h = as.numeric(inDat$logC[genes[1], which(kpart1 == 1)])
)
stanOptions$bim.cv$fit <- NULL
stanOptions$bim.cv$fit <- do.call(stan, stanOptions$bim.cv)

## k-fold CV unimodal for compilation
stanOptions$uni.cv <- stanOptions$bim.cv
stanOptions$uni.cv$file <- level1.uni.cv.stanFn
stanOptions$uni.cv$fit <- NULL
stanOptions$uni.cv$fit <- do.call(stan, stanOptions$uni.cv)

################# Full run ##########################
start <- Sys.time()
geneEfit <- lapply(
  genes,
  function(x, stanOptions, logC) {

    curNsamples <- ncol(logC)

    print(paste0("NOW ON ",x))
    stanOptions$bim$data$logE <- as.numeric(logC[x, ])
    stanOptions$uni$data$logE <- as.numeric(logC[x, ])

    # full data fit
    print(paste0("FULL FIT bimodal for ",x,"!"))
    print(paste0("-> USING STAN FILE: ", stanOptions$bim$file))
    # curBimFit <- do.call(stan, stanOptions$bim) ;
    curBimFit <- vb(
      object = readRDS("script/STAN/bimodal.rds"),
      data = list(
        nSamples = nSamples,
        nCells = nCells,
        nDrivers = nDrivers,
        logE = as.numeric(inDat$logC[genes[1],]),
        cell = as.integer(factor(inDat$cell, levels=unique(inDat$cell))),
        driver = as.integer(factor(inDat$driver, levels=unique(inDat$driver)))
      ),
      tol_rel_obj = 1e-4
    )
    print("-> extracting FULL FIT bimodal parameters!")
    bimPars <- inferState.getStanFitPars(
      data = stanOptions$bim$data,
      fit = curBimFit,
      mode = "bimodal"
    )
    print("-> DONE")

    print(paste0("FULL FIT unimodal for ",x,"!"))
    print(paste0("-> USING STAN FILE: ", stanOptions$uni$file))
    # curUniFit <- do.call(stan, stanOptions$uni) ;
    curUniFit <- vb(
      object = readRDS("script/STAN/unimodal.rds"),
      data = list(
        nSamples = nSamples,
        nCells = nCells,
        nDrivers = nDrivers,
        logE = as.numeric(inDat$logC[genes[1],]),
        cell = as.integer(factor(inDat$cell, levels=unique(inDat$cell))),
        driver = as.integer(factor(inDat$driver, levels=unique(inDat$driver)))
      ),
      tol_rel_obj = 1e-4
    )
    uniPars <- inferState.getStanFitPars(
      data = stanOptions$uni$data,
      fit  = curUniFit,
      mode ="unimodal"
    )

    # 1. paritionK code
    kparts <- inferState.partitionK(
      n = curNsamples,
      k = numK,
      # skipThese = c(
      #   curNsamples,
      #   curNsamples - 1),
      labels = as.integer(as.factor(inDat$driver))
    )
    bim.elpd_i <- vector(length = curNsamples, "numeric")
    uni.elpd_i <- vector(length = curNsamples, "numeric")

    for (k in 1:numK) {
      print(paste("on K partition : ",k," for gene ",x))

      hList <- which(kparts == k)
      tList <- which(kparts != k)

      stanOptions$bim.cv$data <- list(
        logE_h <- as.numeric(logC[x, hList]),
        logE_t <- as.numeric(logC[x, tList]),
        nSamples_h <- length(hList),
        nSamples_t <- length(tList)
      )

      stanOptions$uni.cv$data <- stanOptions$bim.cv$data

      print("calling bim.cv")
      # bimFit.cv <- do.call(stan, stanOptions$bim.cv) ;
      bimFit.cv <- vb(
        readRDS("script/STAN/bimodal_CV.rds"),
        data = list(
          logE_h <- as.numeric(logC[x, hList]),
          logE_t <- as.numeric(logC[x, tList]),
          nSamples_h <- length(hList),
          nSamples_t <- length(tList)
        ),
        tol_rel_obj = 1e-4
      )

      print("calling uni.cv")
      # uniFit.cv <- do.call(stan, stanOptions$uni.cv) ;
      uniFit.cv <- vb(
        readRDS("script/STAN/unimodal_CV.rds"),
        data = list(
          logE_h <- as.numeric(logC[x, hList]),
          logE_t <- as.numeric(logC[x, tList]),
          nSamples_h <- length(hList),
          nSamples_t <- length(tList)
        ),
        tol_rel_obj = 1e-4
      )

      print("extracting log lik")
      bim.kll <- extract_log_lik(bimFit.cv,
                                  parameter_name="log_lik_h")
      uni.kll <- extract_log_lik(uniFit.cv,
                                  parameter_name="log_lik_h")

      for (h in 1:length(hList)) {
        # elpd should be log column mean (logsum - logrow)...
        # but it was implemented as row means in the original code
        # ref: https://github.com/stan-dev/loo/blob/0284928f9c9966329cae8cc252be29c5f406641f/R/elpd.R#L50
        # bim.elpd_i[hList[h]] <- logSumExp(bim.kll[h,]) -
        bim.elpd_i[hList[h]] <- logSumExp(bim.kll[, h]) -
          log(nrow(bim.kll)) ;

        # uni.elpd_i[hList[h]] <- logSumExp(uni.kll[h,]) -
        uni.elpd_i[hList[h]] <- logSumExp(uni.kll[, h]) -
          log(nrow(uni.kll)) ;
      }
    }

    bim.elpd <- sum(bim.elpd_i)
    uni.elpd <- sum(uni.elpd_i)
    diff.elpd <- bim.elpd - uni.elpd
    print(paste0("bim.elpd=",bim.elpd))
    print(paste0("uni.elpd=",uni.elpd))

    bim.elpd.se  <- sqrt(curNsamples * var(bim.elpd_i))
    uni.elpd.se  <- sqrt(curNsamples * var(uni.elpd_i))
    diff.elpd.se <- sqrt(curNsamples * var(bim.elpd_i - uni.elpd_i))

    return(
      list(
        bimPars = bimPars,
        uniPars = uniPars,
        bim.elpd_i = bim.elpd_i,
        uni.elpd_i = uni.elpd_i,
        bim.elpd = bim.elpd,
        uni.elpd = uni.elpd,
        bim.elpd.se = bim.elpd.se,
        uni.elpd.se = uni.elpd.se,
        diff.elpd = diff.elpd,
        diff.elpd.se = diff.elpd.se
      )
    )
  },
  stanOptions = stanOptions,
  logC = upm
)
print(Sys.time() - start)

names(geneEfit) <- genes
curRes <- list(geneEfit = geneEfit,
                geneList = geneList,
                inDat = inDat)

# saveRDS(curRes, snakemake@output[[1]])
