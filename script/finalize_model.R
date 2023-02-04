#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(loo))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(txtplot))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(R.utils))
source("script/utils.R")

curRes <- readRDS(snakemake@input[['main_mod']])
## Load metadata csv
metadata <- read.csv(snakemake@input[['meta']], row.names = 1)

if (!"pseudobulk" %in% colnames(metadata)) {
  metadata$pseudobulk <- 1
}

## Extract values to be modeled from the input object
## This is more meaningful as Davis et al. used the object for other
## purposes.
## In my case, keeping this is just to avoid massively changing the data
## processing logic of the modeling code.
upm <- curRes$inDat$O
inDat <- curRes$inDat
geneEfit <- curRes$geneEfit
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

## Defining probability matrices
# Goal is to populate these probability matrices:
#### bimodal probability matrices
## Cell-level
p_on_bimodal_gc <- matrix(nrow = nGenes,
                          ncol = nCells,
                          dimnames = list(genes, unique(cells)))
## Driver (pseudobulk)-level
p_on_bimodal_gd <- matrix(nrow = nGenes,
                          ncol = nDrivers,
                          dimnames = list(genes,unique(drivers)))
## Sample-level (individual columns)
p_on_bimodal_gs <- matrix(nrow = nGenes,
                          ncol = nSamples,
                          dimnames = list(genes,samples))
#### Metrics to compare the fit of bimodal vs unimodal model
### elpd is a Bayesian metric for the fit from cross validation.
## Differences and the standard error of differences of elpds
diff_elpd_g <- setNames(vector(length = nGenes), genes)
diff_elpd_se_g <- setNames(vector(length = nGenes), genes)

#### Log probability: Calculated by taking logarithm of the p_* matrices
logp_on_unimodal_g <- setNames(vector(length = nGenes), genes)
logp_on_gs <- matrix(nrow = nGenes,
                      ncol = nSamples,
                      dimnames = list(genes,samples))
logp_on_gd <- matrix(nrow = nGenes,
                      ncol = nDrivers,
                      dimnames = list(genes,unique(drivers)))
logp_on_gc <- matrix(nrow = nGenes,
                      ncol = nCells,
                      dimnames = list(genes,unique(cells)))

## Fill result table and print warnings if elpd is missing
for (g in seq(nGenes)) {
  if (! "diff.elpd" %in% names(geneEfit[[genes[g]]])) {
    print(paste0("SKIPPING: ",geneList[genes[g]]," since missing diff.elpd esitmate"))
    next;
  }
  diff_elpd_g[g] <- geneEfit[[genes[g]]]$diff.elpd
  diff_elpd_se_g[g] <- geneEfit[[genes[g]]]$diff.elpd.se
  p_on_bimodal_gc[g,] <- exp(geneEfit[[genes[g]]]$bimPars$p_on_bimodal_gc)
  p_on_bimodal_gd[g,] <- exp(geneEfit[[genes[g]]]$bimPars$p_on_bimodal_gd)
  p_on_bimodal_gs[g,] <- exp(geneEfit[[genes[g]]]$bimPars$p_on_bimodal_gs)
}

## Part 2. Model global distributions of p(E|off) and p(E|on)
## from bimodal genes (diff_elpd > 2 * diff_elpd_se)

t.bimGenes <- names(diff_elpd_g)[which(diff_elpd_g > 2 * diff_elpd_se_g)]
print(paste0("-> n bimodal genes=",length(t.bimGenes)))
bim.logC <- logC[t.bimGenes,]
bim.p_on_gc <- p_on_bimodal_gc[t.bimGenes,]

coarse.p_on_bimodal_gs <- matrix(nrow = nGenes,
                                  ncol = nSamples,
                                  dimnames = list(genes,samples))

## TODO: This might be replaced with
# coarse.p_on_bimodal_gs <- p_on_bimodal_gc[ , inDat$cell]
for (s in 1:nSamples) {
  coarse.p_on_bimodal_gs[,s] <- p_on_bimodal_gc[,inDat$cell[s]]
}
###############################################################

bim.p_on_gs <- coarse.p_on_bimodal_gs[t.bimGenes, , drop = FALSE]

real.onE <- as.vector(
  unlist(
    lapply(
      1:nrow(bim.logC),
      function(x, t1, t2){
        t1[x, which(t2[x, ] > 0.9)]
      },
      t1 = bim.logC, t2 = bim.p_on_gs
    )
  )
)

txtdensity(real.onE, width=80)

real.offE <- as.vector(
  unlist(
    lapply(
      1:nrow(bim.logC),
      function(x,t1,t2){
        t1[x,which(t2[x, ] < 0.1)]
      },
      t1 = bim.logC, t2 = bim.p_on_gs
    )
  )
)

txtdensity(real.offE, width=80)

print("FITTING ON DISTRIBUTION")
uniEonfit <- Mclust(real.onE,G = 1)
print("- >DONE")

uniEonfit$para$sd <- sqrt(uniEonfit$para$variance$sigmasq)

print("FITTING off DISTRIBUTION")
uniEofffit <- Mclust(real.offE, G = 1)
print("- >DONE")

uniEofffit$para$sd <- sqrt(uniEofffit$para$variance$sigmasq)
pi_on.uni <- length(real.onE) / (length(real.onE) + length(real.offE))
logpi_on.uni <- log(pi_on.uni)
logpi_off.uni <- log1p(-1 * pi_on.uni)

print(paste0("uniEon: mean=",uniEonfit$para$mean," sd=",uniEonfit$para$sd))
print(paste0("uniEoff: mean=",uniEofffit$para$mean," sd=",uniEofffit$para$sd))

## Part 2b.  Estimate p(on|unimodal) for all gene
i <- 1
for (g in geneList) {
  logp_numer <- 0 ; logp_denom <- 0 ;
  for (s in samples) {
    t.on <- logpi_on.uni + dnorm(logC[g, s],
                                  mean = uniEonfit$para$mean,
                                  sd   = uniEonfit$para$sd,
                                  log  = TRUE)
    logp_numer <- logp_numer + t.on

    t.off <- logpi_off.uni + dnorm(logC[g,s],
                                    mean    = uniEofffit$para$mean,
                                    sd      = uniEofffit$para$sd,
                                    log     = TRUE)
    logp_denom <- logp_denom + t.off

  }
  logp_on_unimodal_g[g] <- logp_numer - logSumExp(c(logp_denom, logp_numer))
  if (i %% 500 == 0) { print(paste0("done with ",i," genes")) }
  i <- i + 1
}

## Part 3. Estimate P(on_gs) by choosing unimodal or bimodal
i <- 1
for (g in geneList) {

  if (diff_elpd_g[g] > 0 &
      diff_elpd_g[g] > 2 * diff_elpd_se_g[g]) {

    geneEfit[[g]]$exprType <- "bimodal"

    logp_on_gs[g,] <- log(p_on_bimodal_gs[g,]) ;
    logp_on_gd[g,] <- log(p_on_bimodal_gd[g,]) ;
    logp_on_gc[g,] <- log(p_on_bimodal_gc[g,]) ;

  } else {

    geneEfit[[g]]$exprType <- "unimodal"

    logp_on_gs[g,] <- logp_on_unimodal_g[g] ;
    logp_on_gd[g,] <- logp_on_unimodal_g[g] ;
    logp_on_gc[g,] <- logp_on_unimodal_g[g] ;

  }
  if (i %% 500 == 0) { print(paste0("done with ",i," genes")) }
  i <- i + 1
}

final_model <- list(
  uniFit = list(uniEonfit   = uniEonfit,
                uniEofffit  = uniEofffit,
                pi_on.uni   = pi_on.uni),
  inDat = inDat,
  logC = logC,
  geneEfit = geneEfit,
  p_on_bimodal_gc = p_on_bimodal_gc,
  p_on_bimodal_gd = p_on_bimodal_gd,
  p_on_bimodal_gs = p_on_bimodal_gs,
  logp_on_unimodal_g = logp_on_unimodal_g,
  logp_on_gc = logp_on_gc,
  logp_on_gd = logp_on_gd,
  logp_on_gs = logp_on_gs
)

## Save model
saveRDS(
  final_model,
  snakemake@output[['model']]
)

# Part 4. Making the table #
Modeled_expression = matrix(nrow = length(genes), ncol = nCells + 9)
row.names(Modeled_expression) = genes
colnames(Modeled_expression) = c(unique(inDat$cell),
                                  "diff.elpd",
                                  "diff.elpd.se",
                                  "bimodal.on_level",
                                  "bimodal.off_level",
                                  "bimodal.sd",
                                  "bimodal.pi",
                                  "unimodal.level",
                                  "unimodal.sd",
                                  "exprType"
)

for (gene in rownames(Modeled_expression)) {
  if (geneEfit[[gene]]$diff.elpd <= geneEfit[[gene]]$diff.elpd.se*2) {
    #gene is unimodal
    Modeled_expression[gene,] = c(round(exp(logp_on_gc[gene,]), digits = 2),
                                  round(final_model$geneEfit[[gene]]$diff.elpd, digits = 2),
                                  round(final_model$geneEfit[[gene]]$diff.elpd.se, digits = 2),
                                  round(final_model$geneEfit[[gene]]$bimPars$mu_on, digits = 2),
                                  round(final_model$geneEfit[[gene]]$bimPars$mu_off, digits = 2),
                                  round(final_model$geneEfit[[gene]]$bimPars$sd_on, digits = 2),
                                  round(final_model$geneEfit[[gene]]$bimPars$pi_on, digits = 2),
                                  round(final_model$geneEfit[[gene]]$uniPars$mu1, digits = 2),
                                  round(final_model$geneEfit[[gene]]$uniPars$sd1, digits = 2),
                                  geneEfit[[gene]]$exprType)
  } else {
    #gene is bimodal
    Modeled_expression[gene,] = c(round(p_on_bimodal_gc[gene,], digits = 2),
                                  round(final_model$geneEfit[[gene]]$diff.elpd, digits = 2),
                                  round(final_model$geneEfit[[gene]]$diff.elpd.se, digits = 2),
                                  round(final_model$geneEfit[[gene]]$bimPars$mu_on, digits = 2),
                                  round(final_model$geneEfit[[gene]]$bimPars$mu_off, digits = 2),
                                  round(final_model$geneEfit[[gene]]$bimPars$sd_on, digits = 2),
                                  round(final_model$geneEfit[[gene]]$bimPars$pi_on, digits = 2),
                                  round(final_model$geneEfit[[gene]]$uniPars$mu1, digits = 2),
                                  round(final_model$geneEfit[[gene]]$uniPars$sd1, digits = 2),
                                  geneEfit[[gene]]$exprType)
  }
}


saveRDS(object = Modeled_expression,
        file = snakemake@output[["prob_mtx"]])
