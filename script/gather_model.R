#!/usr/bin/env Rscript
obj_list <- lapply(snakemake@input[['indi_mod']], readRDS)

## Load metadata csv
metadata <- read.csv(snakemake@input[['meta']], row.names = 1)
if (!"pseudobulk" %in% colnames(metadata)) {
  metadata$pseudobulk <- 1
}

geneEfit <- lapply(
  obj_list, function(x) {
    x[["geneEfit"]]
  }
)

geneEfit <- Reduce(append, geneEfit)

geneList <- lapply(
  obj_list, function(x) {
    x[["geneList"]]
  }
)

geneList <- Reduce(c, geneList)

inDat <- lapply(
  obj_list, function(x) {
    x[["inDat"]]
  }
)

O <- lapply(
  inDat, function(x) {x[["O"]]}
)

O <- do.call(rbind, O)

cell <- lapply(
  inDat, function(x) {x[["cell"]]}
)

cell <- cell[[1]]

merged_inDat <- list(
  nGenes = nrow(O),
  nCells = length(unique(metadata$cell)),
  nSamples = length(unique(metadata$sample)),
  nDrivers = length(unique(metadata$pseudobulk)),
  O = O,
  logC = O,
  cell = cell,
  driver = metadata$pseudobulk,
  sample = colnames(O)
)

curRes <- list(geneEfit = geneEfit,
               geneList = geneList,
               inDat = merged_inDat)
saveRDS(curRes, snakemake@output[[1]])