#!/usr/bin/env Rscript
obj_list <- lapply(snakemake@input[['indi_mod']], readRDS)

message("Done loading objects.")
metadata <- read.csv(snakemake@input[['meta']], row.names = 1)
if (!"pseudobulk" %in% colnames(metadata)) {
  metadata$pseudobulk <- 1
}

message("Done loading metadata")
print(head(metadata))

geneEfit <- lapply(
  obj_list, function(x) {
    x[["geneEfit"]]
  }
)

geneEfit <- Reduce(append, geneEfit)
message("Done gathering fit objects")

geneList <- lapply(
  obj_list, function(x) {
    x[["geneList"]]
  }
)

geneList <- Reduce(c, geneList)

message("Done gathering gene list")

inDat <- lapply(
  obj_list, function(x) {
    x[["inDat"]]
  }
)

O <- lapply(
  inDat, function(x) {x[["O"]]}
)

O <- do.call(rbind, O)


message("Done gathering expression matrices")

cell <- lapply(
  inDat, function(x) {x[["cell"]]}
)

cell <- cell[[1]]


message("Done gathering cell barcodes")
message("Constructing output object")

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

message("Begin exporting final object")
saveRDS(curRes, snakemake@output[[1]])
