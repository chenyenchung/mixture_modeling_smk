#!/usr/bin/env Rscript

split_num <- 200
## Can I get split num from snakemake? No.
## splitgather is not exposed in the snakemake object

main <- read.csv(snakemake@input[[1]], row.names = 1)

binsize <- round(nrow(main) / split_num, digits = 0) + 1
bin_id <- rep(1:split_num, binsize)[seq(nrow(main))]
main <- split(main, bin_id)

print("Done splitting")

for (i in seq(split_num)) {
# print(i)
  write.csv(main[[i]],
            paste0("data/", snakemake@wildcards[["sample"]], "/splitted/", i, "-of-", split_num, ".csv")
            )
}
