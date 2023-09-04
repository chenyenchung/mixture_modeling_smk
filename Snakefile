## This pipeline expects you to provide lognormalized expression matrix
## in which each row is a gene while each column is a cluster.
##
## Additionally, you need a two-column csv (with row.names but is not
## used by the script)
##
## - Column "cell" is the cell type that you annotate.
##   This corresponds to individual cell types in Davis et al.
## - Column "sample" can be optionally set as libraries if you prefer
##   to split a cluster by libraries (or any pseudobulking you prefer.
##   I am still figuring out how pseudobulking influences model 
##   performance).
##
## File structure expected under your working directory:
## 
## - result: For storing results
## - int: Intermediate files of less importance. Some are noted
##        by temp() and will be removed by Snakemake automatically.
## - data: Under which each directory is expected to be correspond
##         to a dataset and two csv (see above) files are required
##         to run the pipeline.
##
## - script: Scripts that are called by Snakemake. Can be used on
##           there own but will need some adjustments.
## - log: Running logs
##
## Running environment:
##
## Snakemake should be able to use whatever that is in your $PATH.
## You will need Stan and Rstan and a couple of R packages.
## (Please see the R scripts under script/).
## For each step in Snakemake (rule xxx), the envmodules element
## is LMOD modules on NYU HPC. You might need to modify them to
## corresponding ones in your environment. If you don't use
## LMOD modules, they are safe to remove if you have those
## softwares installed.


## Snakemake uses Python syntax, and the directories under data/
## to run is defined by the following list (SAMPLE).
## (i.e., if you run the pipeline as is now, it will go to
## data/ol_atlas and try to run the whole pipeline.
## Multiple directories can be passed and processed together
## (e.g., ['dir1', 'dir2'...])

SAMPLE = ['ol_atlas']


scattergather:
  split = 200
rule all:
  input: expand("result/{sample}/fit_model.rds", sample=SAMPLE)

rule split:
  input: "data/{sample}/lognorm.csv",
  output: temp(scatter.split("data/{{sample}}/splitted/{scatteritem}.csv"))
  envmodules: "r/gcc/4.2.0"
  resources:
    mem_mb = 2000,
    runtime =20 
  threads: 1
  script: "script/split_n.R"
    
rule stan_model:
  input:
    mtx = "data/{sample}/splitted/{scatteritem}.csv",
    meta = "data/{sample}/metadata.csv"
  output: temp("int/{sample}_{scatteritem}/stan_model.rds")
  envmodules: "r/gcc/4.2.0"
  params:
    nIter = 500,
    numK = 10,
  resources:
    mem_mb = 12000,
    runtime = 180
  threads: 2
  script: "script/mixture_modeling.R"

rule gather_out:
  input:
    indi_mod = gather.split("int/{{sample}}_{scatteritem}/stan_model.rds"),
    meta = "data/{sample}/metadata.csv"
  output: "int/{sample}/stan_main.rds"
  envmodules: "r/gcc/4.2.0"
  resources:
    mem_mb = 2000,
    runtime = 20
  threads: 1
  script: "script/gather_model.R"

rule finalize_model:
  input:
    main_mod = rules.gather_out.output,
    meta = "data/{sample}/metadata.csv"
  output:
    prob_mtx = "result/{sample}/modeled_prob.rds",
    model = "result/{sample}/fit_model.rds"
  envmodules: "r/gcc/4.2.0"
  resources:
    mem_mb = 8000,
    runtime = 60
  threads: 1
  script: "script/finalize_model.R"
