SAMPLE = ['dsx', 'dsx_sp']

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
