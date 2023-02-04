inferState.partitionK <- function(n, k, labels, skipThese) {
  # Purpose: Partitions a set of points for k-fold cross-validation
  #
  # returns: array of k-partition assignment

  fullList <- 1:n
  if (!missing(skipThese)) {
    fullList <- setdiff(fullList, skipThese)
  } else {
    skipThese <- c()
  }

  kparts <- vector(mode="numeric",length=n)

  shuffledI <- sample(fullList, length(fullList))
  medSize <- floor(length(fullList) / k)
  for (i in 1:k) {
    t.start <- ((i - 1) * medSize) + 1
    for (j in t.start:(t.start + medSize - 1)) {
      kparts[shuffledI[j]] <- i
    }
  }

  if (length(fullList) %% k > 0) {
    for (j in (k * medSize + 1):length(fullList)) {
      kparts[shuffledI[j]] <- k
    }
  }


  # If points of a particuilar label are all found in same parition, shuffle with
  # another class point to reblanace

  keepLabels <- labels[setdiff(fullList, skipThese)]

  allGood <- FALSE
  while (!allGood) {
    allGood <- TRUE
    for (curLabel in unique(keepLabels)) {
      parts <- unique(kparts[labels == curLabel])
      if (length(parts) == 1) {

        # swap 2 points: 1 from this label, 2 from another label, in another class.
        i <- sample(which(labels == curLabel))[1]
        j <- sample(which(labels != curLabel & kparts != parts[1]))[1]
        part.t <- kparts[i]
        kparts[i] <- kparts[j]
        kparts[j] <- part.t

        allGood <- FALSE
      }
    }
  }

  for (i in skipThese) {
    kparts[i] <- k + 1 }

  return(kparts)

}

inferState.getStanFitPars <- function(data, fit, mode) {
  # Purpose: extracts STAN fit parameters

  curPars <- list()

  if (mode == "cv") {

    curPars$log_lik_t <- extract_log_lik(fit, parameter_name="log_lik_t")
    curPars$log_lik_h <- extract_log_lik(fit, parameter_name="log_lik_h")

  } else if (mode == "bimodal") {

    stansum <- summary(fit)$summary ;
    #      curPars$fullSum <- stansum

    curPars$p_on_bimodal_gc           <- stansum[
      paste0("pon_gc[",1:data$nCells,"]"),"50%"]
    curPars$p_on_bimodal_gd           <- stansum[
      paste0("pon_gd[",1:data$nDrivers,"]"),"50%"]
    curPars$p_on_bimodal_gs           <- stansum[
      paste0("pon_gs[",1:data$nSamples,"]"),"50%"]
    curPars$mu_on                     <- stansum["mu_on","50%"]
    curPars$mu_off                    <- stansum["mu_off","50%"]
    curPars$sd_on                     <- stansum["sd_on","50%"]
    curPars$sd_off                    <- stansum["sd_off","50%"]
    curPars$pi_on                     <- stansum["pi_on","50%"]

    curPars$rhat.p_on_bimodal_gc      <- stansum[
      paste0("pon_gc[",1:data$nCells,"]"),"Rhat"]
    curPars$rhat.p_on_bimodal_gd      <- stansum[
      paste0("pon_gd[",1:data$nDrivers,"]"),"Rhat"]
    curPars$rhat.p_on_bimodal_gs      <- stansum[
      paste0("pon_gs[",1:data$nSamples,"]"),"Rhat"]
    curPars$rhat.mu_on                <- stansum["mu_on","Rhat"]
    curPars$rhat.mu_off               <- stansum["mu_off","Rhat"]
    curPars$rhat.sd_on                <- stansum["sd_on","Rhat"]
    curPars$rhat.sd_off               <- stansum["sd_off","Rhat"]
    curPars$rhat.pi_on                <- stansum["pi_on","Rhat"]

    curPars$neff.p_on_bimodal_gc      <- stansum[
      paste0("pon_gc[",1:data$nCells,"]"),"n_eff"]
    curPars$neff.p_on_bimodal_gd      <- stansum[
      paste0("pon_gd[",1:data$nDrivers,"]"),"n_eff"]
    curPars$neff.p_on_bimodal_gs      <- stansum[
      paste0("pon_gs[",1:data$nSamples,"]"),"n_eff"]
    curPars$neff.mu_on                <- stansum["mu_on","n_eff"]
    curPars$neff.mu_off               <- stansum["mu_off","n_eff"]
    curPars$neff.sd_on                <- stansum["sd_on","n_eff"]
    curPars$neff.sd_off               <- stansum["sd_off","n_eff"]
    curPars$neff.pi_on                <- stansum["pi_on","n_eff"]

  } else if (mode == "unimodal") {

    stansum <- summary(fit)$summary ;
    print("DEBUGGING unimodal par extraction!!")
    print(rownames(stansum))

    curPars$mu1       <- stansum["mu1","50%"]
    curPars$sd1       <- stansum["sd1","50%"]
    curPars$rhat.mu1  <- stansum["mu1","Rhat"]
    curPars$rhat.sd1  <- stansum["sd1","Rhat"]

    curPars$neff.mu1  <- stansum["mu1","n_eff"]
    curPars$neff.sd1  <- stansum["sd1","n_eff"]

  }

  return(curPars)
}

show_help <- function() {
  message("Please provide the following options:")
  message("\t--name\t\tA unique name for this run. Used to decide where to save results.")
  message("\t--expMtx\tPath to the expression matrix to analyze. Rows are genes; columns are samples.")
  message("\t--metadata\t(Optional) Path to a metadata table.\n\t\t\tMust contain 2 columns: sample: Sample names and cell: Cluster identity.")
  message("\t--nIter\t\t(Optional) number of iterations for Stan. Default is 500")
  message("\t--numK\t\t(Optional) number of splits for cross-validation. Default is 10")
  message("\t--genes\t\tA txt file in which each line is a gene symbol. If unset, all genes will be analyzed.")
  message("\t--test\t\tOnly run the first 2 genes.")
  message("\t--start\t\tStart the analysis from the n-th gene.")
  message("\t--end\t\tEnd the analysis from the n-th gene.")
  # print(str(cmd_arg))
}
