# The raw mRNA files don't have some features that miRNA do,
# so the functions from the `AgiMicroRna` package need to be tweaked a bit

# Read mRNA files
readmRnaAFE = function(targets, verbose = FALSE) {
  if (!is(targets, "data.frame")) {
    stop("'targets' must be a data.frame")
  }
  dd = AgiMicroRna:::read.agiMicroRna(targets, columns = list(
    #TGS = "gTotalGeneSignal", # no feature named like this
    #TPS = "gTotalProbeSignal", # no feature named like this
    meanS = "gMeanSignal",
    medianS = "gMedianSignal", # added this one
    procS = "gProcessedSignal"),
    other.columns = list(
      #IsGeneDetected = "gIsGeneDetected", # no GENE DETECTED, CANNOT DO FILTERING!!!
      IsSaturated = "gIsSaturated",
      IsFeatNonUnifOF = "gIsFeatNonUnifOL",
      IsFeatPopnOL = "gIsFeatPopnOL",
      BGKmd = "gBGMedianSignal",
      BGKus = "gBGUsed"),
    annotation = c("ControlType", "ProbeName", "GeneName"), verbose = TRUE)
  if (verbose) {
    cat("", "\n")
    cat("  uRNAList:", "\n")
    #cat("\tdd$TGS:\t\t'gTotalGeneSignal' ", "\n")
    #cat("\tdd$TPS:\t\t'gTotalProbeSignal' ", "\n")
    cat("\tdd$meanS:\t'gMeanSignal' ", "\n")
    cat("\tdd$medianS:\t'gMedianSignal' ", "\n")
    cat("\tdd$procS:\t'gProcessedSignal' ", "\n")
    cat("", "\n")
  }
  return(dd)
}

# from `tgsMicroRna` function, tweak to use meanS
meanSignalMrna = function(dd, offset = 0, half = TRUE, verbose = TRUE) {
  if (!is(dd, "uRNAList")) {
    stop("'input' must be a uRNAList")
  }

  TT = length(dd$genes$ProbeName)
  uniqueProbe = unique(dd$genes$ProbeName)
  LUP = length(uniqueProbe)
  uniqueGene = unique(dd$genes$GeneName)
  LUG = length(uniqueGene)
  if (verbose) {
    cat("\n")
    cat("GETTING Agilent Feature Extraction Mean Signal", "\n")
    cat("\n")
    cat("\tTotal Probes:\t", TT, "\n")
    cat("\tUnique Probe: \t", LUP, "\n")
    cat("\tUnique Gene: \t", LUG, "\n")
    cat("\n")
  }
  ug = which(duplicated(dd$genes$GeneName) == FALSE)

  dd_res = dd[ug, ] # take only one gene (All the replicated genes have the same estimated meanS?)

  n_samples = nrow(dd$targets)

  # we are adjusting negative values on the Mean Signal
  # in miRNA data, this function does the correction on the TGS (Total Gene Signal)

  if (half) { # better this...
    for (i in 1:n_samples) {
      index = which(dd_res$meanS[, i] < 0.5)
      dd_res$meanS[index, i] = 0.5
    }
  }
  else { # than this
    min = min(dd_res$meanS)
    for (i in 1:n_samples) {
      dd_res$meanS[, i] = dd_res$meanS[, i] + (abs(min) + offset)
    }
  }

  return(dd_res)
}

# from `tgsNormalization` function, tweak to use `meanS`
meanSignal_Normalization = function (dd, NORMmethod = "quantile", targets, verbose = FALSE) {

  if (!is(dd, "uRNAList")) {
    stop("'input' must be a uRNAList")
  }

  if (NORMmethod != "none" && NORMmethod != "quantile" && NORMmethod !=
      "scale") {
    stop("NORMmethod should be one of 'none', 'quantile','scale'")
  }

  if (NORMmethod != "none") {
    ddNORM = dd
    exprsNORM = limma::normalizeBetweenArrays(dd$meanS, method = NORMmethod)
    if (NORMmethod == "quantile") {
      ddNORM$meanS = log2(exprsNORM)
    }
    else {
      ddNORM$meanS = exprsNORM
    }
    rm(exprsNORM)
  }
  else { # if no normalization, just take the log2
    ddNORM = dd
    ddNORM$meanS = log2(dd$meanS)
  }
  if (verbose) {
    cat("------------------------------------------------------",
      "\n")
    cat("\tNORMMALIZATION:\t", NORMmethod, "\n")
    cat("\tOUTPUT in log-2 scale", "\n")
    cat("------------------------------------------------------",
      "\n")
  }

  return(ddNORM)
}

















