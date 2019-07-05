initFiles1 <- function (pattern = "Chrom", snpInfoDir, signalFile, mc.cores = 1) {
  if (missing(snpInfoDir) | missing(signalFile)) {
    stop('[Enrichment:initFiles] argument(s) missing.', call. = FALSE)
  } else {}
  snpInfoDir <- .checkFilePath(snpInfoDir)
  FILES <- list.files(snpInfoDir, pattern = ".bim")
  .checkSnpInfoDir(snpInfoDir)
  .checkSignalFile(signalFile)
  tmpDir <- gsub("\\\\", "/", tempdir())
  dir.create(paste0(tmpDir, "/snpEnrichment/"), showWarnings = FALSE)
  cat("All files are ready for chromosome:\n  ")
  resParallel <- mclapply2(X = seq(22), mc.cores = min(22, mc.cores), FUN = function (iChr) {
    # newPattern <- unlist(strsplit(grep(paste0(pattern, iChr, "[^0-9]*.bim"), FILES, value = TRUE), ".bim"))[1]
    newPattern <- gsub(".bim", "", grep(paste0(pattern, iChr, "[^0-9]"), FILES, value = TRUE))
    err1 <- try(.writeSignal(pattern = newPattern, snpInfoDir = snpInfoDir, signalFile = signalFile), silent = TRUE)
    err2 <- try(.writeFreq(pattern = newPattern, snpInfoDir = snpInfoDir), silent = TRUE)
    cat(iChr, " ", sep = "")
    if (class(err1)=="try-error" | class(err2)=="try-error") {
      return(invisible("ERROR"))
    } else {
      return(invisible())
    }
  })
  # if (any(unlist(resParallel)=="ERROR")) {
  #   stop("[Enrichment:initFiles] initialize files failed.", call. = FALSE)
  # } else {}
  cat("\n\n")
  return(resParallel)
}

library(parallel)
setGeneric(name = "enrichSNP", def = function (List, Table, EnrichmentRatio, Z, PValue, Resampling) {standardGeneric("enrichSNP")})
setGeneric(name = "chromosome", def = function (Data, LD, eSNP, xSNP) {standardGeneric("chromosome")})
setGeneric(name = "enrichment", def = function (Loss, Call, eSNP, xSNP, Chromosomes) {standardGeneric("enrichment")})
setGeneric(name = "doLDblock", def = function (object, mc.cores = 1) {standardGeneric("doLDblock")})
setGeneric(name = "excludeSNP", def = function (object, excludeFile, mc.cores = 1) {standardGeneric("excludeSNP")})
setGeneric(name = "computeER", def = function (object, sigThresh = 0.05, mc.cores = 1) {standardGeneric("computeER")})
setGeneric(name = "reSample", def = function (object, nSample = 100, empiricPvalue = TRUE, MAFpool = c(0.05, 0.10, 0.2, 0.3, 0.4, 0.5), mc.cores = 1, onlyGenome = TRUE, ...) {standardGeneric("reSample")})
setGeneric(name = "compareEnrichment", def = function (object.x, object.y, pattern = "Chrom", nSample = 100, empiricPvalue = TRUE, mc.cores = 1, onlyGenome = TRUE) {standardGeneric("compareEnrichment")})
setGeneric(name = "is.enrichment", def = function (object) {standardGeneric("is.enrichment")})
setGeneric(name = "is.chromosome", def = function (object) {standardGeneric("is.chromosome")})
setGeneric(name = "is.EnrichSNP", def = function (object) {standardGeneric("is.EnrichSNP")})
setGeneric(name = "reset", def = function (object, i) {standardGeneric("reset")})
setGeneric(name = "getEnrichSNP", def = function (object, type = "eSNP") {standardGeneric("getEnrichSNP")})

setMethod(f = "reSample", signature = "ANY", definition = function (object, nSample, sigThresh, MAFpool, mc.cores) {
  if (!(is.enrichment(object) & is.chromosome(object))) {
    stop('[Method:reSample] not available for "', class(object), '" object.', call. = FALSE)
  } else {}
})
setMethod(f = "excludeSNP", signature = "ANY", definition = function (object, excludeFile, mc.cores = 1) {
  if (!is.enrichment(object)) {
    stop('[Method:excludeSNP] not available for "', class(object), '" object.', call. = FALSE)
  } else {}
})
setMethod(f = "reset", signature = "ANY", definition = function (object, i) {
  if (!(is.enrichment(object) & is.chromosome(object))) {
    stop('[Method:reset] not available for "', class(object), '" object.', call. = FALSE)
  } else {}
})
setMethod(f = "getEnrichSNP", signature = "ANY", definition = function (object, type = "eSNP") {
  if (!(is.enrichment(object))) {
    stop('[Method:getEnrichSNP] not available for "', class(object), '" object.', call. = FALSE)
  } else {}
})


.verbose <- function (expr) {return(invisible(capture.output(expr)))}


GC <- function (verbose = getOption("verbose"), reset = FALSE) {
  while (!identical(gc(verbose, reset)[, 4], gc(verbose, reset)[, 4])) {}
  return(gc(verbose, reset))
}


.checkFilePath <- function (path) {
  # END <- unlist(regmatches(path, regexec(".$", path)))
  # START <- unlist(regmatches(path, regexec("^.", path)))
  # if (START != "/") {
  # path <- paste0("/", path)
  # } else {}
  # if (END != "/") {
  # path <- paste0(path, "/")
  # } else {}
  path <- gsub("/*$", "/", path)
  return(path)
}


maxCores <- function (mc.cores = 1) {
  if (Sys.info()[["sysname"]] == "Linux") {
    nbCores <- detectCores()
    mc.cores.old <- mc.cores
    if (file.exists("/proc/meminfo")) {
      memInfo <- readLines("/proc/meminfo")
      sysMemFree <- memInfo[grep('^MemFree:', memInfo)]
      sysMemCached <- memInfo[grep('^Cached:', memInfo)]
      sysMemAvailable <- 0.95*(as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", sysMemFree)) + as.numeric(gsub("[^0-9]*([0-9]*)", "\\1", sysMemCached)))
      sysProc <- as.numeric(unlist(strsplit(system(paste("ps v", Sys.getpid()), intern = TRUE)[2], " +"))[8])
      mc.cores <- max(min(as.numeric(mc.cores), floor(sysMemAvailable/sysProc)), 1)
      if (mc.cores > nbCores) {
        mc.cores <- nbCores
      } else {}
      if (mc.cores != mc.cores.old) {
        warning(paste0('To avoid memory overload "mc.cores" was decreased to "', mc.cores, '".'), call. = FALSE)
      } else {}
    } else {
      mc.cores <- ifelse(mc.cores.old>nbCores, nbCores, mc.cores.old)
    }
  } else {
    mc.cores <- 1
  }
  return(mc.cores)
}


mclapply2 <- function (X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE, mc.allow.recursive = TRUE) {
  if (Sys.info()[["sysname"]] != "Linux") {
    mc.cores <- 1
  } else {
    mc.cores <- min(detectCores(), mc.cores)
  }
  return(mclapply(X = X, FUN = FUN, ...,
                  mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.silent = mc.silent,
                  mc.cores = maxCores(mc.cores), mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive))
}


.writeSignal <- function (pattern, snpInfoDir, signalFile) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  if (length(unlist(strsplit(readLines(signalFile, n = 1), split = "\t")))>1) {
    signal <- read.delim(file = signalFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                         colClasses = c("character", "numeric"), na.string = c("NA", ""),
                         check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"))
  } else {
    if (length(unlist(strsplit(readLines(signalFile, n = 1), split = " ")))>1) {
      signal <- read.delim(file = signalFile, header = TRUE, sep = " ", stringsAsFactors = FALSE,
                           colClasses = c("character", "numeric"), na.string = c("NA", ""),
                           check.names = FALSE, strip.white = TRUE, col.names =  c("SNP", "PVALUE"))
    } else {
      stop('[Enrichment:initFiles] only " " and "\t" are allowed as columns separator in Signal file.', call. = FALSE)
    }
  }
  
  chrom.bim <- read.delim(file = paste0(snpInfoDir, pattern, ".bim"), header = FALSE, stringsAsFactors = FALSE,
                          colClasses = c("numeric", "character", "NULL", "NULL", "NULL", "NULL"), na.string = c("NA", ""),
                          check.names = FALSE, strip.white = TRUE, col.names =  c("CHR", "SNP", "", "", "", ""))
  signal <- merge(signal, chrom.bim, by = "SNP")[, c(3, 1, 2)]
  eval(parse(text = paste0('write.table(signal, file = "', tmpDir, '/snpEnrichment/', pattern,
                           '.signal", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")')))
  return(invisible())
}


.readSignal <- function (pattern) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  signal <- read.delim(file = paste0(tmpDir, "/snpEnrichment/", pattern, ".signal"), header = TRUE, stringsAsFactors = FALSE,
                       colClasses = c("NULL", "character", "numeric"), na.string = c("NA", ""),
                       check.names = FALSE, strip.white = TRUE, col.names =  c("", "SNP", "PVALUE"))
  return(signal)
}


.writeFreq <- function (pattern, snpInfoDir) {
  tmpDir <- gsub("\\\\", "/", tempdir())
  IN <- paste0(snpInfoDir, pattern)
  OUT <- paste0(tmpDir, "/snpEnrichment/", pattern)
  plinkData <- read.plink(bed = paste0(IN, ".bed"), bim = paste0(IN, ".bim"), fam = paste0(IN, ".fam"), select.snps = .readSignal(pattern)[, "SNP"])
  plinkFreq <- col.summary(plinkData$genotypes)
  plinkFreq <- cbind(snp.name = rownames(plinkFreq), MAF = plinkFreq[, "MAF"])
  plinkRes <- merge(plinkData$map, plinkFreq, by = "snp.name")
  plinkRes <- plinkRes[, c("chromosome", "snp.name", "position", "MAF")]
  plinkRes[, "MAF"] <- as.numeric(as.character(plinkRes[, "MAF"]))
  colnames(plinkRes) <- c("CHR", "SNP", "POS", "MAF")
  write.table(plinkRes, paste0(OUT, ".all"), row.names = FALSE, sep = "\t")
  return(invisible())
}


.checkSnpInfoDir <- function (snpInfoDir) {
  snpInfoDir <- .checkFilePath(snpInfoDir)
  if (length(list.files(snpInfoDir, pattern = "*.bim"))!=22) {
    stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.bim")), " 'bim' files found when 22 is needed."), call. = FALSE)
  } else {}
  if (length(list.files(snpInfoDir, pattern = "*.bed"))!=22) {
    stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.bed")), " 'bed' files found when 22 is needed."), call. = FALSE)
  } else {}
  if (length(list.files(snpInfoDir, pattern = "*.fam"))!=22) {
    stop(paste0("[Enrichment:initFiles] Only ", length(list.files(snpInfoDir, pattern = "*.fam")), " 'fam' files found when 22 is needed."), call. = FALSE)
  } else {}
  return(invisible())
}


.checkSignalFile <- function (signalFile) {
  if (!file.exists(signalFile)) {
    stop(paste0("[Enrichment:initFiles] ", signalFile, " doesn't exist."), call. = FALSE)
  } else {}
  return(invisible())
}


.checkSnpListDir <- function (snpListDir, pattern) {
  snpListDir <- .checkFilePath(snpListDir)
  if (length(list.files(snpListDir, pattern = pattern))==0) {
    stop(paste0("[Enrichment:readEnrichment] No snp list file found when at least one is needed."), call. = FALSE)
  } else {}
  return(invisible())
}
