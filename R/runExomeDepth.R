#' @title Run ExomeDepth on mcCNV counts object
#' @description Wrapper to run the ExomeDepth algorithm on a
#' [counts object][validObjects]
#' @param counts data.table [counts object][validObjects]
#' @param transProb passed to [ExomeDepth::CallCNVs] 'transition.probability'
#' @param cnvLength passed to [ExomeDepth::CallCNVs] 'expected.CNV.length'
#' @param ... arguments passed to [mclapply][parallel::mclapply]
#'
#' @details
#' Runs ExomeDepth using the default parameters, then maps the call information
#' back to the counts object.
#' Will use [parallel::mclapply()] & [parallel::mcmapply()] to parallelize
#' the computation when available.
#'
#' @import data.table
#' @import mcCNV
#' @importFrom ExomeDepth select.reference.set CallCNVs
#' @importClassesFrom ExomeDepth ExomeDepth
#' @export

runExomeDepth <- function(counts, transProb = 1e-4, cnvLength = 5e4) {

  cmat <- cnvCountsToMatrix(counts)
  int <- unique(counts[ , .(seqnames, start, end)])
  int[ , intName := sprintf("%s:%d-%d", seqnames, start, end)]
  setkey(int, seqnames, start, end)
  sbjVec <- colnames(cmat)

  getRef <- function(sbj) {
    try(select.reference.set(test.counts = cmat[ , sbj],
                             reference.counts = cmat[ , setdiff(sbjVec, sbj)],
                             bin.length = int$end - int$start + 1,
                             n.bins.reduced = min(1e4, nrow(cmat))))
  }

  refList <- lapply(sbjVec, getRef)
  names(refList) <- sbjVec

  failed <- sapply(refList, is, 'try-error')
  if (all(failed)) stop("ExomeDepth failed for all samples.")
  if (any(failed)) {
    refList <- refList[!failed]
    warning('ExomeDepth failed for: ', setdiff(sbjVec, names(refList)))
    sbjVec <- sbjVec[!failed]
  }

  calcCN <- function(sbj) {
    ref <- rowSums(cmat[ , refList[[sbj]]$reference.choice, drop = FALSE])
    cn <- new('ExomeDepth',
              test = cmat[ , sbj],
              reference = ref,
              formula = 'cbind(test, reference) ~ 1')
    cn <- CallCNVs(x = cn,
                   chromosome = int$seqnames,
                   start = int$start,
                   end = int$end,
                   name = int$intName,
                   transition.probability = transProb,
                   expected.CNV.length = cnvLength)
    cn
  }

  cnList <- lapply(sbjVec, calcCN)
  names(cnList) <- sbjVec

  xpndCNV <- function(sbj) {
    d <- as.data.table(cnList[[sbj]]@CNV.calls)
    if (nrow(d) == 0) {
      return(data.table(int[0],
                        type = character(),
                        nexons = numeric(),
                        BF = numeric(),
                        reads.expected = integer(),
                        reads.observed = integer(),
                        varID = character(),
                        subject = character()))
    }
    i1 <- d[ , unlist(mapply(seq, start.p, end.p, SIMPLIFY = FALSE))]
    d <- d[d[ , rep(.I, nexons)]]
    d <- d[ , .(type, nexons, BF, reads.expected, reads.observed, varID = id)]
    d <- cbind(int[i1], d)
    d[ , subject := sbj]
    d[]
  }

  calls <- lapply(sbjVec, xpndCNV)
  calls <- rbindlist(calls)

  makeCorTbl <- function(sbj) {
    tbl <- as.data.table(refList[[sbj]]$summary.stats)
    tbl[ , subject := sbj]
    tbl[ , selected := ref.samples %in% refList[[sbj]]$reference.choice]
    ref <- rowSums(cmat[ , refList[[sbj]]$reference.choice, drop = FALSE])
    v <- cor(cmat[ , sbj], ref)
    tbl[ , overallCor := v]
    tbl[ , overallPhi := unique(cnList[[sbj]]@phi)]
    tbl[]
  }

  correlations <- rbindlist(lapply(sbjVec, makeCorTbl))
  setcolorder(correlations, "subject")

  list(calls = calls[], correlations = correlations[])

}
