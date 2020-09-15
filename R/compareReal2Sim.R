#' @title Calculate interval statistics for real & simulated pool based on real
#' @description Calculate interval statistics for real & simulated pool based
#' on the interval created from the real counts
#' @param p character of length 1, pool name
#' @param seed seed value passed to cnvSimPool
#' @import mcCNV
#' @import data.table
#' @export

compareReal2Sim <- function(p, seed = 1234) {
  data("subjectMeta", envir = environment())
  cts <- subsetCounts(subjectMeta[pool == p, subject])
  int <- cts[ , .(medDep = median(molCount)), by = .(seqnames, start, end)]
  int[ , captureProb := medDep/sum(medDep)]
  rng <- subjectMeta[pool == p, range(totalMolCount)]
  sim <- cnvSimPool(countRange = rng, interval = int, seed = seed)
  stats <- list()
  stats$real <- calcIntStats(cts)
  stats$sim  <- calcIntStats(sim)
  stats
}
