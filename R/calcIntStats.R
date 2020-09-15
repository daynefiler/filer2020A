#' @title Calculate mean, var, size factor, phi, & shrunken phi
#' @description Calculate mean, var, size factor, phi, & shrunken phi
#' @param counts [counts object][mcCNV::validObjects]
#' @import mcCNV
#' @import data.table
#' @export

calcIntStats <- function(counts) {
  stopifnot(cnvValidCounts(counts))
  counts <- copy(counts)
  counts[ , use := molCount > 10]
  counts[ , intName := sprintf("%s-%s:%s", seqnames, start, end)]
  setindex(counts, intName)
  setindex(counts, subject)
  counts[ , geomn := exp(mean(log(molCount[use]))), by = intName]
  counts[ , sf := median(molCount/geomn, na.rm = TRUE), by = subject]
  counts[ , mn := mean(molCount[use]/sf[use]), by = intName]
  counts[ , vr :=  var(molCount[use]/sf[use]), by = intName]
  counts[ , phi := (vr - mn)/(mn^2)]
  keep <- c("mn", "vr", "phi")
  int <- counts[ , .SD[1], by = .(seqnames, start, end), .SDcols = keep]
  setnames(int, keep, c("intMean", "intVar", "intPhi"))
  int[]
}
