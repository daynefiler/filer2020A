#' @title Calculate log-normal fit of capture probabilities
#' @description Calculate log-normal fit of capture probabilities
#' @param p character of length 1, pool name
#' @param plt logical of length 1, draw the histogram w/ fit when TRUE
#' @import mcCNV
#' @import data.table
#' @importFrom MASS fitdistr
#' @import graphics
#' @export

fitLogNorm <- function(p, plt = TRUE) {
  data("subjectMeta", envir = environment())
  cts <- subsetCounts(subjectMeta[pool == p, subject])
  int <- cts[ , .(medDep = median(molCount)), by = .(seqnames, start, end)]
  int[ , captureProb := medDep/sum(medDep)]
  lnFit <- fitdistr(int[medDep > 0, captureProb], 'lognormal')
  if (plt) {
    hist(int$captureProb,
         breaks = 1000,
         xlim = quantile(int$captureProb, c(0, 0.99)),
         freq = FALSE,
         border = "darkgrey",
         col = "lightgrey",
         main = "",
         axes = FALSE,
         ann = FALSE)
    axis(side = 1); title(xlab = "Probability of capture")
    xvls <- seq(par('usr')[1], par('usr')[2], length.out = 500)
    lines(x = xvls,
          y = dlnorm(x = xvls,
                     meanlog = lnFit$estimate['meanlog'],
                     sdlog = lnFit$estimate['sdlog']),
          lwd = 2)
  }
  lnFit$estimate
}
