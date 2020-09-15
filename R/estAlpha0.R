#' @importFrom Rfast Digamma Trigamma colmeans rowmeans rowsums Lgamma Log
#' @import mcCNV
#' @export

estAlpha0 <- function(subjects, its = 5e3, tol = 1e-03, returnAllAlpha = FALSE) {

  data(subjectCounts, envir = environment())
  x <- cnvCountsToMatrix(subsetCounts(subjects))
  keep <- apply(x, MARGIN = 1, function(x) all(x >= 5 & x <= 2000))
  x <- x[keep, ]
  x <- t(sweep(x, MARGIN = 2, STATS = Rfast::colsums(x), FUN = "/"))
  ## code taken from Rfast::diri.nr2
  dm <- dim(x)
  n <- dm[1]
  p <- dm[2]
  m <- Rfast::colmeans(x)
  zx <- t(Log(x))
  down <- -sum(m * (Rfast::rowmeans(zx) - log(m)))
  sa <- 0.5 * (p - 1)/down
  a1 <- sa * m
  gm <- rowsums(zx)
  z <- n * Rfast::Digamma(sa)
  g <- z - n * Digamma(a1) + gm
  qk <- -n * Trigamma(a1)
  b <- sum(g/qk)/(1/z - sum(1/qk))
  a2 <- a1 - (g - b)/qk
  ## code added in
  i <- 1L
  tl <- numeric(length = its)
  a0 <- numeric(length = its)
  ll <- numeric(length = its)
  repeat {
    ## code taken from Rfast::diri.nr2
    a1 <- a2
    z <- n * digamma(sum(a1))
    g <- z - n * Digamma(a1) + gm
    qk <- -n * Trigamma(a1)
    b <- sum(g/qk)/(1/z - sum(1/qk))
    a2 <- a1 - (g - b)/qk
    ## code added in
    tl[i] <- sum(abs(a2 - a1))
    cat("iteration:", i, ", tol:", tl[i], "\n")
    a0[i] <- sum(a2)
    ll[i] <- n*Lgamma(a0[i]) - n*sum(Lgamma(a2)) + sum(zx*(a2 - 1))
    if (i == its || tl[i] < tol) break
    i <- i + 1L
  }
  names(a2) <- colnames(x)
  res <- list(N = length(a2),
              a0 = a0[i],
              a0vec = a0[1:i],
              loglik = ll[1:i],
              tol = tl[1:i])
  if (returnAllAlpha) res$a <- a2
  res

}
