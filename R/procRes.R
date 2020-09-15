#' @title Process call summaries
#' @description Takes simulation study results, and calculates confusion matrix
#' statistics based.
#' @param x results of simulation study
#' @import data.table
#' @importFrom dlfUtils calcMCC
#' @export

procRes <- function(x) {
  if (!"dep"   %in% names(x)) x[ , dep := NA]
  x[ , N := as.numeric(N)]
  xSmry <- x[ ,
              .(tp = sum(N[ actVar &  prdVar], na.rm = TRUE),
                fp = sum(N[!actVar &  prdVar], na.rm = TRUE),
                tn = sum(N[!actVar & !prdVar], na.rm = TRUE),
                fn = sum(N[ actVar & !prdVar], na.rm = TRUE)),
              by = .(dep, rep)]
  xSmry[ , fpr := fp/(fp + tn)]
  xSmry[ , tpr := tp/(tp + fn)]
  xSmry[ , spc := tn/(fp + tn)]
  xSmry[ , fdr := fp/(fp + tp)]
  xSmry[ , ba  := (tpr + spc)/2]
  xSmry[ , ppv := tp/(tp + fp)]
  xSmry[ , npv := tn/(tn + fn)]
  xSmry[ , dor := (tp/fp)/(fn/tn)]
  xSmry[ , mcc := calcMCC(tp, tn, fp, fn)]
  xMn <- xSmry[ , lapply(.SD, mean), by = dep]
  xSD <- xSmry[ , lapply(.SD, sd),   by = dep]
  list(mnDat = xMn, sdDat = xSD, all = xSmry)
}
