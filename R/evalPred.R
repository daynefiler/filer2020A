#' @title Evaluate predictions
#' @description Evaluate model predictions
#' @param x logical, prediction of TRUE (positive) or FALSE (negative)
#' @param y logical, the "truth" set to compare x to
#' @importFrom dlfUtils calcMCC
#' @importFrom caret confusionMatrix
#' @export

evalPred <- function(x, y) {
  stopifnot(is.logical(x))
  stopifnot(is.logical(y))
  stopifnot(length(x) == length(y))
  tb <- table(x, y)
  mcc <- calcMCC(tp = tb["TRUE",  "TRUE"],
                 tn = tb["FALSE", "FALSE"],
                 fp = tb["TRUE",  "FALSE"],
                 fn = tb["FALSE", "TRUE"])
  fdr <- tb["TRUE",  "FALSE"]/(tb["TRUE",  "TRUE"] + tb["TRUE",  "FALSE"])
  cm <- confusionMatrix(tb, positive = "TRUE")
  out <- c(MCC = mcc,
           TPR = unname(cm$byClass["Sensitivity"]),
           FDR = fdr,
           PPV = unname(cm$byClass["Pos Pred Value"]),
           BalAcc = unname(cm$byClass["Balanced Accuracy"]))
  signif(out, 3)
}



