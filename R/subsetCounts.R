#' @title Make counts object from subjectCounts
#' @description Make [counts object][mcCNV::validObjects] from [subjectCounts]
#' @param subjects vector of subjects from [subjectCounts]
#' @return mcCNV counts object for the given subjects
#' @import data.table
#' @export

subsetCounts <- function(subjects) {
  data("subjectCounts", envir = environment())
  data("subjectMeta", envir = environment())
  cap <- unique(subjectMeta[subject %in% subjects, capture])
  cols <- c("seqnames", "start", "end")
  sub <- subjectCounts[capture == cap, .SD, .SDcols = c(cols, subjects)]
  sub <- melt(sub, id.vars = cols,
              measure.vars = subjects,
              variable.name = "subject",
              value.name = "molCount")
  sub[ , nCoverMult := NA_integer_]
  sub[]
}
