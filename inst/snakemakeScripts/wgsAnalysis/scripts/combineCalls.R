#!/usr/bin/env Rscript

"Usage: combineCalls.R <sbj> <pytor> <erds> <int> <out>

Options:
  -h --help     show this help page
  sbj           character giving subject name
  pytor         file path to cnvpytor 'sbj.calls' file
  erds          file path to erds 'sbj.events' file
  int           file path to mcCNV interval object
  out           file path for output file (RDS object)
" -> doc

require(data.table, quietly = TRUE)
require(docopt, quietly = TRUE); opt <- docopt::docopt(doc)

## Read & process cnvpytor calls
pCalls <- fread(opt$pytor, select = 1:2, col.names = c("cnvpytor", "reg"))
regNms <- c("seqnames", "start", "end")
pCalls[ , c(regNms) := tstrsplit(reg, ":|-", type.convert = TRUE)]
pCalls <- pCalls[seqnames %in% 1:22]
pCalls[ , reg := NULL]
old <- c("deletion", "duplication")
pCalls[ , cnvpytor := c("del", "dup")[match(cnvpytor, old)]]

## Read & process erds calls
eCalls <- fread(opt$erds, select = c(1:3, 5))
setnames(eCalls, c(regNms, "erds"))
eCalls <- eCalls[seqnames %in% 1:22]
eCalls[ , erds := tolower(erds)]

## Map calls to interval
dat <- readRDS(opt$int)
dat <- dat[ , .(seqnames, start, end)]
setkeyv(dat, regNms)

map2dat <- function(x) {
  setkeyv(x, regNms)
  dat <- foverlaps(dat, x, mult = "first", nomatch = NA)
  dat[ , c("start", "end") := NULL]
  setnames(dat, c("i.start", "i.end"), c("start", "end"))
  dat[]
}
dat <- map2dat(eCalls)
dat <- map2dat(pCalls)
dat[ , subject := opt$sbj]
dat <- dat[!is.na(cnvpytor) | !is.na(erds)]
saveRDS(dat, file = opt$out)

