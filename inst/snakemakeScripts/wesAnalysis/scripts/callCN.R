#!/usr/bin/env Rscript

"Usage: callCN.R [options] [--width=<width>]... <counts>...

Options:
  -h --help            show this help page
  --prior=<prior>      prior value used in calling copy number [default: 0.2]
  --outfile=<outfile>  the output file path & name [default: pool.counts]
  --width=<width>      calling window width (number of contiguous exons to 
                       aggregate); can provide multiple values -- must be 
                       pos non-zero integers [default: 1]
  counts               path to .RDS count objects
" -> doc

require(mcCNV, quietly = TRUE)
require(docopt, quietly = TRUE); opt <- docopt::docopt(doc)

wid <- as.integer(opt$width)
if (any(is.na(wid)) || any(wid < 1)) {
  stop("Invalid width given, must all be positive, non-zero integers")
}
cat("Gathering counts... ")
cts <- cnvGatherCounts(opt$counts)
cat("done.\n")
res <- cnvCallCN(counts = cts, 
                 prior = as.numeric(opt$prior),
                 width = wid, 
                 delta = 20L, 
                 iterations = 30L, 
                 verbose = TRUE)
cat("Writing outputs... ")
saveRDS(res, file = opt$outfile)
cat("done.\n")
