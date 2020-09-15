#!/usr/bin/env Rscript

"Usage: getCounts.R <bamfile> <intfile> <outfile> 

Options:
  -h --help   show this help page
  bamfile     the path to the bamfile for processing
  intfile     the path to the interval object
  outfile     the output file path & name
" -> doc

require(mcCNV, quietly = TRUE)
require(docopt, quietly = TRUE); opt <- docopt::docopt(doc)

cts <- cnvGetCounts(bamfile = opt$bamfile, 
                    interval = readRDS(opt$intfile), 
                    verbose = TRUE)
cat("Writing output file...")
saveRDS(cts, file = opt$outfile)
cat("done.\n")
