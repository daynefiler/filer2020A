#!/usr/bin/env Rscript

"Usage: exomeDepth.R [options] <counts>...

Options:
  -h --help                show this help page
  --transProb=<transProb>  transition probaility [default: 0.0001]
  --cnvLength=<cnvLength>  expected variant length [default: 50000]
  --outfile=<outfile>      the output file path [default: pool.exomeDepth]
  counts                   path to .RDS count objects
" -> doc

require(mcCNV, quietly = TRUE)
require(filer2020A, quietly = TRUE)
require(docopt, quietly = TRUE); opt <- docopt::docopt(doc)

cat("Gathering counts... ")
cts <- cnvGatherCounts(opt$counts)
cat("done.\n")
res <- runExomeDepth(counts = cts, 
                     transProb = as.numeric(opt$transProb),
                     cnvLength = as.numeric(opt$cnvLength))
cat("Writing outputs... ")
saveRDS(res, file = opt$outfile)
cat("done.\n")
