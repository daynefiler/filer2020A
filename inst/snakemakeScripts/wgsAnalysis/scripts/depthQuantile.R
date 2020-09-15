#!/usr/bin/env Rscript

"Usage: depthQuantile.R <input> <output>

Options:
  -h --help     show this help page
  input         path to input samtools depth file
  output        path to output file
" -> doc

require(data.table, quietly = TRUE)
require(docopt, quietly = TRUE); opt <- docopt::docopt(doc)

dat <- fread(opt$input, select = 3, col.names = "depth")
qnt <- quantile(dat$depth, probs = seq(0, 1, 0.01))
write.table(qnt, opt$output, col.names = FALSE)
