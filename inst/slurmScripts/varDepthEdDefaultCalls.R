##----------------------------------------------------------------------------##
## Script to perform ExomeDepth algorithm w/ defaults on varDepth simulation
##----------------------------------------------------------------------------##

library(mcCNV)
library(rslurm)
library(data.table)
library(stringr)
library(filer2020A)

## Directory for storing the data
wd <- getwd()

pars <- data.table(fl = Sys.glob(file.path(wd, "varDepthCounts/*.counts")))
pars[ , dep := as.integer(sub("d", "", str_extract(fl, "d[0-9]{3}")))]
## Set up the file system
odir <- file.path(wd, "varDepthEdDefaultCalls")
if (dir.exists(odir)) unlink(odir, recursive = TRUE, force = TRUE)
dir.create(odir)

pars[ , odir := odir]

doCalc <- function(fl, dep, odir) {
  dat <- readRDS(fl)
  ## Setup output
  ifl <- basename(fl)
  rep <- as.integer(sub("r", "", str_extract(ifl, "r[0-9]{4}")))
  ofmt <- sub(".counts$", ".edDefaultCallSummary", ifl)
  edFl <- file.path(odir, ofmt)
  ## ExomeDepth calls
  edCalls <- try(runExomeDepth(counts = dat))
  edFail <- is(edCalls, 'try-error')
  if (!edFail) {
    setkey(edCalls$calls, subject, seqnames, start, end)
    setkey(dat, subject, seqnames, start, end)
    edRes <- edCalls$calls[dat][ , .N, by = .(actCN, type)]
    edRes[ , dep := dep]
    edRes[ , rep := rep]
    saveRDS(edRes, file = edFl)
  }
  !edFail
}

slurm_apply(f = doCalc,
            params = pars,
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "varDepthEdDefaultCalls",
            slurm_options = list(mem = 12000,
                                 array = sprintf("0-%d%%%d",
                                                 nrow(pars) - 1,
                                                 1000),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "10-00:00:00"))




