##----------------------------------------------------------------------------##
## Script to perform mcCNV algorithm on varDepth simulation
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
odir <- file.path(wd, "varDepthMcCalls")
if (dir.exists(odir)) unlink(odir, recursive = TRUE, force = TRUE)
dir.create(odir)

pars[ , odir := odir]

doCalc <- function(fl, dep, odir) {
  dat <- readRDS(fl)
  ## Setup output
  ifl <- basename(fl)
  rep <- as.integer(sub("r", "", str_extract(ifl, "r[0-9]{4}")))
  ofmt <- sub(".counts$", ".mcCallSummary", ifl)
  mcFl <- file.path(odir, ofmt)
  ## mcCNV calls
  mcCalls <- try(cnvCallCN(counts = dat, verbose = TRUE))
  mcFail <- is(mcCalls, 'try-error')
  if (!mcFail) {
    mcRes <- mcCalls[ , .N, by = .(actCN, CN, passFilter)]
    mcRes[ , dep := dep]
    mcRes[ , rep := rep]
    saveRDS(mcRes, file = mcFl)
  }
  !mcFail
}

slurm_apply(f = doCalc,
            params = pars,
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "varDepthMcCalls",
            slurm_options = list(mem = 12000,
                                 array = sprintf("0-%d%%%d",
                                                 nrow(pars) - 1,
                                                 1000),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "10-00:00:00"))




