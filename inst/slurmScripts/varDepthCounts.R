##----------------------------------------------------------------------------##
## Script to create the varDepth simulated counts
##----------------------------------------------------------------------------##

library(mcCNV)
library(rslurm)
library(data.table)
library(filer2020A)

## Directory for storing the data
wd <- getwd()

data(subjectMeta)
pl1Int <- subsetCounts(subjectMeta[pool == "Pool1", subject])
pl1Int <- pl1Int[ , .(ttl = sum(molCount)), by = .(seqnames, start, end)]
pl1Int[ , captureProb := ttl/sum(ttl)]
pl1Int[ , ttl := NULL]

deps <- as.integer(seq(5, 100, 5)) ## Sequencing depths

## Set up the file system
odir <- file.path(wd, "varDepthCounts")
if (dir.exists(odir)) unlink(odir, recursive = TRUE, force = TRUE)
dir.create(odir)

pars <- data.table(expand.grid(dep = deps, rep = seq(200)))
set.seed(1234)
pars[ , seed := sample(1e6, .N)]
pars[ , odir := odir]

simPool <- function(dep, rep, seed, odir) {
  wndw <- as.integer(c(dep - 0.3*dep, dep + 0.3*dep)*1e6)
  cnt <- try(cnvSimPool(nSubjects = 16L,
                        countRange = wndw,
                        interval = pl1Int,
                        seed = seed,
                        variantWidth = 1L))
  fname <- sprintf("varDepth_d%0.3d_r%0.4d.counts", dep, rep)
  saveRDS(cnt, file = file.path(odir, fname))
  !is(cnt, 'try-error')
}

slurm_apply(f = simPool,
            add_objects = "pl1Int",
            params = pars,
            nodes = nrow(pars),
            cpus_per_node = 1,
            jobname = "varDepthCounts",
            slurm_options = list(mem = 8000,
                                 array = sprintf("0-%d%%%d",
                                                 nrow(pars) - 1,
                                                 1000),
                                 'cpus-per-task' = 1,
                                 error =  "%A_%a.err",
                                 output = "%A_%a.out",
                                 time = "2-00:00:00"))


