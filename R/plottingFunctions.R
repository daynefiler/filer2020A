#' @title Plot procRes results
#' @description Plot procRes results
#' @param res simulation result summary, output from [procRes]
#' @param stat character, the statistic from [procRes] to plot
#' @param ylab character, text for ylab
#' @param h numeric, value to draw dashed horizontal line
#' @importFrom ggthemes tableau_seq_gradient_pal
#' @import graphics
#' @export

pltStat <- function(res, stat, ylab = toupper(stat), h = NULL) {
  mndat <- res$mnDat
  sddat <- res$sdDat
  par(mar = c(4, 4, 1, 1))
  plot.new()
  plot.window(ylim = c(0, 1), xlim = c(5, 100))
  abline(h = h, lty = "dashed", col = "darkgrey")
  abline(v = seq(5, 100, 5), lty = "dotted", col = "lightgrey")
  dep <- mndat[order(dep)][ , dep]
  sval <- mndat[order(dep)][ , get(stat)]
  SD   <- sddat[order(dep)][ , get(stat)]
  lines(dep, sval, col = "darkgray", lwd = 1.5)
  polygon(x = c(dep, rev(dep)),
          y = c(sval + 3*SD, rev(sval - 3*SD)),
          col = col2alpha("darkgray", alpha = 0.25),
          border = NA)
  C <- seq(40, 250, 30)
  # axis(side = 1, at = C*76002546/200E6, labels = paste0(C, "x"),
  #      tcl = -0.8, mgp = c(3, 2, 0), col = "gray50", col.axis = "gray50")
  # axis(side = 1, at = seq(10, 100, by = 10), tcl = -0.4)
  # mtext("Depth", side = 1, line = 4)
  # axis(side = 2)
  # mtext(ylab, side = 2, line = 3)
  axis(side = 1, at = seq(10, 100, by = 10))
  axis(side = 2)
  title(xlab = "Depth/sample (millions)", ylab = ylab)
}

#' @title Plot best prior across sequencing depth & variant width
#' @description Plot best prior across sequencing depth & variant width
#' @param dat data.table w/ prior, factor giving sequencing depth (dep), and
#' factor giving variant width (fwid)
#' @import graphics
#' @importFrom dlfUtils line2user
#' @export

pltBestPrior <- function(dat) {
  dat <- copy(dat)
  dat[ , dep := as.numeric(levels(dep))[dep]]
  bluePal <- tableau_gradient_pal()(seq(0, 1, length.out = 5))
  palette(colorRampPalette(bluePal)(length(levels(dat$fwid))))
  par(mar = c(4, 4, 1, 1), oma = c(0, 0, 2, 0))
  plot.new()
  plot.window(ylim = range(dat$prior),
              xlim = range(dat$dep), log = "y")
  abline(v = seq(5, 100, 5), lty = "dotted", col = "lightgrey")
  for (i in levels(dat$fwid)) {
    with(dat[fwid == i],
         points(dep, prior, col = i, lwd = 2, type = "o", pch = 16))
  }
  axis(side = 1, at = seq(10, 100, by = 10))
  axis(side = 2)
  mtext("Depth/sample (millions)", side = 1, line = 3)
  mtext("Prior", side = 2, line = 3)
  legend(x = mean(par('usr')[1:2]),
         y = line2user(1, 3, outer = TRUE),
         horiz = TRUE,
         pch = 16,
         lwd = 2,
         legend = unique(dat$fwid),
         bty = "n",
         xpd = NA,
         xjust = 0.5,
         col = unique(dat$fwid))
}

#' @title Plot prior by depth performance matrix
#' @description Plot prior by depth performance matrix
#' @param mndat object$mnDat, where object is output from [procRes]
#' @param wid variant width to plot
#' @param stat character, the statistic from [procRes] to plot
#' @importFrom ggthemes tableau_gradient_pal
#' @import graphics
#' @import plot.matrix
#' @export

pltPriorMatrix <- function(mnDat, wid = 1, stat = "mcc") {
  mat <- dcast(prior ~ dep, data = mnDat[fwid == wid], value.var = stat)
  setorder(mat, -prior)
  mat <- as.matrix(mat, rownames = "prior")
  par(mar = c(4, 4, 1, 4))
  bluePal <- tableau_gradient_pal()(seq(0, 1, length.out = 5))
  plot(mat,
       breaks = 40,
       xlab = "Depth/sample (millions)",
       ylab = "Prior",
       main = "",
       axis.col = NULL,
       axis.row = NULL,
       border = "white",
       fmt.key="%.2f",
       key = list(cex.axis = 0.6, tick = FALSE),
       col = colorRampPalette(rev(c(bluePal[5], "gray99"))))
  axis(side = 2,
       at = rev(seq_along(rownames(mat))),
       labels = rownames(mat),
       las = 2,
       tick = FALSE,
       cex.axis = 0.6)
  axis(side = 1,
       at = seq_along(colnames(mat)),
       labels = colnames(mat),
       tick = FALSE,
       cex.axis = 0.6)
}

#' @title Plot variance by mean contours
#' @description Plot variance by mean contours
#' @param dat subset of [realMeanVar]
#' @param grpVec factor, same length as \code{nrow(dat)}; defines groups for
#' plotting
#' @param nbin number of bins in each dimension to aggregate data
#' @param ncontour number of contours to plot
#' @param colVec vector of colors, same length as \code{levels(grpVec)}
#' @import graphics
#' @importFrom dlfUtils alpha2opaque log2dDen col2alpha
#' @export

pltMnVrCont <- function(dat,
                        grpVec,
                        colVec = palette()[unique(grpVec)],
                        nbin = 100,
                        ncontour = 20,
                        lgnd = TRUE) {
  dat <- copy(dat)
  if ('grpVec' %in% names(dat)) dat[ , grpVec := NULL]
  nLvl <- length(levels(grpVec))
  grpVec <- grpVec[!is.na(dat$intVar)]
  dat <- copy(dat)[!is.na(intVar)]
  q <- c(0.001, 0.999)
  mvr <- with(dat, c(quantile(intMean, q), quantile(intVar, q)))
  denLst <- vector(mode = "list", length = nLvl)
  colPal <- function(col) {
    aseq <- seq(0.5, 1, length.out = ncontour)
    sapply(col2alpha(col, aseq), alpha2opaque)
  }
  for (i in seq_along(levels(grpVec))) {
    denLst[[i]]$d <- with(dat[grpVec == levels(grpVec)[i]],
                          log2dDen(intMean, intVar, lims = mvr, nbin = nbin))
    denLst[[i]]$c <- colVec[i]
    denLst[[i]]$cp <- colPal(colVec[i])
    denLst[[i]]$l <- with(dat[grpVec == levels(grpVec)[i]],
                          lm(log10(intVar) ~ log10(intMean)))
    denLst[[i]]$m <- density(dat[grpVec == levels(grpVec)[i], log10(intMean)],
                             from = log10(mvr[1]), to = log10(mvr[2]))
    denLst[[i]]$v <- density(dat[grpVec == levels(grpVec)[i], log10(intVar)],
                             from = log10(mvr[3]), to = log10(mvr[4]))
  }
  addContour <- function(d) {
    contour(d$d,
            col = d$cp,
            add = TRUE,
            nlevels = ncontour,
            axes = FALSE,
            drawlabels = FALSE)
  }
  layout(matrix(c(2, 1, 4, 3), ncol = 2), heights = c(1, 5), widths = c(5, 1))
  par(mar = c(4, 4, 0.1, 0.1))
  plot.new(); plot.window(c(1, 100), c(1, 100))
  sapply(denLst, addContour)
  par(usr = log10(mvr))
  addLine <- function(d) abline(d$l, col = d$c, lwd = 2, lty = "dotted")
  sapply(denLst, addLine)
  xticks <- axisTicks(usr = log10(mvr[1:2]), log = TRUE)
  yticks <- axisTicks(usr = log10(mvr[3:4]), log = TRUE)
  axis(side = 1, at = log10(xticks), labels = xticks)
  axis(side = 2, at = log10(yticks), labels = yticks)
  title(xlab = "Mean per exon", ylab = "Variance per exon")
  par(mar = c(0.1, 4, 0.1, 0.1))
  abline(0, 1, lty = "dashed", col = "gray30", lwd = 2)
  plot.new()
  mnDenMax <- Reduce(max, lapply(denLst, function(x) x$m$y))
  par(usr = c(log10(mvr)[1:2], -0.1, mnDenMax*1.04))
  sapply(denLst, function(d) lines(d$m, col = d$c, lwd = 2))
  par(mar = c(4, 0.1, 0.1, 0.1))
  plot.new()
  vrDenMax <- Reduce(max, lapply(denLst, function(x) x$v$y))
  par(usr = c(-0.1, vrDenMax*1.04, log10(mvr)[3:4]))
  sapply(denLst, function(d) lines(x = d$v$y, y = d$v$x, col = d$c, lwd = 2))
  par(mar = rep(0, 4))
  plot.new()
  if (lgnd) {
    legend(x = "center",
           lwd = 4,
           col = colVec,
           legend = levels(grpVec),
           bty = "n",
           cex = 0.75)
  }
}

getQuantSub <- function(sub) {
  t1 <- !is.na(sub$vr)
  t2 <- sub[ , vr < quantile(vr, 0.999, na.rm = TRUE)]
  t3 <- sub[ , vr > quantile(vr, 0.001, na.rm = TRUE)]
  t4 <- sub[ , mn < quantile(mn, 0.999, na.rm = TRUE)]
  t5 <- sub[ , mn > quantile(mn, 0.001, na.rm = TRUE)]
  t1 & t2 & t3 & t4 & t5
}

#' @title Plot ECDF curves comparing two vectors
#' @description Plot ECDF curves comparing two vectors
#' @param v1 vector 1 values, gray solid line
#' @param v2 vector 2 values, black dotted line
#' @import graphics
#' @importFrom stats ecdf
#' @export

pltECDF <- function(v1, v2, q = 0.99, xlab) {
  v1ECDF <- ecdf(v1); v2ECDF <- ecdf(v2)
  xlim <- quantile(c(v1, v2), c(0, q), na.rm = TRUE)
  par(mar = c(4, 4, 1, 1))
  plot.new()
  plot.window(xlim, c(0, 1))
  abline(h = seq(0, 1, 0.1), lwd = 0.5, col = "gray90")
  xvals <- seq(xlim[1], xlim[2], length.out = 500)
  lines(xvals, v1ECDF(xvals), col = "gray70", lwd = 2)
  lines(xvals, v2ECDF(xvals), lty = "dotdash", lwd = 2)
  axis(1); axis(2)
  title(ylab = "Proportion of observations", xlab = xlab)
  ## Compute KS-test statistic
  mxPt <- which.max(abs(v1ECDF(xvals) - v2ECDF(xvals)))
  abline(v = xvals[mxPt], lty = "dashed", col = "darkgray")
  tst <- c(max(abs(v1ECDF(xvals) - v2ECDF(xvals))), 1.358*sqrt(1000/500^2))
  tst <- signif(tst, 3)
  text(x = grconvertX(0.95, "npc"), y = grconvertY(0.05, "npc"),
       paste0("test stat: ", tst[1], "\ncrit val: ", tst[2]),
       adj = c(1, 0))
}

#' @title Plot number of false calls by depth
#' @description Plot number of false calls by depth
#' @param res data.table w/ 'ACT', 'anyCNV', 'PRO', 'dep', & 'N' columns; e.g.
#' [sRes$noVariants$mcCNV]
#' @import beeswarm
#' @importFrom dlfUtils col2alpha
#' @export

pltFalseCalls <- function(res) {
  sub <- res[!ACT &  anyCNV & !PRO]
  sub[ , col := "darkgray"]
  sub[as.logical(dep %% 10), col := "darkorange"]
  dev.new()
  bsplt <- beeswarm(N ~ dep, data = sub,
                    log = TRUE,
                    do.plot = TRUE,
                    method = "center",
                    corral = "gutter",
                    spacing = 0.3,
                    cex = 0.5,
                    pch = 16)
  dev.off()
  par(mar = c(4, 4, 1, 1) + 0.1)
  plot.new()
  plot.window(xlim = range(bsplt$x), ylim = range(bsplt$y), log = 'y')
  abline(v = seq(0.5, 20.5, 1), lty = "dashed", col = "gray80")
  abline(h = 10^(0:5), lwd = 1.25, col = "gray80")
  abline(h = 5*10^(0:5), col = "gray80", lwd = 0.75)
  points(y ~ x, data = bsplt,
         cex = 0.5,
         col = col2alpha(sub$col, 1),
         pch = 16)
  axis(at = c(1, 5)*10^rep(0:5, each = 2),
       side = 2, tick = FALSE, las = 2,
       line = -1,
       labels = expression(1, 5, 10, 50, 100, 500, 1%*%10^3,
                           5%*%10^3, 1%*%10^4, 5%*%10^4, 1%*%10^5, 5))
  axis(side = 1, at = 1:20, labels = seq(5, 100, 5), tick = FALSE)
  title(xlab = "Depth/sample (millions)", ylab = "False-positive calls")
  boxplot(N ~ dep, data = sub,
          log = "y",
          add = TRUE,
          axes = FALSE,
          ann = FALSE,
          col = NA,
          border = "gray30",
          outline = FALSE,
          boxwex = 0.5,
          lwd = 0.7)
}

#' @title Plot summary statistics comparing two fits
#' @description Plot summary statistics comparing two fits
#' @param xMn x-axis data; object$mnDat, where object is output from [procRes]
#' @param yMn y-axis data; object$mnDat, where object is output from [procRes]
#' @param stat character, the statistic from [procRes] to plot
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param lbl character, the label for the upper left corner
#' @importFrom ggthemes tableau_gradient_pal
#' @import graphics
#' @export

pltStatCompare <- function(xRes, yRes, stat, xlab = "X", ylab = "Y",
                           lbl = toupper(stat)) {
  rng <- range(c(xRes$mnDat[ , get(stat)], yRes$mnDat[ , get(stat)]))
  par(mar = c(4, 4, 1, 1) + 0.1)
  plot.new()
  plot.window(rng, rng)
  abline(0, 1, lty = "dashed")
  points(x = xRes$mnDat[order(dep), get(stat)],
         y = yRes$mnDat[order(dep), get(stat)],
         col = "darkgray",
         type = "c",
         pch = 16)
  text(x = xRes$mnDat[order(dep), get(stat)],
       y = yRes$mnDat[order(dep), get(stat)],
       labels = xRes$mnDat[order(dep), dep],
       col = "darkgray",
       cex = 0.75,
       font = 2)
  axis(side = 1)
  axis(side = 2)
  text(grconvertX(0.05, "npc"), grconvertY(0.95, "npc"), adj = c(0, 1), lbl)
  title(xlab = xlab, ylab = ylab)
}

#' @title Plot subject stat by pool
#' @description Plot subject stat, e.g. totalMolCount, by pool
#' @param stat character giving the stat to plot, e.g. totalMolCount
#' @param poolVec character vector giving the pools to plot
#' @param logY logical, use log y-axis when TRUE
#' @param ylab (optional) alternate y-axis label; defaults to stat
#' @import data.table
#' @import graphics
#' @import beeswarm
#' @importFrom dlfUtils col2alpha
#' @export

pltSubjectStatByPool <- function(stat,
                                 poolVec = NULL,
                                 logY = TRUE,
                                 ylab = stat) {
  data(subjectCorr, envir = environment())
  data(subjectMeta, envir = environment())
  setkey(subjectCorr, subject)
  setkey(subjectMeta, subject)
  subjectCorr <- subjectMeta[subjectCorr]
  keep <- c("subject", "pool", "capture", "multiplexCapture", stat)
  if (is.null(poolVec)) poolVec <- unique(subjectMeta$pool)
  rd <- subjectCorr[pool %in% poolVec, unique(.SD), .SDcols = keep]
  setorder(rd, capture, multiplexCapture, pool)
  rd[ , GRP := .GRP, by = pool]
  dev.new()
  bs <- beeswarm(get(stat) ~ GRP,
                 cex = 0.8,
                 corral = "wrap",
                 log = logY,
                 data = rd)
  dev.off()
  bs <- as.data.table(bs)
  bs <- bs[ , .(x, y, GRP = as.integer(x.orig), stat = y.orig)]
  setkey(bs, GRP)
  setkey(rd, GRP)
  bs <- unique(rd[ , .(GRP, pool, capture, multiplexCapture)])[bs]
  par(mar = c(1, 4, 1, 1) + 0.1)
  plot.new()
  plot.window(xlim = range(bs$x),
              ylim = range(bs$y),
              log = ifelse(logY, "y", ""))
  points(y ~ x,
         data = bs,
         pch = ifelse(multiplexCapture, 16, 17),
         col = col2alpha('gray40'),
         cex = 0.8)
  boxplot(y ~ GRP,
          data = bs,
          frame = FALSE,
          add = TRUE,
          ann = FALSE,
          axes = FALSE,
          col = "transparent",
          border = "gray30",
          outline = FALSE,
          boxwex = 0.4)
  text(x = rd[ , as.integer(unique(GRP))] - 0.2,
       y = grconvertY(0.05, from = "nfc"),
       labels = rd[ , unique(pool)],
       srt = 90,
       adj = c(0, 0),
       cex = 0.75,
       xpd = NA)
  axis(side = 2)
  title(ylab = ylab)
  cap <- bs[ , .(x = length(unique(GRP))), by = capture]
  cap[ , x := cumsum(x) + 0.5]
  abline(v = cap$x, col = "gray20", lty = "dotted")
  text(x = cap$x - 0.2,
       y = grconvertY(0.95, "nfc"),
       labels = cap$capture,
       srt = 90,
       adj = c(1, 0),
       cex = 0.75,
       xpd = NA)
}

#' @title Plot alpha0 estimates by pool
#' @description Plot alpha0 estimates by pool
#' @param a0tbl data.table containing alpha0 values and other summary info
#' @param mcCol color for multiplexed capture
#' @param icCol color for independent capture
#' @import data.table
#' @import graphics
#' @import dlfUtils
#' @export

pltAlpha0 <- function(a0tbl, mcCol = "darkorange", icCol = "darkblue") {
  par(mar = c(4, 4.5, 2, 1) + 0.1)
  plot.new()
  plot.window(xlim = range(a0tbl[ , c(mnCount, mxCount)]),
              ylim = range(c(10, a0tbl$aMn)),
              log = "y")
  a0tbl[ ,
        lines(x = c(mnCount, mxCount),
              y = c(aMn, aMn),
              lwd = 1.5,
              col = ifelse(mc, mcCol, icCol)),
        by = pool]
  points(aMn ~ mdCount, data = a0tbl,
         pch = ifelse(idt, 15, 17),
         col = ifelse(mc, mcCol, icCol))
  axis(side = 1)
  axis(side = 2)
  title(ylab = expression(hat(alpha)[0]/N),
        xlab = "Total counts (median/range)")
  legend(x = mean(par('usr')[1:2]), y = line2user(2, 3),
         xjust = 0.5,
         legend = c("AGL", "IDT", "MC", "IC"),
         lwd = c(NA, NA, 3, 3),
         col = c("black", "black", mcCol, icCol),
         pch = c(17, 15, NA, NA),
         horiz = TRUE,
         bty = "n",
         xpd = NA)
}

