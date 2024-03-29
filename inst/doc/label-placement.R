## ----setup, echo=FALSE, results='hide', warning=FALSE, message=FALSE------------------------------
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  background = '#F7F7F7',
  fig.align = 'center',
  dev = 'png',
  comment = "#>"
)
options(width = 100, stringsAsFactors = FALSE, timeout = 600)

# keep examples from using more than 2 cores
data.table::setDTthreads(Sys.getenv("OMP_THREAD_LIMIT", unset = 2))

library(aqp, quietly = TRUE)
.bluecolors <- hcl.colors(n = 25, palette = 'Blues')[-25]

## ----fig.width=6, fig.height=5, echo = FALSE------------------------------------------------------
data("osd")
par(mar = c(0, 0, 0, 0))
plotSPC(osd, name = NA, fixLabelCollisions = FALSE, cex.names = 1.75, hz.depths = TRUE, depth.axis = FALSE, cex.id = 0.66, name.style = 'center-center')

## ----fig.width=8, fig.height=5--------------------------------------------------------------------
library(aqp)

# soil profile template
x <- quickSPC(
  "SPC:Oi|Oe|AAAAA|E1|E2|Bhs1Bhs1Bhs1|Bhs2|CCCCCC|RRRRRRRRRR", 
  interval = 3
)

# set variability in horizon thickness
horizons(x)$thick.sd <- 3
# setup horizon designations
horizons(x)$nm <- factor(x$name, levels = c('Oi', 'Oe', 'A', 'E1', 'E2', 'Bhs1', 'Bhs2', 'C', 'R'))

# simulate 3 realizations of original profile
set.seed(10101)
s <- perturb(x, n = 3, thickness.attr = 'thick.sd', min.thickness = 2)

# create 3 copies of each simulation
s <- duplicate(s, times = 3)

## ----fig.width=8, fig.height=5, echo = FALSE------------------------------------------------------
.a <- list(cex.names = 0.8, color = 'nm', name = NA, show.legend = FALSE, hz.depths = TRUE, depth.axis = FALSE, cex.id = 0.66, col.palette = .bluecolors)

options(.aqp.plotSPC.args = .a)

par(mar = c(0, 0, 0, 0))
plotSPC(s, fixLabelCollisions = FALSE)
title('Example Data', line = -1.5)

## ----fig.width=8, fig.height=5, echo = FALSE------------------------------------------------------
par(mar = c(0, 0, 0, 0))
# select "Electrostatic" label placement method
# adjust default charge density to 0.5
plotSPC(s, fixOverlapArgs = list(method = 'E', q = 0.5))
title('Label Adjustment by Electrostatic Simulation', line = -1.5)

## ----fig.width=8, fig.height=5, echo = FALSE------------------------------------------------------
par(mar = c(0, 0, 0, 0))
plotSPC(s, fixOverlapArgs = list(method = 'S'))
title('Label Adjustment by Simulated Annealing', line = -1.5)

options(.aqp.plotSPC.args = NULL)

## ----fig.width=8, fig.height=5--------------------------------------------------------------------
x <- quickSPC(
  "SPC:Oi|Oe|AAA|E1|E2|E3|BhsBhsBhsBhs|Bt1|Bt2|Bt3Bt3|CCCCCC|Ab1|Ab2|2C2C2C2C2C2C|2Cr|2R2R2R2R2R2R2R2R", 
  interval = 1
)

x$z <- as.numeric(x$hzID)

## ----fig.width=8, fig.height=5--------------------------------------------------------------------
# pretty colors
.bluecolors <- hcl.colors(n = 25, palette = 'Blues')[-25]

# plotSPC arguments
.a <- list(
  width = 0.2, 
  max.depth = 40,
  hz.depths = TRUE, 
  name.style = 'center-center', 
  cex.names = 1.5,
  name = NA,
  depth.axis = FALSE, 
  color = 'z',
  show.legend = FALSE,
  print.id = FALSE,
  col.palette = .bluecolors
)

# set plotSPC default arguments
options(.aqp.plotSPC.args = .a)

# wrapper function to test label collision solutions
testIt <- function(x, ...) {
  
  # make sketches
  plotSPC(x, ...)
  
  # a normalized index of label adjustment
  .LAI <- get('last_spc_plot', envir = aqp.env)$hz.depth.LAI
  .LAI <- ifelse(is.na(.LAI), 0, .LAI)
  
  # annotate with label adjustment index
  .txt <- sprintf("LAI: %0.3f", .LAI)
  mtext(.txt, side = 1, at = 1, line = -1.5, cex = 0.8)
}

## ----fig.width=8, fig.height=5--------------------------------------------------------------------
par(mar = c(1, 0, 0, 0), mfcol = c(1, 8))

testIt(x, fixLabelCollisions = FALSE)
title('No\nAdjustment', line = -3.5, adj = 0.5)

testIt(x, fixOverlapArgs = list(method = 'S'))
title('SANN\nsim 1', line = -3.5, adj = 0.5)

testIt(x, fixOverlapArgs = list(method = 'S'))
title('SANN\nsim 2', line = -3.5, adj = 0.5)

testIt(x, fixOverlapArgs = list(method = 'S'))
title('SANN\nsim 3', line = -3.5, adj = 0.5)

testIt(x, fixOverlapArgs = list(method = 'E', q = 1.5))
title('Electrostatic\nq = 1.5', line = -3.5, adj = 0.5)

testIt(x, fixOverlapArgs = list(method = 'E', q = 1))
title('Electrostatic\nq = 1', line = -3.5, adj = 0.5)

testIt(x, fixOverlapArgs = list(method = 'E', q = 0.5))
title('Electrostatic\nq = 0.5', line = -3.5, adj = 0.5)

testIt(x, fixOverlapArgs = list(method = 'E', q = 0.25))
title('Electrostatic\nq = 0.25', line = -3.5, adj = 0.5)

## ----fig.width=8, fig.height=5--------------------------------------------------------------------
library(aqp)
library(soilDB)

s <- c('inks' , 'pardee', 'clarksville', 'palau', 'hao', 'inks', 'eheuiki', 'puaulu', 'zook', 'cecil')
x <- fetchOSD(s)

par(mar = c(0, 0, 0, 3))

.args <- list(width = 0.3, name.style = 'center-center', hz.depths = TRUE, cex.names = 1)
options(.aqp.plotSPC.args = .args)

plotSPC(x, fixOverlapArgs = list(method = 'E', q = 1.25), max.depth = 151)
plotSPC(x, fixOverlapArgs = list(method = 'E', q = 1), max.depth = 151)
plotSPC(x, fixOverlapArgs = list(method = 'E', q = 0.75), max.depth = 151)
plotSPC(x, fixOverlapArgs = list(method = 'E', q = 0.5), max.depth = 151)
plotSPC(x, fixOverlapArgs = list(method = 'E', q = 0.25), max.depth = 151)

plotSPC(x, fixOverlapArgs = list(method = 'S'), max.depth = 151)

## -------------------------------------------------------------------------------------------------
x <- c(1, 2, 3, 3.4, 3.5, 5, 6, 10)

overlapMetrics(x, thresh = 0.5)

## ----message=TRUE---------------------------------------------------------------------------------
# vector of positions, typically labels but could be profile sketch alignment on the x-axis
s <- c(1, 2, 2.3, 4, 5, 5.5, 7)

# simulated annealing, solution is non-deterministic
fixOverlap(s, thresh = 0.5, method = 'S')

# electrostatics-inspired simulation of particles
# solution is deterministic
fixOverlap(s, thresh = 0.5, method = 'E')

## ----fig.width=8.5, fig.height=5------------------------------------------------------------------
evalMethods <- function(x, thresh, q, ...) {
  
  cols <- hcl.colors(n = 9, palette = 'Zissou 1', rev = TRUE)
  cols <- colorRampPalette(cols)(length(x))
  
  z <- fixOverlap(x, thresh = thresh, method = 'E', maxIter = 100, trace = TRUE, q = q)
  .n <- nrow(z$states)
  
  op <- par(mar = c(0, 2, 2, 0.5), bg = 'black', fg = 'white')
  layout(matrix(c(1, 2, 3, 4), ncol = 2, nrow = 2), heights = c(0.33, 0.66))
  
  plot(seq_along(z$cost), z$cost, las = 1, type = 'b', axes = FALSE, cex = 0.66, xlim = c(1, .n))
  mtext(text = sprintf("Converged (%s): %s", .n, z$converged), at = 0, side = 3, line = 0, cex = 0.75, font = 3, adj = 0)
  matplot(rbind(x, z$states), type = 'l', lty = 1, las = 1, axes = FALSE, col = cols, lwd = 1)
  
  points(x = rep(1, times = length(x)), y = x, cex = 0.66, pch = 16, col = cols)
  points(x = rep(.n + 1, times = length(x)), y = z$x, cex = 0.66, pch = 16, col = cols)
  
  text(x = 1, y = x, col = cols, labels = seq_along(x), cex = 0.66, font = 2, pos = 2)
  text(x = .n + 1, y = z$x, col = cols, labels = seq_along(x), cex = 0.66, font = 2, pos = 4)
  
  axis(side = 2, at = unique(x), labels = round(unique(x), 1), col.axis = par('fg'), las = 1, cex.axis = 0.6)
  title(main = 'Electrostatic Simulation', line = 1, col.main = 'white')
  
  ## SANN_1D doesn't always preserve rank ordering
  ##  ->> not designed to use unsorted input
  ##  ->> maybe impossible with ties in x?
  
  z <- fixOverlap(x, thresh = thresh, method = 'S', trace = TRUE, maxIter = 1000)
  .n <- nrow(z$states)
  
  plot(seq_along(z$stats), z$stats, las = 1, type = 'b', axes = FALSE, cex = 0.66, xlim = c(1, .n))
  mtext(text = sprintf("Converged (%s): %s", .n, z$converged), at = 0, side = 3, line = 0, cex = 0.75, font = 3, adj = 0)
  
  matplot(z$states, type = 'l', lty = 1, las = 1, axes = FALSE, col = cols)
  
  points(x = rep(1, times = length(x)), y = z$states[1, ], cex = 0.66, pch = 16, col = cols)
  points(x = rep(.n, times = length(x)), y = z$x, cex = 0.66, pch = 16, col = cols)
  
  text(x = 1, y = z$states[1, ], col = cols, labels = seq_along(x), cex = 0.66, font = 2, pos = 2)
  text(x = .n, y = z$x, col = cols, labels = seq_along(x), cex = 0.66, font = 2, pos = 4)
  
  axis(side = 2, at = unique(x), labels = round(unique(x), 1), col.axis = par('fg'), las = 1, cex.axis = 0.6)
  
  title(main = 'Simulated Annealing', line = 1, col.main = 'white')
  
  # reset graphics state
  par(op)
  layout(1)
  
}

## ----fig.width=8.5, fig.height=5------------------------------------------------------------------
# explore effect of charge (q)
# too large -> chaos
x <- c(0, 2, 5, 12, 18, 20, 35, 40, 50, 56, 90, 120, 145, 150)

# just about right, very few perturbations required
evalMethods(x, thresh = 5, q = 1.1)

# ok, but now most label positions are affected
evalMethods(x, thresh = 5, q = 1.8)

# too high, wasting time on more iterations
evalMethods(x, thresh = 5, q = 3)

# far too high, wasting more time with little gain
evalMethods(x, thresh = 5, q = 4)

# chaos and failure to converge
evalMethods(x, thresh = 5, q = 5)

## ----eval = FALSE---------------------------------------------------------------------------------
#  # threshold too large
#  evalMethods(x, thresh = 10, q = 3)
#  
#  
#  # large threshold
#  x <- c(0, 5, 12, 18, 20, 35, 40, 55, 90, 120, 145, 150)
#  evalMethods(x, thresh = 9, q = 2)
#  
#  # single iteration enough
#  x <- c(0, 3, 20, 35, 40, 55, 90, 120, 145, 150)
#  evalMethods(x, thresh = 6, q = 1)
#  
#  # clusters
#  x <- sort(c(0, jitter(rep(10, 3)), jitter(rep(25, 3)), jitter(rep(90, 3)), 150))
#  evalMethods(x, thresh = 6, q = 3)
#  evalMethods(x, thresh = 6, q = 2)
#  
#  
#  ## impact of scale / offset
#  x <- c(0, 5, 12, 18, 20, 35, 40, 50, 120, 145, 150)
#  
#  # works as expected
#  evalMethods(x, thresh = 5, q = 1.1)
#  
#  # works as expected, as long as threshold is scaled
#  evalMethods(x / 10, thresh = 5 / 10, q = 1.1)
#  
#  # works as expected, as long as threshold is scaled
#  evalMethods(x * 10, thresh = 5 * 10, q = 1.1)
#  
#  
#  # all work as expected, threshold not modified
#  evalMethods(x + 10, thresh = 5, q = 1.1)
#  evalMethods(x + 100, thresh = 5, q = 1.1)
#  evalMethods(x + 1000, thresh = 5, q = 1.1)
#  
#  # works as expected
#  x <- c(315, 325, 341, 353, 366, 374, 422)
#  fixOverlap(x, thresh = 9.7, q = 1, method = 'E')
#  evalMethods(x, thresh = 9.7, q = 1)
#  
#  
#  x <- c(1.0075, 1.1200, 1.3450, 1.6450, 1.8700, 1.8825)
#  fixOverlap(x, thresh = 0.05442329, q = 1)
#  evalMethods(x, thresh = 0.05442329, q = 1)

## ----fig.width=8, fig.height=6--------------------------------------------------------------------
tracePlot <- function(x, z, cex.axis.labels = 0.85) {
  # setup plot device
  op <- par(mar = c(4, 4, 1, 1), bg = 'black', fg = 'white')
  layout(matrix(c(1,2,3)), widths = 1, heights = c(1,1,2))
  
  # order:
  # B: boundary condition violation
  # O: rank (order) violation
  # +: accepted perturbation
  # -: rejected perturbation
  cols <- c(grey(0.5), grey(0.85), 'royalblue', 'firebrick')
  
  cols.lines <- hcl.colors(9, 'Zissou 1', rev = TRUE)
  cols.lines <- colorRampPalette(cols.lines)(length(x))
  
  # total overlap (objective function) progress
  plot(
    seq_along(z$stats), z$stats, 
    type = 'h', las = 1,
    xlab = 'Iteration', ylab = 'Total Overlap',
    axes = FALSE,
    col = cols[as.numeric(z$log)]
  )
  
  axis(side = 2, cex.axis = cex.axis.labels, col.axis = 'white', las = 1, line = -2)
  mtext('Overlap', side = 2, line = 2, cex = cex.axis.labels, font = 2)
  
  # deviation from original configuration
  plot(
    seq_along(z$stats), z$ssd, 
    type = 'h', las = 1,
    xlab = 'Iteration', ylab = 'Deviation',
    axes = FALSE,
    col = cols[as.numeric(z$log)]
  )
  
  axis(side = 2, cex.axis = cex.axis.labels, col.axis = 'white', las = 1, line = -2)
  mtext('Deviation', side = 2, line = 2, cex = cex.axis.labels, font = 2)
  legend('top', legend = c('boundary\nviolation', 'rank\nviolation', 'accepted\nperturbation', 'rejected\nperturbation'), col = cols, bty = 'n', horiz = TRUE, inset = -0.5, lty = 1, lwd = 2, xpd = NA)
  
  # adjustments at each iteration
  matplot(
    z$states, type = 'l', 
    lty = 1, las = 1, 
    xlab = 'Iteration', ylab = 'x-position',
    axes = FALSE,
    col = cols.lines
  )
  
  axis(side = 2, cex.axis = cex.axis.labels, col.axis = 'white', las = 1, at = x, labels = round(x, 1))
  axis(side = 4, cex.axis = cex.axis.labels, col.axis = 'white', las = 1, at = z$x, labels = round(z$x, 1), line = -2)
  mtext('Position', side = 2, line = 2.5, cex = cex.axis.labels, font = 2)
  
  axis(side = 1, cex.axis = 1, col.axis = 'white', line = 0)
  mtext('Iteration', side = 1, line = 2.5, cex = cex.axis.labels, font = 2)
  
  par(op)
  layout(1)
}

## ----fig.width=8, fig.height=6, message=TRUE------------------------------------------------------
x <- c(0, 1, 2, 2.2, 2.8, 3.5, 6, 8, 10, 10.1, 12.8, 13, 14.8, 15, 15.5)

# fix overlap, return debugging information
set.seed(10101)
z <- fixOverlap(x, thresh = 0.73, method = 'S', trace = TRUE)

# check convergence
z$converged

# inspect algorithm trace
tracePlot(x, z)

# trace log
# B: boundary condition violation
# O: rank (order) violation
# +: accepted perturbation
# -: rejected perturbation
table(z$log)

## ----fig.width=8, fig.height=6, message=TRUE------------------------------------------------------
# fix overlap, return debugging information
set.seed(101010)
x <- sort(runif(10, min = 2.5, max = 3.5))

# widen boundary conditions
z <- fixOverlap(x, thresh = 0.2, trace = TRUE, min.x = 0, max.x = 10, maxIter = 2000, adj = 0.05)

# check convergence
z$converged

# inspect algorithm trace
tracePlot(x, z)

## -------------------------------------------------------------------------------------------------
# reset plotSPC() options
options(.aqp.plotSPC.args = NULL)

