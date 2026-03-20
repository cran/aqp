## ----setup, echo=FALSE, results='hide', warning=FALSE---------------------------------------------
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  background = '#F7F7F7',
  fig.align = 'center',
  dev = 'png',
  dpi = 96,
  optipng = knitr::hook_optipng,
  comment = "#>"
)

# keep examples from using more than 2 cores
data.table::setDTthreads(Sys.getenv("OMP_THREAD_LIMIT", unset = 2))

options(width = 100, stringsAsFactors = FALSE, timeout = 600)

suppressMessages(library(aqp, quietly = TRUE))

## ----echo = FALSE, fig.width=8, fig.height=4------------------------------------------------------
.m <- '10YR 6/6'

data("munsell.spectra.wide")
w <- munsell.spectra.wide[, 1]
s <- munsell.spectra.wide[, .m]

par(mar = c(4.5, 4.5, 2, 1), cex.axis = 0.75, lend = 2)
plot(w, s, xlab = 'Wavelength (nm)', ylab = 'Reflectance', type = 'b', main = .m, cex = 0.8, pch = 16, axes = FALSE)
rect(xleft = 400, ybottom = 0.35, xright = 450, ytop = 0.4, col = parseMunsell(.m), border = 1, lwd = 1)
axis(side = 1, at = seq(380, 730, by = 20))
axis(side = 2, las = 1)

## -------------------------------------------------------------------------------------------------
# Munsell -> hex color
parseMunsell('5PB 4/6')

# Munsell -> sRGB
parseMunsell('5PB 4/6',  return_triplets = TRUE)

# Munsell -> CIELAB
parseMunsell('5PB 4/6',  returnLAB = TRUE)

# hex color -> Munsell
col2Munsell('#476189FF')

# neutral color
parseMunsell('N 5/')

# non-standard notation
getClosestMunsellChip('3.3YR 4.4/6.1', convertColors = FALSE)

## ----fig.width=10, fig.height=4-------------------------------------------------------------------
# example moist soil colors from the Musick soil series
m <- c("7.5YR 2.5/1", "10YR 3/2", "7.5YR 4/3", "2.5YR 3/6", "2.5YR 3/6", 
       "2.5YR 4/6", "5YR 4/6", "7.5YR 5/4")

mm <- parseMunsell(m, convertColors = FALSE)

d.p <- estimateSoilColor(
  hue = mm$hue, 
  value = mm$value, 
  chroma = mm$chroma, 
  method = 'procrustes', 
  sourceMoistureState = 'moist'
)

d.ols <- estimateSoilColor(
  hue = mm$hue, 
  value = mm$value, 
  chroma = mm$chroma, 
  method = 'ols', 
  sourceMoistureState = 'moist'
)

d.p <- sprintf("%s %s/%s", d.p$hue, d.p$value, d.p$chroma)
d.ols <- sprintf("%s %s/%s", d.ols$hue, d.ols$value, d.ols$chroma)

colorContrastPlot(m, d.p, labels = c('Moist', 'Estimated\nDry'), d.cex = 0.9)

# it is the same
# colorContrastPlot(m, d.ols, labels = c('Moist', 'Estimated\nDry'), d.cex = 0.9)

## ----fig.width=6, fig.height=3--------------------------------------------------------------------
data("Ohz.colors")

Ohz.colors$col <- parseMunsell(Ohz.colors$L1.munsell)

op <- par(mfrow = c(2, 1), mar = c(0.5, 0.5, 1.5, 0))

with(
  Ohz.colors[Ohz.colors$state == 'dry', ],
  soilPalette(colors = col, lab = sprintf("%s\n%s", genhz, L1.munsell), lab.cex = 1)
)
title(main = 'Dry Colors')

with(
  Ohz.colors[Ohz.colors$state == 'moist', ],
  soilPalette(colors = col, lab = sprintf("%s\n%s", genhz, L1.munsell), lab.cex = 1)
)
title(main = 'Moist Colors')

# restore original base graphics state
par(op)

