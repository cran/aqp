% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulateColor.R
\name{simulateColor}
\alias{simulateColor}
\title{Simulate Soil Colors}
\usage{
simulateColor(
  method = c("dE00", "proportions", "mvnorm"),
  n,
  parameters,
  SPC = NULL
)
}
\arguments{
\item{method}{simulation method, see details}

\item{n}{number of simulated colors per group}

\item{parameters}{a \code{list}, format depends on \code{method}:
\itemize{
\item \code{proportions}: output from \code{\link[=aggregateColor]{aggregateColor()}}
\item \code{dE00}: formatted as \code{list(m = '7.5YR 3/3', thresh = 5, hues = c('7.5YR'))}
\item \code{mvnorm}: formatted as \code{list(hvc = x)}
}

Where \code{m} is a single representative Munsell chip, \code{thresh} is a threshold specified in CIE2000 color contrast (dE00), \code{hues} is a vector of allowed Munsell hues, and \code{x} is a \code{data.frame} representing columns of Munsell hue, value, and chroma having at least 3 rows.}

\item{SPC}{\code{SoilProfileCollection}, attempt to modify \code{SPC} with simulated colors}
}
\value{
a \code{list}, unless \code{SPC} is specified, then a \code{SoilProfileCollection} object
}
\description{
Simulate plausible soil colors based on several possible parameterization of a "range in characteristics" (RIC). Soil color RIC can be specified by a list of parameters:
\itemize{
\item soil color proportions, as output from \code{\link[=aggregateColor]{aggregateColor()}} -- \code{method = 'proportions'}
\item most likely Munsell color, CIE2000 threshold, and vector of acceptable hues -- \code{method = 'dE00'}
\item \code{data.frame} of Munsell hue, value, and chroma representing observed soil colors -- \code{method = 'mvnorm'}
}
}
\examples{

# restrict examples to 2 cores
data.table::setDTthreads(Sys.getenv("OMP_THREAD_LIMIT", unset = 2))

# m: representative or most likely color
# thresh: dE00 threshold
# hues: allowed Munsell hues
p <- list(
  'A' = list(m = '7.5YR 3/3', thresh = 5, hues = c('7.5YR')),
  'BA' = list(m = '7.5YR 4/4', thresh = 8, hues = c('7.5YR')),
  'Bt1' = list(m = '7.5YR 4/4', thresh = 8, hues = c('5YR', '7.5YR')),
  'Bt2' = list(m = '5YR 4/5', thresh = 8, hues = c('5YR', '7.5YR')),
  'Bt3' = list(m = '10YR 4/6', thresh = 10, hues = c('10YR', '7.5YR')),
  'Cr' = list(m = '2.5G 6/2', thresh = 15, hues = c('2.5G', '2.5GY', '2.5BG'))
)

# simulate
(cols <- simulateColor(method = 'dE00', n = 10, parameters = p))

# preview
previewColors(parseMunsell(unlist(cols)), method = 'MDS')

# another example, this time using a larger dE00 threshold
p <- list(
  'A' = list(m = '7.5YR 3/3', thresh = 20, hues = c('10YR', '7.5YR', '5YR'))
)

# simulate
set.seed(54654)
cols <- simulateColor(method = 'dE00', n = 200, parameters = p)

# flatten
cols <- unlist(cols)

# tabulate, sort: most frequent color should be 7.5YR 3/3
sort(table(cols), decreasing = TRUE)

# review colors
previewColors(parseMunsell(cols))

# what does a dE00 threshold look like on 3 pages of hue?
contrastChart('7.5YR 3/3', hues = c('10YR', '7.5YR', '5YR'), thresh = 20)

}
\author{
D.E. Beaudette
}
