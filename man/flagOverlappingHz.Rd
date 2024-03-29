% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flagOverlappingHz.R
\name{flagOverlappingHz}
\alias{flagOverlappingHz}
\title{Flag perfectly overlapping horizons within a SoilProfileCollection}
\usage{
flagOverlappingHz(x)
}
\arguments{
\item{x}{a \code{SoilProfileCollection} object}
}
\value{
logical vector with length (and order) matching the horizons of \code{x}
}
\description{
Flag perfectly overlapping horizons within a SoilProfileCollection
}
\examples{

# two overlapping horizons
z <- data.frame(
  id = 'SPC',
  top = c(0, 25, 25, 50, 75, 100, 100),
  bottom = c(25, 50, 50, 75, 100, 125, 125)
)

# init SPC
depths(z) <- id ~ top + bottom

# flag perfectly overlapping horizons
z$.overlapFlag <- flagOverlappingHz(z)

# thematic sketches
plotSPC(z, color = '.overlapFlag', hz.depths = TRUE, 
depth.axis = FALSE, cex.names = 0.85)

}
\author{
D.E. Beaudette
}
