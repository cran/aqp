% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/previewColors.R
\name{previewColors}
\alias{previewColors}
\title{Preview Colors}
\usage{
previewColors(
  cols,
  method = "grid",
  col.order = NULL,
  nrow = ceiling(sqrt(length(cols))),
  ncol = nrow,
  border.col = "black",
  pt.cex = 2,
  pt.pch = 15
)
}
\arguments{
\item{cols}{vector of R colors}

\item{method}{either "grid", "MDS", or "manual", see details}

\item{col.order}{integer vector used to order colors}

\item{nrow}{number of rows used by "grid" method}

\item{ncol}{number of columns used by "grid" method}

\item{border.col}{border color used by "grid" method}

\item{pt.cex}{point scaling factor used by "MDS" method}

\item{pt.pch}{point symbol used by "MDS" method}
}
\value{
When \code{method = "grid" or "manual"} a vector of color order is returned. When \code{method = "MDS"}, the output from \code{MASS::cmdscale}.
}
\description{
Preview colors as a grid arranged according to CIE2000 distances.
}
\details{
Color sorting is based on CIE2000 distances as calculated by \code{farver::compare_colour()}. The "grid" method arranges colors in a rectangular grid with ordering based on divisive hierarchical clustering of the pair-wise distances. Unique colors are used when \code{cols} contains more than 5,000 colors.

The "MDS" method arranges unique colors via classical multidimensional scaling (principal coordinates) via \code{MASS::cmdscale}.

Colors can be manually arranged by supplying a vector of integers to \code{col.order} and setting \code{method='manual'}.
}
\examples{

# example data
data(sp2)

# convert into SoilProfileCollection object
depths(sp2) <- id ~ top + bottom

previewColors(sp2$soil_color)
previewColors(sp2$soil_color, method = 'MDS', pt.cex = 3)

}
\author{
D.E. Beaudette
}
