% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explainPlotSPC.R
\name{explainPlotSPC}
\alias{explainPlotSPC}
\title{Visual Explanation for \code{plotSPC}}
\usage{
explainPlotSPC(x, ...)
}
\arguments{
\item{x}{a \code{SoilProfileCollection} object}

\item{\dots}{arguments passed to \code{\link{plotSPC}}}
}
\value{
a list of internally-used ordering vectors and graphical offsets / scaling factors

Attempt to fix overlapping sketches when using relative horizontal spacing.

This is a very simple optimization algorithm for adjusting horizontal sketch
positions until affected profiles are farther apart than a given threshold.
Rank-ordering and boundary conditions are enforced on the adjustments.
Failure to converge within \code{maxIter} results in an integer sequence.
}
\description{
Create a visual explanation for the many arguments to \code{plotSPC}. Call this function instead of \code{plotSPC}, all objects after \code{x} are passed on to \code{plotSPC}. Nearly all of the figures in the \href{https://ncss-tech.github.io/AQP/aqp/aqp-intro.html}{Introduction to SoilProfileCollection Objects tutorial} are created with this function.
}
\examples{

# sample data
data(sp4)
depths(sp4) <- id ~ top + bottom

# proposed vector of relative positions, overlap likely
pos <- c(1, 1.1, 3, 4, 5, 5.2, 7, 8, 9, 10)

# try it
explainPlotSPC(sp4, name='name', relative.pos=pos)

# attempt to fix using an integer sequence, short-circut will prevent adjustments
explainPlotSPC(sp4, name='name', relative.pos=fixOverlap(1:10, trace=TRUE))

# attempt to adjust using defaults
explainPlotSPC(sp4, name='name', relative.pos=fixOverlap(pos, trace=TRUE))

# attempt to adjust and tinker with defaults
explainPlotSPC(sp4, name='name', relative.pos=fixOverlap(pos, adj = 0.2, trace=TRUE))

# repeatable adjustments
set.seed(10101)
explainPlotSPC(sp4, name='name', relative.pos=fixOverlap(pos, thresh = 0.7, trace=TRUE))

# more complex adjustments required
pos <- c(1, 2, 3, 3.3, 5, 5.1, 5.5, 8, 9, 9.1)

# tinker
explainPlotSPC(sp4, name='name', relative.pos=pos)
explainPlotSPC(sp4, name='name', relative.pos=fixOverlap(pos, trace=TRUE))
explainPlotSPC(sp4, name='name', relative.pos=fixOverlap(pos, thresh = 0.7, trace=TRUE))
explainPlotSPC(sp4, name='name', relative.pos=fixOverlap(pos, thres=0.7, adj = 0.2, trace=TRUE))

# no solution possible given these constraints
explainPlotSPC(sp4, name='name', relative.pos=fixOverlap(pos, thres=1, adj = 0.2, trace=TRUE))


}
\seealso{
\code{\link{plotSPC}}
}
\author{
D.E. Beaudette
}
\keyword{manip}
