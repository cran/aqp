% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getArgillicBounds.R
\name{argillic.clay.increase.depth}
\alias{argillic.clay.increase.depth}
\title{Return upper boundary of argillic horizon}
\usage{
argillic.clay.increase.depth(p, clay.attr = "clay")
}
\arguments{
\item{p}{A single-profile \code{SoilProfileCollection} object.}

\item{clay.attr}{OPTIONAL: horizon attribute name referring to clay content.
default: \code{clay}}
}
\value{
A numeric vector containing top depth of argillic horizon, if
present, or NA.
}
\description{
Returns the top depth of the argillic horizon as a numeric vector.
}
\details{
Uses \code{crit.clay.argillic} to determine threshold clay increase, and
\code{get.increase.matrix} to determine where increase is met within a
vertical distance of 30 cm.
}
\examples{

data(sp1, package = 'aqp')
depths(sp1) <- id ~ top + bottom
site(sp1) <- ~ group

p <- sp1[1]
attr <- 'prop' # clay contents 
foo <- argillic.clay.increase.depth(p, clay.attr = attr)
foo

}
\seealso{
\code{getArgillicBounds}, \code{get.increase.matrix},
\code{crit.clay.argillic}
}
\author{
Andrew Gene Brown
}
\keyword{manip}
