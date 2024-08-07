% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depthOf.R
\name{depthOf}
\alias{depthOf}
\alias{maxDepthOf}
\alias{minDepthOf}
\title{Get top or bottom depths of horizons matching a regular expression pattern}
\usage{
depthOf(
  p,
  pattern,
  FUN = NULL,
  top = TRUE,
  hzdesgn = hzdesgnname(p, required = TRUE),
  no.contact.depth = NULL,
  no.contact.assigned = NA_real_,
  na.rm = TRUE,
  simplify = TRUE
)

maxDepthOf(
  p,
  pattern,
  top = TRUE,
  hzdesgn = hzdesgnname(p, required = TRUE),
  no.contact.depth = NULL,
  no.contact.assigned = NA,
  na.rm = TRUE,
  simplify = TRUE
)

minDepthOf(
  p,
  pattern,
  top = TRUE,
  hzdesgn = hzdesgnname(p, required = TRUE),
  no.contact.depth = NULL,
  no.contact.assigned = NA,
  na.rm = TRUE,
  simplify = TRUE
)
}
\arguments{
\item{p}{a SoilProfileCollection}

\item{pattern}{a regular expression to match in the horizon designation column. See:\code{hzdesgn}}

\item{FUN}{a function that returns a single value, and takes argument \code{na.rm}}

\item{top}{should the top (TRUE) or bottom (FALSE) depth be returned for matching horizons? Default: \code{TRUE}.}

\item{hzdesgn}{column name containing horizon designations. Default: \code{guessHzDesgnName(p)}}

\item{no.contact.depth}{depth to assume that contact did not occur.}

\item{no.contact.assigned}{depth to assign when a contact did not occur.}

\item{na.rm}{logical. Remove \code{NA}? (default: \code{TRUE})}

\item{simplify}{logical. Return single profile results as vector (default: \code{TRUE}) or \code{data.frame} (\code{FALSE})}
}
\value{
a numeric vector containing specified depth(s) of horizons matching a pattern. If \code{length(p) > 1} then a \emph{data.frame} containing profile ID, horizon ID, top or bottom depths, horizon designation and pattern.
}
\description{
The \code{depthOf} family of functions calculate depth of occurrence of a horizon designation pattern, or any other value that can be coerced to character and matched with a regular expression.

If you need all depths of occurrence for a particular pattern, \code{depthOf} is what you are looking for. \code{minDepthOf} and \code{maxDepthOf} are wrappers around \code{depthOf} that return the minimum and maximum depth. They are all set up to handle missing values and missing "contacts" with the target pattern.
}
\examples{
# construct a fake profile
spc <- data.frame(id=1, taxsubgrp = "Lithic Haploxerepts",
                  hzname   = c("A","AB","Bw","BC","R"),
                  hzdept   = c(0,  20, 32, 42,  49),
                  hzdepb   = c(20, 32, 42, 49, 200),
                  clay     = c(19, 22, 22, 21,  NA),
                  texcl    = c("l","l","l", "l","br"),
                  d_value  = c(5,   5,  5,  6,  NA),
                  m_value  = c(2.5, 3,  3,  4,  NA),
                  m_chroma = c(2,   3,  4,  4,  NA))

# promote to SoilProfileCollection
depths(spc) <- id ~ hzdept + hzdepb
hzdesgnname(spc) <- 'hzname'
hztexclname(spc) <- 'texcl'

# multiple horizons contain B
depthOf(spc, "B")

# deepest top depth of horizon containing B
maxDepthOf(spc, "B")

# shallowest top depth
minDepthOf(spc, "B")

# deepest bottom depth
maxDepthOf(spc, "B", top = FALSE)

# deepest bottom depth above 35cm
maxDepthOf(spc, "B", top = FALSE, no.contact.depth = 35)

# assign infinity (Inf) if B horizon does not start within 10cm
minDepthOf(spc, "B", no.contact.depth = 10, no.contact.assigned = Inf)

}
\author{
Andrew G. Brown
}
