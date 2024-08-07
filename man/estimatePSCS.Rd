% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimatePSCS.R
\name{estimatePSCS}
\alias{estimatePSCS}
\title{Estimate boundaries of the particle size control section (U.S Soil Taxonomy;
12th edition)}
\usage{
estimatePSCS(
  p,
  hzdesgn = hzdesgnname(p, required = TRUE),
  clay.attr = hzmetaname(p, "clay", required = TRUE),
  texcl.attr = hztexclname(p, required = TRUE),
  tax_order_field = "tax_order",
  bottom.pattern = "Cr|R|Cd",
  simplify = TRUE,
  ...
)
}
\arguments{
\item{p}{A SoilProfileCollection}

\item{hzdesgn}{Name of the horizon attribute containing the horizon
designation. Default 'hzname'}

\item{clay.attr}{Name of the horizon attribute containing clay contents.
Default 'clay'}

\item{texcl.attr}{Name of the horizon attribute containing textural class
(used for finding sandy textures). Default 'texcl'}

\item{tax_order_field}{Name of the site attribute containing taxonomic
order; for handling PSCS rules for Andisols in lieu of lab data. May be NA
or column missing altogether, in which case Andisol PSC possibility is
ignored.}

\item{bottom.pattern}{Regular expression pattern to match a root-restrictive
contact. Default matches Cr, R or Cd. This argument is passed to both
estimateSoilDepth and getArgillicBounds.}

\item{simplify}{Return a length 2 vector with upper and lower boundary when p has length 1? Default TRUE.}

\item{...}{additional arguments are passed to getArgillicBounds()}
}
\value{
A numeric vector (when \code{simplify=TRUE}) containing the top and bottom depth of the particle
size control section. First value is top, second value is bottom.
If \code{p} contains more than one profile, the result is a data.frame with profile ID plus PSCS top and bottom depths.
}
\description{
Estimates the upper and lower boundary of the particle size control section
by applying a programmatic version of the particle size control section key
from the Keys to Soil Taxonomy (12th edition).
}
\details{
Requires information to identify argillic horizons (clay contents, horizon
designations) with \code{getArgillicBounds()} as well as the presence of
plow layers and surface organic soil material. Any
\code{getArgillicBounds()} arguments may be passed to \code{estimatePSCS}.

Requires information on taxonomic order (to handle andisols).

WARNING: Soils in arenic or grossarenic subgroups, with fragipans, or with
strongly contrasting PSCs may not be classified correctly. The author would
welcome a dataset to develop this functionality for.
}
\examples{

data(sp1, package = 'aqp')
depths(sp1) <- id ~ top + bottom
site(sp1) <- ~ group

# set required metadata
hzdesgnname(sp1) <- 'name'
hztexclname(sp1) <- 'texture'
hzmetaname(sp1, 'clay') <- 'prop'

x <- estimatePSCS(sp1)
x
}
\references{
Soil Survey Staff. 2014. Keys to Soil Taxonomy, 12th ed.
USDA-Natural Resources Conservation Service, Washington, DC.
}
\seealso{
\code{getArgillicBounds}, \code{getSurfaceHorizonDepth}
}
\author{
Andrew Gene Brown
}
\keyword{manip}
