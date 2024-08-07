% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profile_compare.R
\docType{methods}
\name{pc}
\alias{pc}
\alias{pc.SPC}
\alias{profile_compare}
\alias{profile_compare,SoilProfileCollection-method}
\alias{profile_compare,data.frame-method}
\title{Deprecated: Numerical Soil Profile Comparison}
\usage{
pc(
  s,
  vars,
  max_d,
  k,
  filter = NULL,
  sample_interval = NA,
  replace_na = TRUE,
  add_soil_flag = TRUE,
  return_depth_distances = FALSE,
  strict_hz_eval = FALSE,
  progress = "none",
  plot.depth.matrix = FALSE,
  rescale.result = FALSE,
  verbose = FALSE
)

pc.SPC(s, vars, rescale.result = FALSE, ...)

\S4method{profile_compare}{SoilProfileCollection}(s, vars, rescale.result = FALSE, ...)

\S4method{profile_compare}{data.frame}(
  s,
  vars,
  max_d,
  k,
  filter = NULL,
  sample_interval = NA,
  replace_na = TRUE,
  add_soil_flag = TRUE,
  return_depth_distances = FALSE,
  strict_hz_eval = FALSE,
  progress = "none",
  plot.depth.matrix = FALSE,
  rescale.result = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{s}{A SoilProfileCollection}

\item{vars}{Variable names}

\item{max_d}{depth-slices up to this depth are considered in the comparison}

\item{k}{a depth-weighting coeficient, use '0' for no depth-weighting (see
examples below)}

\item{filter}{an index used to determine which horizons (rows) are included
in the analysis}

\item{sample_interval}{use every n-th depth slice instead of every depth
slice, useful for working with > 1000 profiles at a time}

\item{replace_na}{if TRUE, missing data are replaced by maximum
dissimilarity (TRUE)}

\item{add_soil_flag}{The algorithm will generate a 'soil'/'non-soil' matrix
for use when comparing soil profiles with large differences in depth (TRUE).
See details section below.}

\item{return_depth_distances}{return intermediate, depth-wise dissimilarity
results (FALSE)}

\item{strict_hz_eval}{should horizons be strictly checked for internal
self-consistency? (FALSE)}

\item{progress}{'none' (default)}

\item{plot.depth.matrix}{should a plot of the 'soil'/'non-soil' matrix be
returned (FALSE)}

\item{rescale.result}{Rescale result? Default: FALSE}

\item{verbose}{extra debug output (FALSE)}

\item{...}{Additional arguments passed to \code{profile_compare()}}
}
\value{
A dissimilarity matrix object of class 'dissimilarity, dist',
optionally scaled by max(D).
}
\description{
See \code{NCSP()} for the modern interface to numerical soil profile comparison/classification.
}
\details{
Performs a numerical comparison of soil profiles using named properties,
based on a weighted, summed, depth-segment-aligned dissimilarity
calculation. If \code{s} is a \code{\link{SoilProfileCollection}},
site-level variables (2 or more) can also be used. The site-level and
horizon-level dissimilarity matrices are then re-scaled and averaged.

Variability in soil depth can interfere significantly with the calculation
of between-profile dissimilarity-- what is the numerical ``distance'' (or
dissimilarity) between a slice of soil from profile A and the corresponding,
but missing, slice from a shallower profile B? Gower's distance metric would
yield a NULL distance, despite the fact that intuition suggests otherwise:
shallower soils should be more dissimilar from deeper soils. For example,
when a 25 cm deep profile is compared with a 50 cm deep profile, numerical
distances are only accumulated for the first 25 cm of soil (distances from
26 - 50 cm are NULL). When summed, the total distance between these profiles
will generally be less than the distance between two profiles of equal
depth. Our algorithm has an option (setting replace_na=TRUE) to replace NULL
distances with the maximum distance between any pair of profiles for the
current depth slice. In this way, the numerical distance between a slice of
soil and a corresponding slice of non-soil reflects the fact that these two
materials should be treated very differently (i.e. maximum dissimilarity).

This alternative calculation of dissimilarities between soil and non-soil
slices solves the problem of comparing shallow profiles with deeper
profiles. However, it can result in a new problem: distances calculated
between two shallow profiles will be erroneously inflated beyond the extent
of either profile's depth. Our algorithm has an additional option (setting
add_soil_flag=TRUE) that will preserve NULL distances between slices when
both slices represent non-soil material. With this option enabled, shallow
profiles will only accumulate mutual dissimilarity to the depth of the
deeper profile.

Note that when the \code{add_soil_flag} option is enabled (default), slices
are classified as 'soil' down to the maximum depth to which at least one of
variables used in the dissimilarity calculation is not NA. This will cause
problems when profiles within a collection contain all NAs within the
columns used to determine dissimilarity. An approach for identifying and
removing these kind of profiles is presented in the examples section below.

A notice is issued if there are any NA values within the matrix used for
distance calculations, as these values are optionally replaced by the max
dissimilarity.

Our approach builds on the work of (Moore, 1972) and the previously
mentioned depth-slicing algorithm.
}
\examples{
# See ?NCSP for examples
}
\references{
\itemize{
\item D.E. Beaudette, P. Roudier, A.T. O'Geen, Algorithms for quantitative pedology: A toolkit for soil scientists, Computers & Geosciences, Volume 52, 2013, Pages 258-268, ISSN 0098-3004, \doi{10.1016/j.cageo.2012.10.020}.
\item Moore, A.; Russell, J. & Ward, W. Numerical analysis of soils: A comparison of three soil profile models with field classification. Journal of Soil Science, 1972, 23, 194-209.
}
}
\seealso{
\code{\link[=dice]{dice()}}, \code{\link[cluster:daisy]{cluster::daisy()}}
}
\author{
Dylan E. Beaudette
}
\keyword{manip}
\keyword{methods}
