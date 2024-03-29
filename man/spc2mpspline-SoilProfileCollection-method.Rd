% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spc2mpspline.R
\name{spc2mpspline,SoilProfileCollection-method}
\alias{spc2mpspline,SoilProfileCollection-method}
\alias{spc2mpspline}
\title{SoilProfileCollection wrapper for \code{mpspline2::mpspline()}}
\usage{
\S4method{spc2mpspline}{SoilProfileCollection}(
  object,
  var_name = NULL,
  method = c("est_1cm", "est_icm", "est_dcm"),
  pattern = "R|Cr|Cd|qm",
  hzdesgn = NULL,
  ...
)
}
\arguments{
\item{object}{A SoilProfileCollection}

\item{var_name}{Column name in \verb{@horizons} slot of \code{object} containing numeric values to spline}

\item{method}{Options include "est_1cm" (default; 1cm estimates), "est_icm" (estimates over original layer boundaries), "est_dcm" (estimates over constant interval, specified with \code{d} argument to \code{mpspline3::mpspline()}). Default value for \code{d} is \code{c(0, 5, 15, 30, 60, 100, 200)}.}

\item{pattern}{Regex pattern to match for bottom of profile (passed to \code{minDepthOf()}) default: "R|Cr|Cd|qm"; only used if \code{hzdesgn} is specified}

\item{hzdesgn}{Column name in \verb{@horizons} slot of \code{object} containing horizon designations default: \code{NULL}}

\item{...}{Additional arguments to \code{mpspline2::mpspline()}}
}
\value{
A SoilProfileCollection with 1cm slices. Spline variables are in columns prefixed with "spline_" and RMSE/RMSE_IQR are in columns prefixed with "rmse_". If any profiles were removed from the collection, their profile IDs are stored in \code{attr(result, 'removed')}.
}
\description{
Generate mass-preserving splines for any numeric attribute in a SoilProfileCollection using \code{mpspline2::mpspline()}. mpspline2 implements the method described in Bishop et al. (1999). Currently this function only works with a single \code{var_name} at a time.
}
\details{
This function now relies on the missing data checks provided by the mpspline2 package. See \code{attr(..., 'removed')} to see whole profiles that were removed from the set. Horizons containing \code{NA} in the property of interest are dropped with a message.

Data completeness is assessed and the input SoilProfileCollection is filtered and truncated to create a container for the results from \code{mpspline2::mpspline()}.
}
\examples{
if (requireNamespace("mpspline2")) {
  data(sp1)
  depths(sp1) <- id ~ top + bottom

  res <- spc2mpspline(sp1, "prop")

  plotSPC(res[1:5,], color = "prop_spline", divide.hz = FALSE)
}
}
\references{
Bishop, T.F.A., McBratney, A.B., Laslett, G.M. (1999) Modelling soil attribute depth functions with equal-area quadratic smoothing splines. Geoderma 91(1–2), pp. 27-45. \doi{https://doi.org/10.1016/S0016-7061(99)00003-8}

O'Brien, Lauren (2022). mpspline2: Mass-Preserving Spline Functions for Soil Data. R package version 0.1.6. \url{https://cran.r-project.org/package=mpspline2}
}
\author{
Andrew G. Brown
}
