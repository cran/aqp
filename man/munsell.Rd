% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-documentation.R
\docType{data}
\name{munsell}
\alias{munsell}
\title{Munsell to sRGB Lookup Table for Common Soil Colors}
\format{
A data.frame with 8825 rows.
\itemize{
\item \code{hue}: Munsell Hue, upper case
\item \code{value}: Munsell Value
\item \code{chroma}: Munsell Chroma
\item \code{r}: sRGB "red" value (0-1)
\item \code{g}: sRGB "green" value (0-1)
\item \code{b}: sRGB "blue" value (0-1)
\item \code{L}: CIELAB "L" coordinate
\item \code{A}: CIELAB "A" coordinate
\item \code{B}: CIELAB "B" coordinate
}
}
\source{
Color chip XYZ values: \url{https://www.rit.edu/science/munsell-color-science-lab-educational-resources#munsell-renotation-data}
}
\usage{
data(munsell)
}
\description{
A lookup table of interpolated Munsell color chips for common soil colors.
}
\details{
See \code{munsell2rgb} for conversion examples. Values
are referenced to the D65 standard illuminant.
}
\examples{

data(munsell)

}
\references{
\itemize{
\item Color conversion equations
\itemize{
\item http://www.brucelindbloom.com/index.html?ColorCalcHelp.html
}
\item Methods used to generate this table
\itemize{
\item http://dx.doi.org/10.1016/j.cageo.2012.10.020
}
}
}
\keyword{datasets}
