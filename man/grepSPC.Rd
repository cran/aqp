% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SoilProfileCollection-methods.R
\name{grepSPC}
\alias{grepSPC}
\alias{grepSPC,SoilProfileCollection-method}
\title{Subset SPC with pattern-matching for text-based attributes}
\usage{
grepSPC(object, attr, pattern, ...)
}
\arguments{
\item{object}{A SoilProfileCollection}

\item{attr}{A character vector (column in object) for matching patterns against.}

\item{pattern}{REGEX pattern to match in \code{attr}}

\item{...}{Additional arguments are passed to \code{grep()}}
}
\value{
A SoilProfileCollection.
}
\description{
\code{grepSPC()} is a shorthand function for subsetting \code{SoilProfileCollection} objects. For example, by \code{filter(grepl(spc, ...))} or \code{filter(stringr::str_detect(spc, ...))}. It provides pattern matching for a single text-based site or horizon level attribute.
}
\author{
Andrew G. Brown.
}
