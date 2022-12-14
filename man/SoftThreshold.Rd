% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasicMath.R
\name{SoftThreshold}
\alias{SoftThreshold}
\title{Soft Thresholding Function}
\usage{
SoftThreshold(x, c = 10)
}
\arguments{
\item{x}{A value which may be a large outlier.}

\item{c}{Maximum Value which can be achieved.}
}
\value{
Soft threshold of input.  Used to trim outliers.
}
\description{
Soft Thresholding Function
}
\examples{
SoftThreshold(100)
}
