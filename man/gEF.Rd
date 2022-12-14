% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EstimateCurvature.R
\name{gEF}
\alias{gEF}
\title{Estimating Function Value}
\usage{
gEF(kappa, dxy, dxz, dyz, dxm)
}
\arguments{
\item{kappa}{Curvature of Embedding Space}

\item{dxy}{Triangle Side Length}

\item{dxz}{Triangle Side Length}

\item{dyz}{Triangle Side Length}

\item{dxm}{Length from point x to the midpoint of y and z}
}
\value{
Value of estimating function. This will equal 0 when triangle can be embedded in the space of curvature kappa
}
\description{
Estimating Function Value
}
\examples{
gEF(0,1,1,1,sqrt(3/4))
}
