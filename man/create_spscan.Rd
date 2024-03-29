% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spscan.test.R
\name{create_spscan}
\alias{create_spscan}
\title{Create spscan object}
\usage{
create_spscan(
  tobs,
  pvalue,
  alpha,
  noc_info,
  coords,
  z,
  N,
  N1,
  x,
  d,
  maxd,
  longlat,
  nsim
)
}
\arguments{
\item{tobs}{The vector of observed test statistics in
descending order for non-overlapping windows.}

\item{pvalue}{The pvalues associated with \code{tobs}.}

\item{alpha}{The significance level to determine whether a cluster is
signficant.  Default is 0.1.}

\item{noc_info}{An object from the
\code{\link[smerc]{noc_nn}} function.}

\item{z}{An indicator vector indicating which
observations are cases.}

\item{N}{Number of event locations}

\item{N1}{Number of case events}

\item{x}{A \code{\link[spatstat.geom]{ppp}} object with marks for the case
and control groups.}

\item{d}{A matrix of intercentroid distances for the
event locations.}

\item{maxd}{The radius of the largest possible cluster to consider.  Default
is \code{NULL}, i.e., half the maximum intercentroid distance.}

\item{longlat}{A logical value indicating whether
Euclidean distance (\code{FALSE}) or Great Circle
(WGS84 ellipsoid, \code{FALSE}) should be used. Default
is \code{FALSE}, i.e., Euclidean distance.}

\item{nsim}{The number of simulations from which to compute the p-value.  A
non-negative integer.  Default is 499.}
}
\value{
An \code{spscan} object.
}
\description{
Create spscan object
}
\keyword{internal}
