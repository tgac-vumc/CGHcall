\name{chromosomes}
\alias{chromosomes}
\alias{bpstart}
\alias{bpend}

\title{Retrieve feature position data from cgh objects.}
\description{
  These generic functions access the position data 
  stored in the featureData of an object derived from the \code{\link{cghRaw-class}}, 
  \code{\link{cghSeg-class}} or \code{\link{cghCall-class}}.
}
\usage{
chromosomes(object)
bpstart(object)
bpend(object)
}
\arguments{
  \item{object}{Object derived from class \code{cghRaw}, \code{cghSeg}, or \code{cghCall}}
}
\value{
  \code{chromosomes} returns a vector of chromosome numbers;
  \code{bpstart} returns a vector of basepair start positions;
  \code{bpend} returns a vector of basepair end positions;
}
\author{Sjoerd Vosse}

\seealso{\code{\link{cghRaw-class}}, \code{\link{cghSeg-class}}, \code{\link{cghCall-class}}}

\keyword{manip}
