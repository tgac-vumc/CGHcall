\name{copynumber}
\alias{copynumber}
\alias{copynumber<-}
\alias{segmented}
\alias{segmented<-}
\alias{calls}
\alias{calls<-}

\title{Retrieve copynumber data from cgh objects.}
\description{
  These generic functions access the copynumber values
  of assay data stored in an object derived from the \code{\link{cghRaw-class}}, 
  \code{\link{cghSeg-class}} or \code{\link{cghCall-class}}.
}
\usage{
copynumber(object)
copynumber(object) <- value
segmented(object)
segmented(object) <- value
calls(object)
calls(object) <- value
}
\arguments{
  \item{object}{Object derived from class \code{cghRaw}, \code{cghSeg}, or \code{cghCall}}
  \item{value}{Matrix with rows representing features and columns samples.}
}
\value{
  \code{copynumber} returns a matrix of copynumber values;
}
\author{Sjoerd Vosse}

\seealso{\code{\link{cghRaw-class}}, \code{\link{cghSeg-class}}, \code{\link{cghCall-class}}}

\keyword{manip}
