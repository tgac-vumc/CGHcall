\name{probloss}
\alias{probloss}
\alias{probnorm}
\alias{probgain}
\alias{probamp}
\alias{probloss<-}
\alias{probnorm<-}
\alias{probgain<-}
\alias{probamp<-}

\title{Retrieve call probabilities from a cghCall object.}
\description{
  These generic functions access the call probabilities from assay data
  stored in a object derived from the \code{\link{cghCall-class}}.
}
\usage{
probloss(object)
probloss(object) <- value
probnorm(object)
probnorm(object) <- value
probgain(object)
probgain(object) <- value
probamp(object)
probamp(object) <- value
}
\arguments{
  \item{object}{Object derived from class \code{cghCall}}
  \item{value}{Matrix with rows representing features and columns samples.}
}
\value{
  \code{probloss} returns matrix of call probabilities;
}
\author{Sjoerd Vosse}

\seealso{\code{\link{cghCall-class}}}

\keyword{manip}
