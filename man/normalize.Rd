\name{normalize}
\alias{normalize}
\title{Normalization and cellularity adjustment for arrayCGH data. }
\description{
  This function normalizes arrayCGH data using the global mode or median. It can also adjust for the cellularity of your data.
}
\usage{
normalize(input, method = "median", smoothOutliers = TRUE, ...)
}
\arguments{
  \item{input}{Object of class \code{\link[CGHbase:class.cghRaw]{cghRaw}}. }
  \item{method}{Normalization method, either \code{median}, \code{mode}, or \code{none}. }
  \item{smoothOutliers}{Logical. Indicates whether outliers should be smoothed using the \code{\link[DNAcopy]{smooth.CNA}} function.}
  \item{\dots}{Arguments for \code{\link[DNAcopy]{smooth.CNA}}. }
}
\details{
  The cellularity parameter should be a vector of length n where n is the number of samples in your dataset. The vector is recycled if there are not enough values in it, or truncated if there are too many. For more information on the correction we refer to section 1.6 of the supplementary information for van de Wiel et al. 2006.
}
\value{
  This function returns a dataframe in the same format as the input with normalized and/or cellularity adjusted log2 ratios.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
  data(Wilting)
  ## Convert to 'cghRaw' object
  cgh <- make_cghRaw(Wilting)
  ## First preprocess the data
  raw.data <- preprocess(cgh)
  ## Simple global median normalization for samples with 75\% tumor cells
  normalized.data <- normalize(raw.data)
}
\keyword{ misc }
