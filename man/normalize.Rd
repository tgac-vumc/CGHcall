\name{normalize}
\alias{normalize}
\title{ Normalization and cellularity adjustment for arrayCGH data. }
\description{
  This function normalizes arrayCGH data using the global mode or median. It can also adjust for the cellularity of your data.
}
\usage{
normalize(input, method = "median", cellularity = 1, smoothOutliers = TRUE, ...)
}
\arguments{
  \item{input}{ Object of class \code{\link{cghRaw}}. }
  \item{method}{ Normalization method, either 'median', 'mode', or 'none'. }
  \item{cellularity}{ A vector of cellularities ranging from 0 to 1 to define the contamination of your sample with healthy cells (1 = no contamination). See details for more information. }
  \item{smoothOutliers}{ Logical. Indicates whether outliers should be smoothed using the \code{\link[DNAcopy]{smooth.CNA}} function.}
  \item{\dots}{ Arguments for \code{\link[DNAcopy]{smooth.CNA}}. }
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
  ## Convert to \code{\link{cghRaw}} object
  cgh <- make_cghRaw(Wilting)
  ## First preprocess the data
  raw.data <- preprocess(cgh)
  ## Simple global median normalization for samples with 75\% tumor cells
  perc.tumor <- rep(0.75, 3)
  normalized.data <- normalize(raw.data, cellularity=perc.tumor)
}
\keyword{ misc }
