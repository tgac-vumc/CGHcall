\name{postsegnormalize}
\alias{postsegnormalize}
\title{Post-segmentation normalization}
\description{
  This function normalizes arrayCGH data after segmentation in order to find a better 0-level. 
  }
\usage{
postsegnormalize(input, inter=c(-0.1,0.1))
}
\arguments{
  \item{input}{ Object of class \code{\link{cghSeg}}. }
 \item{inter}{Interval in which the function should search for the normal level.
 }
\details{
 This function recursively searches for the interval containing the most segmented data, decreasing the interval length in each recursion.
 The recursive search makes the post-segmentation normalization robust against local maxima. This function is particularly
 useful for profiles for which, after segmentation, the 0-level does not coincide with many segments. It is more or 
 less harmless to other profiles. We advise to keep the search interval (inter) small, in particular at the positive (gain)
 side to avoid that the 0-level is set to a common gain level. 
 }
\value{
  This function returns a cghSeg object in the same format as the input with post-segmentation-normalized adjusted log2 ratios and segmented values.
}
\author{Mark van de Wiel }
\examples{
  data(Wilting)
  ## Convert to \code{\link{cghRaw}} object
  cgh <- make_cghRaw(Wilting)
  ## First preprocess the data
  raw.data <- preprocess(cgh)
  ## Simple global median normalization for samples with 75\% tumor cells
  perc.tumor <- rep(0.75, 3)
  normalized.data <- normalize(raw.data, cellularity=perc.tumor)  
  ## Segmentation with slightly relaxed significance level to accept change-points.
  ## Note that segmentation can take a long time.
  \dontrun{segmented.data <- segmentData(normalized.data, alpha=0.02)}
  \dontrun{postsegnormalized.data <- postsegnormalize(segmented.data)}
  }
\keyword{ misc }
