\name{segmentData}
\alias{segmentData}
\title{ Breakpoint detection for arrayCGH data. }
\description{
  A wrapper function to run existing breakpoint detection algorithms on arrayCGH data. Currently only DNAcopy is implemented.
}
\usage{
segmentData(input, method = "DNAcopy", ...)
}
\arguments{
  \item{input}{ Object of class \code{\link[CGHbase:class.cghRaw]{cghRaw}}. }
  \item{method}{ The method to be used for breakpoint detection. Currently only 'DNAcopy' is supported, which will run the \code{\link[DNAcopy]{segment}} function.}
  \item{\dots}{ Arguments for \code{\link[DNAcopy]{segment}}. }
}
\details{
	See \code{\link[DNAcopy]{segment}} for details on the algorithm.
}
\value{
  This function returns a dataframe in the same format as the input with segmented arrayCGH data.
}
\references{ 
	Venkatraman, A.S., Olshen, A.B. (2007). A faster circulary binary segmentation algorithm for the analysis of array CGH data. \emph{Bioinformatics, 23}, 657-663.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
  data(WiltingNorm)
  \dontrun{segmented.data <- segmentData(WiltingNorm, alpha=0.02)}
}
\keyword{ misc }
