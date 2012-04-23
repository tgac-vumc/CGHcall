\name{segmentData}
\alias{segmentData}
\title{ Breakpoint detection for arrayCGH data. }
\description{
  A wrapper function to run existing breakpoint detection algorithms on arrayCGH data. Currently only DNAcopy is implemented.
}
\usage{
segmentData(input, clen=10, relSDlong=3, method = "DNAcopy", ...)
}
\arguments{
  \item{input}{Object of class \code{\link[CGHbase:class.cghRaw]{cghRaw}}. }
  \item{clen}{Boundary for short vs long segments, in number of features}
  \item{relSDlong}{Relative undo sd for long segments. See details. }
  \item{method}{The method to be used for breakpoint detection. Currently only DNAcopy is supported, which will run the \code{\link[DNAcopy]{segment}} function.}
  \item{\dots}{Arguments for \code{\link[DNAcopy]{segment}}. }
}
\details{
	See \code{\link[DNAcopy]{segment}} for details on the algorithm. About \code{clen} and \code{relSDlong}: 
these are only relevant when \link[DNAcopy] option \code{undo.splits}=\code{sdundo} is set, in combination with \link[DNAcopy] option \code{undo.SD}. 
\code{relSDlong} provides the undo sd for long segments, which equals \code{undo.SD/relSDlong}. \code{undo.SD} is then used for short segments. 
In the example below, short segments are considered to contain 
less or equal to \code{clen}=10 features. The example below undoes splits for two consecutive short segments if these are less than \code{undo.SD}=3 sd apart, 
while it undoes splits for two long segments if these are less than \code{undo.SD/relSDlong}=3/3=1 sd apart. If, for two consecutive segements, one is short and one is long, 
splits are undone in the same way as for two short segments.
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
  \dontrun{segmented.data <- segmentData(WiltingNorm, alpha=0.02,clen=10,relSDlong=3,undo.SD=3,undo.splits="sdundo")}
}
\keyword{ misc }
