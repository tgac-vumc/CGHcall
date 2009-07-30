\name{plot.summary}
\alias{plot.summary}
\title{ Visualization of aCGH profiles. }
\description{
  This function creates a summary plot for aCGH profiles.
}
\usage{
	plot.summary(x, y, ...)
}
\arguments{
  \item{x}{ An object of class \code{\link{cghCall}}.}
  \item{y}{ This argument is not used and should be missing. }
  \item{\dots}{ Arguments \code{plot}. }
}
\details{
  We find plotted on the x-axis the array probes sorted by chromosomal position. The vertical bars represent the average probability that the positions they cover are gained (green bars) or lost (red bars). The green bars represent gains, the red bars represent losses. When 4 levels have been used for calling, amplifications are indicated with a blue tickmark at the top of the plot.  
}
\value{
  This function creates a plot.
}
\references{
  Mark A. van de Wiel, Kyung In Kim, Sjoerd J. Vosse, Wessel N. van Wieringen, Saskia M. Wilting and Bauke Ylstra. CGHcall: calling aberrations for array CGH tumor profiles. \emph{Bioinformatics, 23}, 892-894.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
	 \dontrun{
   data(Wilting)
   rawcgh <- make_cghSeg(Wilting)
   normalized <- normalize(rawcgh)
   segmented <- segmentData(normalized)
   called <- CGHcall(segmented)
   plot.summary(called)
   }
}
\keyword{ misc }
