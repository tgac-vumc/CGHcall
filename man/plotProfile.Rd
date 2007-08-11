\name{plotProfile}
\alias{plotProfile}
\alias{plotSummary}
\title{ Visualization of aCGH profiles. }
\description{
  This function creates a plot for aCGH profiles, including raw log2ratios, segments and calls.
}
\usage{
plotProfile(CGHcall.result, samples="all", export="no")
plotSummary(CGHcall.result, samples="all", export="no")
}
\arguments{
  \item{CGHcall.result}{ The list of results that is returned by \code{\link[CGHcall]{CGHcall}}.}
  \item{samples}{ Which samples to plot. Either 'all' or a numerical vector. }
  \item{export}{ Whether a postscript file should be exported. Either 'no' or a filename.}
}
\details{
  For \code{plotProfile} we find plotted on the x-axis the array probes sorted by chromosomal position. The black dots represent the normalized log2ratios (values are on the y-axis on the right side). The horizontal blue lines are the segments. The vertical bars represent the probability that the probes they cover are aberrations. The green bars represent gains, the red bars represent losses. If a bar crosses the middle horizontal line (probability of 0.5), the probes in this region are said to be aberrated. When 4 levels have been used for calling, amplifications are indicated with a blue tickmark at the top of the plot.
  
  For \code{plotSummary} we find plotted on the x-axis the array probes sorted by chromosomal position. The vertical bars represent the average probability that the positions they cover are gained (green bars) or lost (red bars). The green bars represent gains, the red bars represent losses. When 4 levels have been used for calling, amplifications are indicated with a blue tickmark at the top of the plot.  
}
\value{
  This function creates a plot.
}
\references{
  Mark A. van de Wiel, Kyung In Kim, Sjoerd J. Vosse, Wessel N. van Wieringen, Saskia M. Wilting and Bauke Ylstra. CGHcall: calling aberrations for array CGH tumor profiles. \emph{Bioinformatics, 23}, 892-894.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
  data(Wilting)
  ## First preprocess the data
  raw.data <- preprocess(Wilting, type="dataframe")
  ## Simple global median normalization for samples with 75\% tumor cells
  perc.tumor <- rep(0.75, 3)
  normalized.data <- normalize(raw.data, cellularity=perc.tumor)
  ## Segmentation with slightly relaxed significance level to accept change-points.
  ## Note that segmentation can take a long time.
  \dontrun{segmented.data <- segmentData(normalized.data, alpha=0.02)}
  ## Call aberrations
  \dontrun{result <- CGHcall(normalized.data, segmented.data)}
  ## Create a plot for the first sample and export it to a file
  \dontrun{plotProfile(result, samples=1, export='plot_1.ps')}
  ## Create a summary plot including the first and third sample
  \dontrun{plotSummary(result, samples=c(1,3), export='summary_1_3.ps')}  
}
\keyword{ misc }
