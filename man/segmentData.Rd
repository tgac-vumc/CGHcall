\name{segmentData}
\alias{segmentData}
\title{ Breakpoint detection for arrayCGH data. }
\description{
  A wrapper function to run existing breakpoint detection algorithms on arrayCGH data. Currently only DNAcopy is implemented.
}
\usage{
segmentData(input, type = "dataframe", method = "DNAcopy", ...)
}
\arguments{
  \item{input}{ Either the name of a file or a dataframe. See details for the format. }
  \item{type}{ What kind of data format is used as input? Either 'dataframe' or 'file'. }
  \item{method}{ The method to be used for breakpoint detection. Currently only 'DNAcopy' is supported, which will run the \code{\link[DNAcopy]{segment}} function.}
  \item{\dots}{ Arguments for \code{\link[DNAcopy]{segment}}. }
}
\details{
  The input should be either a dataframe or a tabseparated textfile (textfiles must contain a header). The first three columns should contain the name, chromosome and position in bp for each array target respectively. The chromosome and position column must contain numbers only. Following these is a column with normalized log2 ratios for each of your samples. If the input type is a textfile, missing values should be represented as 'NA' or an empty field.
}
\value{
  This function returns a dataframe in the same format as the input with segmented arrayCGH data.
}
\references{ 
	Venkatraman, A.S., Olshen, A.B. (2007). A faster circulary binary segmentation algorithm for the analysis of array CGH data. \emph{Bioinformatics, 23}, 657-663.
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
}
\keyword{ misc }
