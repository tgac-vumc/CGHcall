\name{normalize}
\alias{normalize}
\title{ Normalization and cellularity adjustment for arrayCGH data. }
\description{
  This function normalizes arrayCGH data using the global mode or median. It can also adjust for the cellularity of your data.
}
\usage{
normalize(input, type = "dataframe", method = "median", cellularity = 1, smoothOutliers = TRUE, ...)
}
\arguments{
  \item{input}{ Either the name of a file or a dataframe. See details for the format. }
  \item{type}{ What kind of data format is used as input? Either 'dataframe' or 'file'. }
  \item{method}{ Normalization method, either 'median', 'mode', or 'none'. }
  \item{cellularity}{ A vector of cellularities ranging from 0 to 1 to define the contamination of your sample with healthy cells (1 = no contamination). See details for more information. }
  \item{smoothOutliers}{ Logical. Indicates whether outliers should be smoothed using the \code{\link[DNAcopy]{smooth.CNA}} function.}
  \item{\dots}{ Arguments for \code{\link[DNAcopy]{smooth.CNA}}. }
}
\details{
  The input should be either a dataframe or a tabseparated textfile (textfiles must contain a header). The first three columns should contain the name, chromosome and position in bp for each array target respectively. The chromosome and position column must contain numbers only. Following these is a column with log2 ratios for each of your samples. If the input type is a textfile, missing values should be represented as 'NA' or an empty field.

  The cellularity parameter should be a vector of length n where n is the number of samples in your dataset. The vector is recycled if there are not enough values in it, or truncated if there are too many. For more information on the correction we refer to section 1.6 of the supplementary information for van de Wiel et al. 2006.
}
\value{
  This function returns a dataframe in the same format as the input with normalized and/or cellularity adjusted log2 ratios.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
  data(Wilting)
  ## First preprocess the data
  raw.data <- preprocess(Wilting, type="dataframe")
  ## Simple global median normalization for samples with 75\% tumor cells
  perc.tumor <- rep(0.75, 3)
  normalized.data <- normalize(raw.data, cellularity=perc.tumor)
}
\keyword{ misc }
