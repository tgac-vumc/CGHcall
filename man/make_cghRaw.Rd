\name{make_cghRaw}
\alias{make_cghRaw}
\title{ Convert a dataframe or textfile to an object of class cghRaw. }
\description{
  This function converts a dataframe of appropriote format to an object of class cghRaw.
}
\usage{
make_cghRaw(input)
}
\arguments{
  \item{input}{ Either a dataframe or character string containing a filename. See details for the format.}
}
\details{
  The input should be either a dataframe or a tabseparated textfile (textfiles must contain a header). The first four columns should contain the name, chromosome and the start and end position in bp for each array target respectively. The chromosome and position column must contain numbers only. Following these is a column with log2 ratios for each of your samples. If the input type is a textfile, missing values should be represented as 'NA' or an empty field.
}
\value{
  This function returns an object of class \code{\link{cghRaw-class}}.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
  data(Wilting)
  ## Convert to \code{\link{cghRaw}} object
  cgh <- make_cghRaw(Wilting)
}
\keyword{ misc }
