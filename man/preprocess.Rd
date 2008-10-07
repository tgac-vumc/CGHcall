\name{preprocess}
\alias{preprocess}
\title{ Preprocess arrayCGH data }
\description{
	This function preprocesses your aCGH data so it can be processed by other functions without errors.
}
\usage{
preprocess(input, maxmiss = 30, nchrom = 22, ...)
}
\arguments{
  \item{input}{ Object of class \code{\link[CGHbase:class.cghRaw]{cghRaw}}. }
  \item{maxmiss}{ Maximum percentage of missing values per row. }
  \item{nchrom}{ Number of chromosomes. }
  \item{\dots}{ Arguments for \code{\link[impute]{impute.knn}} from the impute package. }
}
\details{
  This function performs the following actions on arrayCGH data:
  \itemize{
    \item Filter out data with missing position information.
    \item Remove data on chromosomes larger than nchrom.
    \item Remove rows with more than maxmiss percentage missing values.
    \item Imputes missing values using the \code{\link[impute]{impute.knn}} function from the impute package.
  }
}
\value{
  This function returns a dataframe in the same format as the input with missing values imputed.
}
\references{ 
  Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor Hastie, Robert Tibshirani, David Botstein, and Russ B. Altman (2001). Missing value estimation methods for DNA microarrays. \emph{Bioinformatics, 17}, 520-525.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
  data(WiltingRaw)
  preprocessed <- preprocess(WiltingRaw)
}
\keyword{ misc }
