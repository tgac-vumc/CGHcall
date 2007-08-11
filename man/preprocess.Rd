\name{preprocess}
\alias{preprocess}
\title{ Preprocess arrayCGH data }
\description{
  This function performs the following actions on arrayCGH data:
  \itemize{
    \item Filter out data with missing position information.
    \item Remove data on chromosomes larger than nchrom.
    \item Remove rows with more than maxmiss percentage missing values.
    \item Imputes missing values using the \code{\link[impute]{impute.knn}} function from the impute package.
  }
}
\usage{
preprocess(input, type = "file", maxmiss = 30, nchrom = 22, ...)
}
\arguments{
  \item{input}{ Either the name of a file or a dataframe. See details for the format. }
  \item{type}{ What kind of data format is used as input? Either 'dataframe' or 'file'. }
  \item{maxmiss}{ Maximum percentage of missing values per row. }
  \item{nchrom}{ Number of chromosomes. }
  \item{\dots}{ Arguments for \code{\link[impute]{impute.knn}} from the impute package. }
}
\details{
  The input should be either a dataframe or a tabseparated textfile (textfiles must contain a header). The first three columns should contain the name, chromosome and position in bp for each array target respectively. The chromosome and position column must contain numbers only. Following these is a column with raw log2 ratios for each of your samples. If the input type is a textfile, missing values should be represented as 'NA' or an empty field.

}
\value{
  This function returns a dataframe in the same format as the input with missing values imputed.
}
\references{ 
  Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor Hastie, Robert Tibshirani, David Botstein, and Russ B. Altman (2001). Missing value estimation methods for DNA microarrays. \emph{Bioinformatics, 17}, 520-525.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
  data(Wilting)
  raw.data <- preprocess(Wilting, type="dataframe")
}
\keyword{ misc }
