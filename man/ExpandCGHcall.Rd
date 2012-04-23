\name{ExpandCGHcall}
\alias{ExpandCGHcall}
\title{Expands result fron CGHcall to CGHcall object. }
\description{
  Expands result from \code{\link{CGHcall}} function to CGHcall object.
}
\usage{
ExpandCGHcall(listcall,inputSegmented, digits=3, divide=4, memeff = FALSE, fileoutpre="Callobj_",CellularityCorrectSeg=TRUE)
}
\arguments{
  \item{listcall}{List object; output of function \code{\link{CGHcall}}}
  \item{inputSegmented}{ An object of class \code{\link{cghSeg}} }
  \item{digits}{Number of decimal digits to be saved in the resulting call object. Allows for saving storage space}
  \item{divide}{Number of batches to divide the work load in. Larger values saves memory, but requires more computing time}
  \item{memeff}{When set to TRUE, memory efficient mode is used: results are written in batches to multiple external files. If FALSE, one output object is provided.}
  \item{fileoutpre}{Only relevant when memeff=TRUE. Define prefix for output file names}
  \item{CellularityCorrectSeg}{If TRUE, corrects segmented and normalized values for cellularity as well}
}
\details{This function is new in version 2.7.0. It allows more memory efficient handling of large data objects. If R crashes because of memory problem, we advise 
to set memeff = TRUE and increase the value of divide.  When multiple files are output (in case of memeff=TRUE) the function combine may be used to combine CGHcall objects.
}

\seealso{\code{\link{CGHcall}}, \code{\link{cghCall-class}}}
\value{
    An object of class \code{\link{cghCall-class}} either as one object (when memeff = FALSE) or as multiple objects stored in .Rdata files in the working directory (when memeff = FALSE)}
\references{
  Mark A. van de Wiel, Kyung In Kim, Sjoerd J. Vosse, Wessel N. van Wieringen, Saskia M. Wilting and Bauke Ylstra. CGHcall: calling aberrations for array CGH tumor 
  profiles. \emph{Bioinformatics, 23}, 892-894.
}
\author{ Sjoerd Vosse & Mark van de Wiel }
\examples{
  data(Wilting)
  ## Convert to \code{\link{cghRaw}} object
  cgh <- make_cghRaw(Wilting)
  print(cgh)
  ## First preprocess the data
  raw.data <- preprocess(cgh)
  ## Simple global median normalization for samples with 75\% tumor cells
  perc.tumor <- rep(0.75, 3)
  normalized.data <- normalize(raw.data)
  ## Segmentation with slightly relaxed significance level to accept change-points.
  ## Note that segmentation can take a long time.
  \dontrun{segmented.data <- segmentData(normalized.data, alpha=0.02)}
  \dontrun{postsegnormalized.data <- postsegnormalize(segmented.data)}
  ## Call aberrations
  \dontrun{result <- CGHcall(postsegnormalized.data, cellularity=perc.tumor)}
  \dontrun{result <- ExpandCGHcall(result,postsegnormalized.data)}
}
\keyword{ misc }
