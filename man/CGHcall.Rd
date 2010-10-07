\name{CGHcall}
\alias{CGHcall}
\title{ Calling aberrations for array CGH tumor profiles. }
\description{
  Calls aberrations for array CGH data using a six state mixture model.
}
\usage{
CGHcall(inputSegmented, prior = "auto", nclass = 4, organism = "human", robustsig="yes", nsegfit=3000, maxnumseg=100, minlsforfit=0.5)
}
\arguments{
  \item{inputSegmented}{ An object of class \code{\link{cghSeg}} }
  \item{prior}{ Options are \code{all}, \code{not all}, or \code{auto}. See details for more information. }
  \item{nclass}{ The number of levels to be used for calling. Either \code{3} (loss, normal, gain) or \code{4} (including amplifications). }
  \item{organism}{ Either \code{human} or \code{other}. This is only used for chromosome arm information when \code{prior} is set to \code{all} or \code{auto} (and samplesize > 20). }
  \item{robustsig}{Options are \code{yes} or \code{no}. \code{yes} enforces a lower bound on the standard deviation of the normal segments}
  \item{nsegfit}{Maximum number of segments used for fitting the mixture model. Posterior probabilities are computed for all segments}
  \item{maxnumseg}{Maximum number of segments per profile used for fitting the model}
  \item{minlsforfit}{Minimum length of the segment (in Mb) to be used for fitting the model}
}
\details{
    Please read the article and the supplementary information for detailed information on the algorithm. 
    The parameter \code{prior} states how the data is used to determine the prior probabilities. When set to \code{all}, the probabilities are determined using the entire genome of each sample. 
    When set to \code{not all} probabilites are determined per chromosome for each sample when \code{organism} is set to \code{other} or per chromosome arm when \code{organism} is \code{human}. 
    The chromosome arm information is taken from the March 2006 version of the UCSC database. When \code{prior} is set to \code{auto}, the way probabilities are determined depends on the sample size. The entire genome is used when the sample size is smaller than 20, otherwise chromosome (arm) information is used.
    Please note that CGHcall uses information from all input data to determine the aberration probabilities. 
    When for example triploid or tetraploid tumors are observed, we advise to run CGHcall separately on 
    those (groups of) samples. Note that robustsig = \code{yes} enforces the sd corresponding to 
    the normal segments to be at least half times the pooled gain/loss sd. Use of "nsegfit" significantly lower computing 
    time with respect to previous CGHcall versions without much accuracy loss. Moreover, "maxnumseg" decreases the 
    impact on the results of profiles with inferior segmentation results. Finally, "minlsforfit" decreases the impact
    of very small aberations (potentially CNVs rather than CNAs) on the fit of the model. Note that always a result for all 
    segments is produced. IN MOST CASES, CGHcall SHOULD BE FOLLOWED BY FUNCTION ExpandCGHcall. 
}
\value{
    This function return a list with five components:
  \item{posteriorfin2}{Matrix containing call probabilities for each segment. First column denotes profile number, followed by k columns with aberration probabilities for each sample, where k is the number of levels used for calling (\code{nclass}).}
  \item{nclone }{Number of clone or probes}
  \item{nc}{Number of samples}
  \item{Number of classes used}
  \item{regionsprof}{Matrix containing information about the segments, 4 colums: profile, start probe, end probe, segmented value}

}

\seealso{\code{\link{ExpandCGHcall}}}

\references{
  Mark A. van de Wiel, Kyung In Kim, Sjoerd J. Vosse, Wessel N. van Wieringen, Saskia M. Wilting and Bauke Ylstra. CGHcall: calling aberrations for array CGH tumor profiles. \emph{Bioinformatics, 23}, 892-894.
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
  normalized.data <- normalize(raw.data, cellularity=perc.tumor)  
  ## Segmentation with slightly relaxed significance level to accept change-points.
  ## Note that segmentation can take a long time.
  \dontrun{segmented.data <- segmentData(normalized.data, alpha=0.02)}
  \dontrun{postsegnormalized.data <- postsegnormalize(segmented.data)}
  ## Call aberrations
  \dontrun{result <- CGHcall(postsegnormalized.data)}
  
  ## Expand to CGHcall object
  \dontrun{result <- ExpandCGHcall(result,postsegnormalized.data)}
}
\keyword{ misc }
