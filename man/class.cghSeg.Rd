\name{cghSeg}
\docType{class}
\alias{class:cghSeg}
\alias{cghSeg}

\alias{cghSeg-class}

\alias{copynumber,cghSeg-method}
\alias{copynumber<-,cghSeg,matrix-method}
\alias{segmented,cghSeg-method}
\alias{segmented<-,cghSeg,matrix-method}

\alias{chromosomes,cghSeg-method}
\alias{bpstart,cghSeg-method}
\alias{bpend,cghSeg-method}

\alias{initialize,cghSeg-method}
\alias{plot,cghSeg,missing-method}

\title{
  Class to contain and describe segmented array comparative genomic hybridization data.
}

\description{
  Container for aCGH data and experimental
  metadata. \code{cghSeg} class is derived from
  \code{\link[Biobase:class.eSet]{eSet}}, and requires a matrix named \code{copynumber} as well as a 
  matrix named \code{segmented} as assayData members of equal dimensions. 
  Furthermore, columns named \code{Chromosome}, \code{Start}, and \code{End} are 
  required as featureData members, containing feature position information.
}
\section{Extends}{
   Directly extends class \code{\link[Biobase:class.eSet]{eSet}}.
}
\section{Creating Objects}{

  \code{new('cghSeg',
    phenoData = [AnnotatedDataFrame],
    experimentData = [MIAME],
    annotation = [character],
    copynumber = [matrix],
    segmented = [matrix],
    featureData = [AnnotatedDataFrame],
    ...)
  }

  An object of class \code{cghSeg} is generally obtained as output
  from \code{\link{segmentData}}.

}
\section{Slots}{
  Inherited from \code{eSet}:
   \describe{
      \item{\code{assayData}:}{Contains matrices with equal
    dimensions, and with column number equal to
    \code{nrow(phenoData)}. \code{assayData} must contain matrices
    \code{copynumber} and \code{segmented} with rows represening array probes
    and columns representing samples. Additional matrices of
    identical size (e.g., representing measurement errors) may
    also be included in
    \code{assayData}. Class:\code{\link[Biobase:class.AssayData]{AssayData-class}}}  
      \item{\code{phenoData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
      \item{\code{featureData}:}{An
	\code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
	with columns 
    \code{Chromosome}, \code{Start}, and \code{End} containing array
    element position  
    data.}
      \item{\code{experimentData}:}{See
	\code{\link[Biobase:class.eSet]{eSet}}} 
      \item{\code{annotation}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
  }
}
\section{Methods}{

  Class-specific methods.
  \describe{
     \item{\code{copynumber(cghSeg)},
       \code{copynumber(cghSeg,matrix)<-}}{Access and 
       set elements named \code{copynumber} in the \code{AssayData-class}
       slot.}
     \item{\code{segmented(cghSeg)},
       \code{segmented(cghSeg,matrix)<-}}{Access and 
       set elements named \code{segmented} in the \code{AssayData-class}
       slot.} 
     \item{\code{chromosomes}, \code{bpstart}, \code{bpend}}{Access the
       chromosomal positions stored in \code{featureData}} 
     \item{plot}{Create a plot containing log2ratios and segments
       ordered by chromosomal position. TWO EXTRA OPTIONS PLUS DEFAULTS:
     
     dotres=10. Every dotres-th log2-ratio is plotted. dotres=1 plots
     all data. However, higher values save a lot of space and  
     allow quicker browsing of the plots. ylimit=c(-2,5): limits of the
     y-axis} 
    }

   See \code{\link[Biobase:class.eSet]{eSet}} for derived methods.   
}

\author{Sjoerd Vosse}

\seealso{
  \code{\link[Biobase:class.eSet]{eSet-class}},
  \code{\link{cghRaw-class}}, \code{\link{cghCall-class}} 
}

\examples{
# create an instance of cghSeg
new("cghSeg")

# create an instance of cghSeg through 'segmentData'
\dontrun{
   data(Wilting)
   rawcgh <- make_cghSeg(Wilting)
   normalized <- normalize(rawcgh)
   segmented <- segmentData(normalized)

   # plot the first sample. Default only every 10th log2-ratio is plotted (dotres=10). Adjust using dotres= option below. 
   plot(segmented[,1])
   # first three chromosomes
   plot(segmented[chromosomes(segmented)<=3,1])

   # get the copynumber values of the third and fourth sample
   log2ratios <- copynumber(segmented[,3:4])

   # get the names of the samples
   sampleNames(segmented)

   # get the names of the array elements
   featureNames(segmented)
}
}

\keyword{classes}
