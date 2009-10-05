\name{cghRaw}
\docType{class}
\alias{class:cghRaw}
\alias{cghRaw}

\alias{cghRaw-class}

\alias{copynumber,cghRaw-method}
\alias{copynumber<-,cghRaw,matrix-method}

\alias{chromosomes,cghRaw-method}
\alias{bpstart,cghRaw-method}
\alias{bpend,cghRaw-method}

\alias{initialize,cghRaw-method}
\alias{plot,cghRaw,missing-method}

\title{
  Class to contain and describe raw or normalized array comparative genomic hybridization data.
}

\description{
  Container for aCGH data and experimental
  metadata. \code{cghRaw} class is derived from
  \code{\link[Biobase:class.eSet]{eSet}}, and requires a matrix named
  \code{copynumber} as 
  assayData member. Furthermore, columns named \code{Chromosome}, \code{Start}, and \code{End}
  are required as featureData members, containing feature position information.
}
\section{Extends}{
   Directly extends class \code{\link[Biobase:class.eSet]{eSet}}.
}
\section{Creating Objects}{

  \code{new('cghRaw',
    phenoData = [AnnotatedDataFrame],
    experimentData = [MIAME],
    annotation = [character],
    copynumber = [matrix],
    featureData = [AnnotatedDataFrame],
    ...)
  }

  \code{make_cghRaw} is a function to convert a dataframe or textfile to
  an object of class \code{cghRaw}. The input should be either a dataframe 
  or a tabseparated textfile (textfiles must contain a header). The first 
  three columns should contain the name, chromosome and position in bp for 
  each array target respectively. The chromosome and position column must 
  contain numbers only. Following these is a column with log2 ratios for 
  each of your samples. If the input type is a textfile, missing values 
  should be represented as 'NA' or an empty field.

}
\section{Slots}{
  Inherited from \code{eSet}:
   \describe{
      \item{\code{assayData}:}{Contains matrices with equal
	dimensions, and with column number equal to
	\code{nrow(phenoData)}. \code{assayData} must contain a matrix
	\code{copynumber} with rows represening array probes
	and columns representing samples. Additional matrices of
	identical size (e.g., representing measurement errors) may
	also be included in
	\code{assayData}. Class:\code{\link[Biobase:class.AssayData]{AssayData-class}}}  
      \item{\code{phenoData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
      \item{\code{featureData}:}{An
	\code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
	with columns 
	\code{Chromosome}, \code{Start}, and \code{End} containing array element position 
	data.}
      \item{\code{experimentData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
      \item{\code{annotation}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
  }
}
\section{Methods}{

  Class-specific methods.
  \describe{
     \item{\code{copynumber(cghRaw)}, \code{copynumber(cghRaw,matrix)<-}}{Access and
       set elements named \code{copynumber} in the \code{AssayData-class}
       slot.}
     \item{\code{chromosomes}, \code{bpstart}, \code{bpend}}{Access the chromosomal positions stored in \code{featureData}}
     \item{plot}{Create a plot containing log2ratios ordered by chromosomal position}
    }

   See \code{\link[Biobase:class.eSet]{eSet}} for derived
   methods. Annotation functionality is not yet supported.    
}

\author{Sjoerd Vosse}

\seealso{
  \code{\link[Biobase:class.eSet]{eSet-class}},
  \code{\link{cghSeg-class}}, \code{\link{cghCall-class}}
  
}

\examples{
# create an instance of cghRaw
new("cghRaw")

# create an instance of cghRaw from a dataframe
data(Wilting)
rawcgh <- make_cghRaw(Wilting)

# plot the first sample
plot(rawcgh[,1])
# first three chromosomes
plot(rawcgh[chromosomes(rawcgh)==1,1])

# get the copynumber values of the third and fourth sample
log2ratios <- copynumber(rawcgh[,3:4])

# get the names of the samples
sampleNames(rawcgh)

# get the names of the array elements
featureNames(rawcgh)
}
\keyword{classes}
