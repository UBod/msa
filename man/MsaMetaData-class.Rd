\name{MsaMetaData-class}
\docType{class}
\alias{class:MsaMetaData}
\alias{MsaMetaData-class}
\alias{MsaMetaData}
\alias{params,MsaMetaData-method}
\alias{version}
\alias{version,MsaMetaData-method}
\title{Class \code{MsaMetaData}}
\description{S4 class for storing metadata about multiple sequence
  alignment results}
\section{Objects}{
  Objects of this virtual class are not be created and used directly. This
  is an auxiliary class used by the classes
  \code{\linkS4class{MsaAAMultipleAlignment}},
  \code{\linkS4class{MsaDNAMultipleAlignment}}, and
  \code{\linkS4class{MsaRNAMultipleAlignment}}
}
\section{Slots}{
The following slots are defined for \code{MsaMetaData} objects:
  \describe{
    \item{\code{version}:}{slot in which information is stored with which
      algorithm the multiple alignment has been computed along with its
      version number.}
    \item{\code{params}:}{list in which the parameters are stored with
      which the multiple alignment algorithm has been executed.}
    \item{\code{call}:}{the matched call with which the object was
      created}
  }
}
\section{Methods}{
  \describe{
    \item{\code{version(object)}:}{accessor to the \code{version} slot}
    \item{\code{params(x)}:}{accessor to the \code{params} slot}
  }
}
\author{Enrico Bonatesta and Christoph Horejs-Kainrath
  <msa@bioinf.jku.at>
}
\references{
  \url{http://www.bioinf.jku.at/software/msa}
  
  U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter
  (2015). msa: an R package for multiple sequence alignment. 
  \emph{Bioinformatics} \bold{31}(24):3997-3999. DOI:
  \doi{10.1093/bioinformatics/btv494}.
}
\seealso{\code{\link{msa}}, \code{\link{msaClustalW}},
  \code{\link{msaClustalOmega}}, \code{\link{msaMuscle}},
  \code{\linkS4class{MsaAAMultipleAlignment}},
  \code{\linkS4class{MsaDNAMultipleAlignment}},
  \code{\linkS4class{MsaRNAMultipleAlignment}}
}
\examples{
## read sequences
filepath <- system.file("examples", "exampleAA.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)

## simple call with default values
myAlignment <- msaClustalOmega(mySeqs)

## show the algorithm version with which the results were created
version(myAlignment)

## show the results
show(myAlignment)

## print the results
print(myAlignment, show="alignment")
print(myAlignment, show=c("alignment", "version"))
print(myAlignment, show="standardParams")
print(myAlignment, show="algParams")
print(myAlignment, show=c("call", "version"))

## show the params
params(myAlignment)
}
\keyword{class}