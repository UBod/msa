\name{msaConvert}
\alias{msaConvert}
\title{Convert Multiple Sequence Alignment for Other Packages}
\description{
  This function converts a multiple sequence alignment object
  to formats used in other sequence analysis packages.
}
\usage{
    msaConvert(x,
               type=c("seqinr::alignment", "bios2mds::align"))
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{MultipleAlignment}}
    (which includes objects of classes
    \code{\linkS4class{MsaAAMultipleAlignment}},
    \code{\linkS4class{MsaDNAMultipleAlignment}}, and
    \code{\linkS4class{MsaRNAMultipleAlignment}})}
  \item{type}{a character string specifying to which type of object
    \code{x} should be converted; currently, the two
    values \code{"seqinr::alignment"} and \code{"bios2mds::align"}
    are supported, i.e. an object of class \code{"alignment"} as
    defined by the \pkg{seqinr} package or an object of class
    \code{"align"} as defined in the \pkg{bios2mds} package.}
  }
\details{
  The function converts \code{x} to the class of object
  as specified by the \code{type} argument. Note that this conversion
  happens independently of the packages \pkg{seqinr} and
  \pkg{bios2mds}. More specifically, lists with the respective
  class attributes are returned without actually loading any of
  the two packages. They need not even be installed. This approach
  has been chosen to avoid abundant dependencies and possible
  incompatibilities. That is also why the standard S3/S4
  mechanism of \code{as}/\code{as.class} functions is not
  used.
}
\value{
  The function returns an object of the class as specified by
  the \code{type} argument.
}
\author{Ulrich Bodenhofer <msa@bioinf.jku.at>
}
\references{
  \url{http://www.bioinf.jku.at/software/msa}
  
  U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter
  (2015). msa: an R package for multiple sequence alignment. 
  \emph{Bioinformatics} \bold{31}(24):3997-3999. DOI:
  \href{http://dx.doi.org/10.1093/bioinformatics/btv494}{10.1093/bioinformatics/btv494}.
}
\seealso{\code{\link{msa}}, \code{\linkS4class{MsaAAMultipleAlignment}},
  \code{\linkS4class{MsaDNAMultipleAlignment}},
  \code{\linkS4class{MsaRNAMultipleAlignment}},
  \code{\linkS4class{MsaMetaData}}
}
\examples{
## read sequences
filepath <- system.file("examples", "exampleAA.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)

## perform multiple alignment
myAlignment <- msa(mySeqs)

## convert to an object of class 'alignment' (package 'seqinr')
msaConvert(myAlignment, "seqinr::alignment")

## convert to an object of class 'align' (package 'bios2mds')
msaConvert(myAlignment, "bios2mds::align")
}
\keyword{manip}
