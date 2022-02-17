\name{msaConvert}
\alias{msaConvert}
\title{Convert Multiple Sequence Alignment for Other Packages}
\description{
  This function converts a multiple sequence alignment object
  to formats used in other sequence analysis packages.
}
\usage{
    msaConvert(x,
               type=c("seqinr::alignment", "bios2mds::align",
                      "ape::AAbin", "ape::DNAbin",
                      "phangorn::phyDat", "bio3d::fasta"))
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{MultipleAlignment}}
    (which includes objects of classes
    \code{\linkS4class{MsaAAMultipleAlignment}},
    \code{\linkS4class{MsaDNAMultipleAlignment}}, and
    \code{\linkS4class{MsaRNAMultipleAlignment}})}
  \item{type}{a character string specifying to which type of object
    \code{x} should be converted; currently, the
    values \code{"seqinr::alignment"}, \code{"bios2mds::align"},
    \code{"ape::AAbin"}, \code{"ape::DNAbin"},
    \code{"phangorn::phyDat"}, and \code{"bio3d::fasta"}.}
  }
\details{
  The function converts \code{x} to the class of object as
  specified by the \code{type} argument. The values possible
  for the \code{type} argument follow the same principle
  \emph{pkg}::\emph{cl}, i.e. \code{x} is converted to class
  \emph{cl} as defined in the \emph{pkg} package.

  The conversions for usage by the packages \pkg{seqinr},
  \pkg{bios2mds}, and \pkg{ape} work independently of these
  packages and do not strictly require these packages. They
  need not even be installed. This approach has been chosen
  to avoid abundant dependencies and possible incompatibilities.
  That is also why the standard S3/S4 mechanism of
  \code{as}/\code{as.class} functions is not used.

  The conversion to the \code{phyDat} class can be done
  easily using the \code{as.phyDat} function from the
  \pkg{phangorn} package. The \code{msaConvert} function
  still provides this conversion for the sake of consistency.
  However, this conversion is just a wrapper function around
  the \code{as.phyDat} function from the \pkg{phangorn}
  package. Thus, the \pkg{phangorn} package needs to be installed.

  The conversion \code{"ape::AAbin"}
  only works for multiple alignments of amino acid sequences,
  while the conversions \code{"ape::DNAbin"} and
  \code{"phangorn::phyDat"} only work for multiple alignments
  of DNA sequences. When converting to \code{"ape::AAbin"},
  gaps/dashes are replaced by \sQuote{X}. Moreover, conversions
  to \code{"ape::DNAbin"} also convert all characters to
  lowercase and replace gaps/dashes by \sQuote{n}.
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
  \doi{10.1093/bioinformatics/btv494}.
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
