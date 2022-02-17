\name{msaConsensusSequence}
\alias{msaConsensusSequence}
\alias{msaConsensusSequence,MultipleAlignment-method}
\alias{msaConsensusSequence,matrix-method}
\title{Computation of Consensus Sequence from Multiple Alignment}
\description{
  This method computes a consensus sequence from a multiple alignment
  or a previously computed consensus matrix. Currently, two different
  ways of these computations are available.
}
\usage{
\S4method{msaConsensusSequence}{matrix}(x, type=c("Biostrings", "upperlower"),
    thresh=c(80, 20), ignoreGaps=FALSE, ...)
\S4method{msaConsensusSequence}{MultipleAlignment}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{MultipleAlignment}}
    (which includes objects of classes
    \code{\linkS4class{MsaAAMultipleAlignment}},
    \code{\linkS4class{MsaDNAMultipleAlignment}}, and
    \code{\linkS4class{MsaRNAMultipleAlignment}}) or a previously
    computed consensus matrix (see details below).}
  \item{type}{a character string specifying how to compute the consensus
    sequence. Currently, types \code{"Biostrings"} and
    \code{"upperlower"} are available (see details below).}
  \item{thresh}{a decreasing two-element numeric vector of numbers
    between 0 and 100 corresponding to the two conservation thresholds.
    Only relevant for \code{type="upperlower"} (see details below),
    otherwise ignored.}
  \item{ignoreGaps}{a logical (default: \code{FALSE})
    indicating whether gaps should be
    considered when computing the consensus sequence. 
    Only relevant for \code{type="upperlower"} (see details below),
    otherwise ignored.}
  \item{...}{when the method is called for a
    \code{\linkS4class{MultipleAlignment}} object, the consensus matrix is
    computed and, including all further arguments, passed
    passed on to the \code{msaConsensusSequence} method with signature
    \code{matrix}. The method with signature \code{matrix} forwards
    additional arguments to the \code{\link{consensusString}} method
    from the \pkg{Biostrings} package if \code{type="Biostrings"}.}
  }
\details{
  The method takes a \code{\linkS4class{MultipleAlignment}} object or a
  previously computed consensus matrix and computes a consensus
  sequence. For \code{type="Biostrings"}, the method
  \code{\link{consensusString}} from the \pkg{Biostrings} package is
  used to compute the consensus sequence. For \code{type="upperlower"},
  two thresholds (argument \code{thresh}, see above) are used to
  compute the consensus sequence:
  \itemize{
    \item{If the relative frequency of the most frequent letter at a
      given position is at least as large as the first threshold (default:
      80\%), then this most frequent letter is used for the consensus
      sequence at this position as it is.}
    \item{If the relative frequency of the most frequent letter at a
      given position is smaller than the first threshold, but at least
      as large as the second threshold (default: 20\%), then this
      most frequent letter is used for the consensus sequence at this
      position, but converted to lower case.}
    \item{If the relative frequency of the most frequent letter in a
      column is even smaller than the second threshold, then a dot is
      used for the consensus sequence at this position.}
    \item{If \code{ignoreGaps=FALSE} (which is the default),
      gaps are treated like all other
      letters except for the fact that obviously no lowercase conversion
      takes place in case that the relative frequency is between the
      two thresholds. If \code{ignoreGaps=TRUE}, gaps are ignored
      completely, and the consensus sequence is computed from the
      non-gap letters only.}
  }

  If the consensus matrix of a multiple alignment of nucleotide
  sequences contains rows labeled \sQuote{+} and/or \sQuote{.}, these
  rows are ignored.
}
\value{
  The function returns a character string with the consensus sequence.
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
  \code{\linkS4class{MsaMetaData}},
  \code{\linkS4class{MultipleAlignment}},
  \code{\link{consensusString}}
}
\examples{
## read sequences
filepath <- system.file("examples", "HemoglobinAA.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)

## perform multiple alignment
myAlignment <- msa(mySeqs)

## regular consensus sequence using consensusString() method from the
## 'Biostrings' package
msaConsensusSequence(myAlignment)

## use the other method
msaConsensusSequence(myAlignment, type="upperlower")

## use the other method with custom parameters
msaConsensusSequence(myAlignment, type="upperlower", thresh=c(50, 20),
                     ignoreGaps=TRUE)

## compute a consensus matrix first
conMat <- consensusMatrix(myAlignment)
msaConsensusSequence(conMat)
}
\keyword{manip}
