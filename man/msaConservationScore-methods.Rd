\name{msaConservationScore}
\alias{msaConservationScore}
\alias{msaConservationScore,MultipleAlignment-method}
\alias{msaConservationScore,matrix-method}
\title{Computation of Conservation Scores from Multiple Alignment}
\description{
  This method computes a vector of conservation scores from a
  multiple alignment or a previously computed consensus matrix.
}
\usage{
\S4method{msaConservationScore}{matrix}(x, substitutionMatrix, gapVsGap=NULL,
    ...)
\S4method{msaConservationScore}{MultipleAlignment}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{MultipleAlignment}}
    (which includes objects of classes
    \code{\linkS4class{MsaAAMultipleAlignment}},
    \code{\linkS4class{MsaDNAMultipleAlignment}}, and
    \code{\linkS4class{MsaRNAMultipleAlignment}}) or a previously
    computed consensus matrix (see details below).}
  \item{substitutionMatrix}{substitution matrix (see
    details below).}
  \item{gapVsGap}{score to use for aligning gaps versus gaps (see
    details below).}
  \item{...}{when the method is called for a
    \code{\linkS4class{MultipleAlignment}} object, the consensus matrix is
    computed and, including all further arguments,
    passed on to the \code{msaConservationScore} method with signature
    \code{matrix}. This method passes all further arguments on to the
    \code{\link{msaConsensusSequence}} method to customize the way
    the consensus sequence is computed.}
  }
\details{
  The method takes a \code{\linkS4class{MultipleAlignment}} object or a
  previously computed consensus matrix and computes the sum of pairwise
  scores for all positions of the alignment. For computing these scores,
  it is compulsory to specify a substitution/scoring matrix. This matrix
  must be provided as a \code{\link{matrix}} object. This can either be
  one of the ready-made matrices provided by the \pkg{Biostrings}
  package (e.g. \code{\link{BLOSUM62}}) or any other hand-crafted
  matrix. In the latter case, the following restrictions apply:
  \itemize{
    \item{The matrix must be quadratic.}
    \item{For reasonable results, the matrix should be symmetric (note
      that this is not checked).}
    \item{Rows and columns must be named and the order of letters/symbols
      in row names and column names must be identical.}
    \item{All letters/symbols occurring in the multiple alignment,
      including gaps \sQuote{-}, must also be found in the row/column
      names of the substitution matrix. For consistency with the
      matrices from the \pkg{Biostrings} package, the
      row/column corresponding to gap penalties
      may also be labeled \sQuote{*} instead of \sQuote{-}.}
  }
  So, nucleotide substitution matrices created by
  \code{\link{nucleotideSubstitutionMatrix}} can be used for multiple
  alignments of nucleotide sequences, but must be
  completed with gap penalty rows and columns (see example below).

  If the consensus matrix of a multiple alignment of nucleotide
  sequences contains rows labeled \sQuote{+} and/or \sQuote{.}, these
  rows are ignored.
  
  The parameter \code{gapVsGap} can be used to control how
  pairs of gaps are scored. If \code{gapVsGap=NULL} (default), the
  corresponding diagonal entry of the substitution matrix is used as is.
  In the BLOSUM matrices, this is usually a positive value, which may
  not make sense under all circumstances. Therefore, the parameter
  \code{gapVsGap} can be set to an alternative value, e.g. 0 for
  ignoring gap-gap pairs.

  The method, in any case, returns a vector of scores that is as long
  as the alignment/consensus matrix has columns. The names of the vector
  entries are the corresponding positions of the consensus sequence of
  the alignment. How this consensus sequence is computed, can be
  controlled with additional arguments that are passed on to the
  \code{\link{msaConsensusSequence}} method.
}
\value{
  The function returns a vector as long as the alignment/consensus
  matrix has columns. The vector is named with the consensus sequence
  (see details above).
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
  \code{\link{msaConsensusSequence}}
}
\examples{
## read sequences
filepath <- system.file("examples", "HemoglobinAA.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)

## perform multiple alignment
myAlignment <- msa(mySeqs)

## compute consensus scores using the BLOSUM62 matrix
data(BLOSUM62)
msaConservationScore(myAlignment, BLOSUM62)

## compute consensus scores using the BLOSUM62 matrix
## without scoring gap-gap pairs and using a different consensus sequence
msaConservationScore(myAlignment, BLOSUM62, gapVsGap=0,
                     type="upperlower")

## compute a consensus matrix first
conMat <- consensusMatrix(myAlignment)
data(PAM250)
msaConservationScore(conMat, PAM250, gapVsGap=0)


## DNA example
filepath <- system.file("examples", "exampleDNA.fasta", package="msa")
mySeqs <- readDNAStringSet(filepath)

## perform multiple alignment
myAlignment <- msa(mySeqs)

## create substitution matrix with gap penalty -8
mat <- nucleotideSubstitutionMatrix(4, -1)
mat <- cbind(rbind(mat, "-"=-8), "-"=-8)

## compute consensus scores using this matrix
msaConservationScore(myAlignment, mat, gapVsGap=0)
}
\keyword{manip}
