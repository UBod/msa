\name{msaCheckNames}
\alias{msaCheckNames}
\title{Check and fix sequence names}
\description{
  This function checks and fixed sequence names of multiple
  alignment objects if they contain characters
  that might lead to LaTeX problems when using
  \code{\link{msaPrettyPrint}}.
}
\usage{
    msaCheckNames(x, replacement=" ", verbose=TRUE)
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{MultipleAlignment}}
    (which includes objects of classes
    \code{\linkS4class{MsaAAMultipleAlignment}},
    \code{\linkS4class{MsaDNAMultipleAlignment}}, and
    \code{\linkS4class{MsaRNAMultipleAlignment}})}
  \item{replacement}{a character string specifying with which
    character(s) potentially problematic characters should be replaced.}
  \item{verbose}{if \code{TRUE} (default), a warning message is shown
    if potentially problematic characters are found. Otherwise,
    the function silently replaces these characters (see details below).}
}
\details{
  The \pkg{Biostrings} package does not impose any restrictions on
  the names of sequences. Consequently, \pkg{msa} also allows all
  possible ASCII strings as sequence (row) names in multiple alignments.
  As soon as \code{\link{msaPrettyPrint}} is used for pretty-printing
  multiple sequence alignments, however, the sequence names are
  interpreted as plain LaTeX source code. Consequently, LaTeX errors
  may arise because of characters or words in the sequence names that LaTeX
  does not or cannot interpret as plain text correctly. This
  particularly includes appearances of special characters and backslash
  characters in the sequence names.

  The \code{msaCheckNames} function takes a multiple alignment object
  and checks sequence names for possibly problematic characters, which
  are all characters but letters (upper and lower case), digits,
  spaces, commas, colons, semicolons, periods, question and exclamation
  marks, dashes, braces, single quotes, and double quotes.
  All other characters are
  considered problematic. The function allows for both checking and
  fixing the sequence names. If called with \code{verbose=TRUE}
  (default), the function prints a warning if a problematic character is
  found. At the same time, regardless of the \code{verbose} argument,
  the function invisibly returns a copy of \code{x} in whose sequence
  names all problematic characters have been replaced by the string
  that is supplied via the \code{replacement} argument (the default is
  a single space).

  In any case, the best solution is to check sequence names carefully
  and to avoid problematic sequence names from the beginning.
}
\value{
  The function invisibly returns a copy of the argument \code{x}
  (therefore, an object of the same class as \code{x}), but
  with modified sequence/row names (see details above).
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
\seealso{\code{\link{msaPrettyPrint}},
  \code{\linkS4class{MsaAAMultipleAlignment}},
  \code{\linkS4class{MsaDNAMultipleAlignment}},
  \code{\linkS4class{MsaRNAMultipleAlignment}}
}
\examples{
## create toy example
mySeqs <- DNAStringSet(c("ACGATCGATC", "ACGACGATC", "ACGATCCCCC"))
names(mySeqs) <- c("Seq. #1", "Seq. \\2", "Seq. ~3")

## perform multiple alignment
myAlignment <- msa(mySeqs)
myAlignment

## check names
msaCheckNames(myAlignment)

## fix names
myAlignment <- msaCheckNames(myAlignment, replacement="", verbose=FALSE)
myAlignment
}
\keyword{manip}
