\name{msaMuscle}
\alias{msaMuscle}
\title{Multiple Sequence Alignment with MUSCLE}
\description{
  This function calls the multiple sequence alignment
  algorithm MUSCLE.
}
\usage{
    msaMuscle(inputSeqs, cluster="default", gapOpening="default",
              gapExtension="default", maxiters="default",
              substitutionMatrix="default",
              type="default", order=c("aligned", "input"),
              verbose=FALSE, help=FALSE, ...)
}
\arguments{
  \item{inputSeqs}{input sequences; see \code{\link{msa}}.
    In the original MUSCLE implementation, this
    parameter is called \code{-in}.}
  \item{cluster}{The clustering method which should be used.
    Possible values are \code{"upgma"}, \code{"upgmamax"},
    \code{"upgmamin"}, \code{"upgmb"}, and \code{"neighborjoining"}.}
  \item{gapOpening}{gap opening penalty; the default 
    is 400 for DNA sequences and 420 for RNA sequences. The default
    for amino acid sequences depends on the profile score settings:
    for the setting \code{le=TRUE}, the default is 2.9, for
    \code{sp=TRUE}, the default is 1,439, and for \code{sv=TRUE},
    the default is 300. Note that these defaults may not be suitable
    if custom substitution matrices are being used. In such a case,
    a sensible choice of gap penalties that fits well to the
    substitution matrix must be made.}
  \item{gapExtension}{gap extension penalty; the default is 0.}
  \item{maxiters}{maximum number of iterations; the default is 16.
    In the original MUSCLE implementation, it is also possible
    to set \code{maxiters} to 0 which leads to an (out of memory) error. 
    Therefore, \code{maxiters=0} is not allowed in \code{msaMuscle}.}
  \item{substitutionMatrix}{substitution matrix for scoring matches and
    mismatches; can be a real matrix or a file name
    If the file interface is used, matrices have to be in NCBI-format.  
    The original MUSCLE implementation also accepts matrices in WU_BLAST
    (AB_BLAST) format, but, due to copyright restrictions,
    this format is not supported by \code{msaMuscle}.}
  \item{type}{type of the input sequences \code{inputSeqs};
    see \code{\link{msa}}.}
  \item{order}{how the sequences should be ordered in the output object
    (see \code{\link{msa}} for more details); the original MUSCLE
    implementation does not allow for preserving the order of input
    sequences. The \code{msaMuscle} function realizes this functionality
    by reverse matching of sequence names. Therefore, the sequences
    need to have unique names. If the sequences do not have
    names or if the names are not unique, the \code{\link{msaMuscle}}
    function assignes generic unique names \code{"Seq1"}-\code{Seqn}
    to the sequences and issues a warning. The choice \code{"input"}
    is not available at all for sequence data that is
    read directly from a FASTA file.}
  \item{verbose}{if \code{TRUE}, the algorithm displays detailed
    information and progress messages.}
  \item{help}{if \code{TRUE}, information about algorithm-specific
    parameters is displayed. In this case, no multiple sequence
    alignment is performed and the function quits after displaying
    the additional help information.}
  \item{...}{further parameters specific to MUSCLE;
    An overview of parameters that are available in this interface
    is shown when calling \code{msaMuscle} with \code{help=TRUE}.
    For more details, see also the documentation of MUSCLE.}
}
\details{This is a function providing the MUSCLE multiple alignment
  algorithm as an R function. It can be used for various types of
  sequence data (see \code{inputSeqs} argument above). Parameters that
  are common to all multiple sequences alignments provided by the
  \pkg{msa} package are explicitly provided by the function and named
  in the same for all algorithms. Most other parameters that are
  specific to MUSCLE can be passed to MUSCLE via additional
  arguments (see argument \code{help} above).

  For a note on the order of output sequences and direct reading from
  FASTA files, see \code{\link{msa}}.
}
\value{
   Depending on the type of sequences for which it was called,
   \code{msaMuscle} returns a \code{\linkS4class{MsaAAMultipleAlignment}}, 
   \code{\linkS4class{MsaDNAMultipleAlignment}}, or
   \code{\linkS4class{MsaRNAMultipleAlignment}} object. 
   If called with \code{help=TRUE}, \code{msaMuscle} returns
   an invisible \code{NULL}.
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

  \url{http://www.drive5.com/muscle/muscle.html}
  
  Edgar, R. C. (2004) MUSCLE: multiple sequence alignment with high 
  accuracy and high throughput.
  \emph{Nucleic Acids Res.} \bold{32}(5):1792-1797. DOI:
  \doi{10.1093/nar/gkh340}.

  Edgar, R. C. (2004) MUSCLE: a multiple sequence alignment method 
  with reduced time and space complexity.
  \emph{BMC Bioinformatics} \bold{5}:113. DOI:
  \doi{10.1186/1471-2105-5-113}.
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

## call msaMuscle with default values
msaMuscle(mySeqs)

## call msaMuscle with custom parameters
msaMuscle(mySeqs, gapOpening=12, gapExtension=3, maxiters=16,
          cluster="upgmamax", SUEFF=0.4, brenner=FALSE,
          order="input", verbose=FALSE)

## call msaMuscle with a custom substitution matrix
data(PAM120)
msaMuscle(mySeqs, substitutionMatrix=PAM120)
}
\keyword{manip}
