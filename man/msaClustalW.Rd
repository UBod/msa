\name{msaClustalW}
\alias{msaClustalW}
\title{Multiple Sequence Alignment with ClustalW}
\description{
  This function calls the multiple sequence alignment
  algorithm ClustalW.
}
\usage{
    msaClustalW(inputSeqs, cluster="default", gapOpening="default", 
                gapExtension="default", maxiters="default", 
                substitutionMatrix="default", type="default", 
                order=c("aligned", "input"), verbose=FALSE,
                help=FALSE, ...)
}
\arguments{
  \item{inputSeqs}{input sequences; see \code{\link{msa}}.
    In the original ClustalW implementation, this
    parameter is called \code{infile}.}
  \item{cluster}{The clustering method which should be used.
    Possible values are \code{"nj"} (default) and \code{"upgma"}.
    In the original ClustalW implementation, this parameter is called 
    \code{clustering}. Please note that \code{cluster="upgma"} leads
    to an unidentified error on Windows with R 4.0.x that even crashes
    the entire R session.}
  \item{gapOpening}{gap opening penalty; the default value for
    nucleotide sequences is 15.0, the default value for 
    amino acid sequences is 10.0.}
  \item{gapExtension}{gap extension penalty; the default value for
    nucleotide sequences is 6.66, the default value for 
    amino acid sequences is 0.2.}
  \item{maxiters}{maximum number of iterations; the default value is 16.
    In the original ClustalW implementation, this parameter is called
    \code{numiters}.}
  \item{substitutionMatrix}{substitution matrix for scoring matches and
    mismatches; can be a real matrix, a file name, or the name of a
    built-in substitution matrix. In the latter case, the choices
    \code{"blosum"}, \code{"pam"}, \code{"gonnet"}, and \code{"id"} are
    supported for amino acid sequences. For aligning nucleotide
    sequences, the choices \code{"iub"} and \code{"clustalw"} are
    possible. The parameter \code{dnamatrix} can also be used instead
    for the sake of backwards compatibility.
    The valid choices for this parameter are \code{"iub"} and
    \code{"clustalw"}. In the original ClustalW implementation, this
    parameter is called \code{matrix}.}
  \item{type}{type of the input sequences \code{inputSeqs};
    see \code{\link{msa}}.}
  \item{order}{how the sequences should be ordered in the output object
    (see \code{\link{msa}}); in the original ClustalW implementation, this
    parameter is called \code{outorder}.}
  \item{verbose}{if \code{TRUE}, the algorithm displays detailed
    information and progress messages.}
  \item{help}{if \code{TRUE}, information about algorithm-specific
    parameters is displayed. In this case, no multiple sequence
    alignment is performed and the function quits after displaying
    the additional help information.}
  \item{...}{further parameters specific to ClustalW;
    An overview of parameters that are available in this interface
    is shown when calling \code{msaClustalW} with \code{help=TRUE}.
    For more details, see also the documentation of ClustalW.}
}
\details{This is a function providing the ClustalW multiple alignment
  algorithm as an R function. It can be used for various types of
  sequence data (see \code{inputSeqs} argument above). Parameters that
  are common to all multiple sequences alignments provided by the
  \pkg{msa} package are explicitly provided by the function and named
  in the same for all algorithms. Most other parameters that are
  specific to ClustalW can be passed to ClustalW via additional
  arguments (see argument \code{help} above).

  For a note on the order of output sequences and direct reading from
  FASTA files, see \code{\link{msa}}.
}
\value{
   Depending on the type of sequences for which it was called,
   \code{msaClustalW} returns a \code{\linkS4class{MsaAAMultipleAlignment}}, 
   \code{\linkS4class{MsaDNAMultipleAlignment}}, or
   \code{\linkS4class{MsaRNAMultipleAlignment}} object. 
   If called with \code{help=TRUE}, \code{msaClustalW} returns
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
 
  \url{http://www.clustal.org/download/clustalw_help.txt}

  Thompson, J. D., Higgins, D. G., and Gibson, T. J. (1994)
  CLUSTAL W: improving the sensitivity of progressive multiple sequence
  alignment through sequence weighting, position-specific gap penalties
  and weight matrix choice.
  \emph{Nucleic Acids Res.} \bold{22}(22):4673-4680. DOI:
  \doi{10.1093/nar/22.22.4673}.
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

## call msaClustalW with default values
msaClustalW(mySeqs)

## call msaClustalW with custom parameters
msaClustalW(mySeqs, gapOpening=1, gapExtension=1, maxiters=16,
            kimura=FALSE, order="input", maxdiv=23)
}
\keyword{manip}
