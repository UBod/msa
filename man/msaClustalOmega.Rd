\name{msaClustalOmega}
\alias{msaClustalOmega}
\title{Multiple Sequence Alignment with ClustalOmega}
\description{
  This function calls the multiple sequence alignment
  algorithm ClustalOmega.
}
\usage{
    msaClustalOmega(inputSeqs, cluster="default",
                    gapOpening="default", gapExtension="default",
                    maxiters="default",  substitutionMatrix="default",
                    type="default", order=c("aligned", "input"),
                    verbose=FALSE, help=FALSE, ...) 
}
\arguments{
  \item{inputSeqs}{input sequences; see \code{\link{msa}}.
    In the original ClustalOmega implementation, this
    parameter is called \code{infile}.}
  \item{cluster}{The cluster size which should be used. The default is 100.  
    In the original ClustalOmega implementation, this parameter is called 
    \code{cluster-size}.}
  \item{gapOpening,gapExtension}{ClustalOmega currently does
    not allow to adjust gap penalties; these arguments are only for
    future extensions and consistency with the other algorithms
    and \code{\link{msa}}. However, setting these parameters to values
    other than \code{"default"} will result in a warning.}
  \item{maxiters}{maximum number of iterations; the default value is 0
    (no limitation). In the original ClustalOmega implementation, this
    parameter is called \code{iterations}.}
  \item{substitutionMatrix}{name of substitution matrix for scoring matches and
    mismatches; can be one of the choices
    \code{"BLOSUM30"}, \code{"BLOSUM40"}, \code{"BLOSUM50"},
    \code{"BLOSUM65"}, \code{"BLOSUM80"}, and \code{"Gonnet"}.
    This parameter is a new feature - the original ClustalOmega
    implementation does not allow for using a custom substitution matrix.}
  \item{type}{type of the input sequences \code{inputSeqs};
    see \code{\link{msa}}.}
  \item{order}{how the sequences should be ordered in the output object
    (see \code{\link{msa}}); in the original ClustalW implementation, this
    parameter is called \code{output-order}.}
  \item{verbose}{if \code{TRUE}, the algorithm displays detailed
    information and progress messages.}
  \item{help}{if \code{TRUE}, information about algorithm-specific
    parameters is displayed. In this case, no multiple sequence
    alignment is performed and the function quits after displaying
    the additional help information.}
  \item{...}{further parameters specific to ClustalOmega;
    An overview of parameters that are available in this interface
    is shown when calling \code{msaClustalOmega} with \code{help=TRUE}.
    For more details, see also the documentation of ClustalOmega.}
}
\details{This is a function providing the ClustalOmega multiple alignment
  algorithm as an R function. It can be used for various types of
  sequence data (see \code{inputSeqs} argument above). Parameters that
  are common to all multiple sequences alignments provided by the
  \pkg{msa} package are explicitly provided by the function and named
  in the same for all algorithms. Most other parameters that are
  specific to ClustalOmega can be passed to ClustalOmega via additional
  arguments (see argument \code{help} above).

  Since ClustalOmega only allows for using built-in amino acid
  substitution matrices, it is hardly useful for multiple alignments
  of nucleotide sequences.

  For a note on the order of output sequences and direct reading from
  FASTA files, see \code{\link{msa}}.
}
\value{
   Depending on the type of sequences for which it was called,
   \code{msaClustalOmega} returns a
   \code{\linkS4class{MsaAAMultipleAlignment}}, 
   \code{\linkS4class{MsaDNAMultipleAlignment}}, or
   \code{\linkS4class{MsaRNAMultipleAlignment}} object.
   If called with \code{help=TRUE}, \code{msaClustalOmega} returns
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

  \url{http://www.clustal.org/omega/README}
  
  Sievers, F., Wilm, A., Dineen, D., Gibson, T. J., Karplus, K., Li, W.,
  Lopez, R., McWilliam, H., Remmert, M., Soeding, J., Thompson, J. D.,
  and Higgins, D. G. (2011) Fast, scalable generation of high-quality
  protein multiple sequence alignments using Clustal Omega.
  \emph{Mol. Syst. Biol.} \bold{7}:539. DOI:
  \doi{10.1038/msb.2011.75}.
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

## call msaClustalOmega with default values
msaClustalOmega(mySeqs)

## call msaClustalOmega with custom parameters
msaClustalOmega(mySeqs, auto=FALSE, cluster=120, dealign=FALSE,
                useKimura=FALSE, order="input", verbose=FALSE)
}
\keyword{manip}
