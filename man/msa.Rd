\name{msa}
\alias{msa}
\title{Unified interface to multiple sequence alignment algorithms}
\description{
  The \code{msa} function provides a unified interface to
  the three multiple sequence alignment algorithms in this package:
  \sQuote{ClustalW}, \sQuote{ClustalOmega}, and \sQuote{MUSCLE}. 
}
\usage{
    msa(inputSeqs, method=c("ClustalW", "ClustalOmega", "Muscle"),
        cluster="default", gapOpening="default",
        gapExtension="default", maxiters="default",
        substitutionMatrix="default", type="default",
        order=c("aligned", "input"), verbose=FALSE, help=FALSE,
        ...)
}
\arguments{
  \item{inputSeqs}{input sequences; this argument can be a character vector,
    an object of class \code{\linkS4class{XStringSet}} (includes the
    classes \code{\linkS4class{AAStringSet}}, \code{\linkS4class{DNAStringSet}},
    and \code{\linkS4class{RNAStringSet}}), or a single character string with a
    file name. In the latter case, the file name is required to have the
    suffix \sQuote{.fa} or \sQuote{.fasta}, and the file must be in
    FASTA format.}
  \item{method}{
    specifies the multiple sequence alignment to be used;
    currently, \code{"ClustalW"}, \code{"ClustalOmega"}, and
    \code{"Muscle"} are supported.}
  \item{cluster}{parameter related to sequence clustering; its
    interpretation and default value depends on the method;
    see \code{\link{msaClustalW}}, \code{\link{msaClustalOmega}}, or
    \code{\link{msaMuscle}} for algorithm-specific information.}
  \item{gapOpening}{gap opening penalty; the defaults are
    specific to the algorithm (see \code{\link{msaClustalW}},
    and \code{\link{msaMuscle}}). Note that the sign of
    this parameter is ignored. The sign is automatically
    adjusted such that the called algorithm penalizes gaps
    instead of rewarding them.}
  \item{gapExtension}{gap extension penalty; the defaults are
    specific to the algorithm (see \code{\link{msaClustalW}},
    and \code{\link{msaMuscle}}). Note that the sign of
    this parameter is ignored. The sign is automatically
    adjusted such that the called algorithm penalizes gaps
    instead of rewarding them.}
  \item{maxiters}{maximum number of iterations; its
    interpretation and default value depends on the method;
    see \code{\link{msaClustalW}}, \code{\link{msaClustalOmega}}, or
    \code{\link{msaMuscle}} for algorithm-specific information.}
  \item{substitutionMatrix}{substitution matrix for scoring matches and
    mismatches; format and defaults depend on the algorithm;
    see \code{\link{msaClustalW}}, \code{\link{msaClustalOmega}}, or
    \code{\link{msaMuscle}} for algorithm-specific information.}
  \item{type}{type of the input sequences \code{inputSeqs}; possible
    values are \code{"dna"}, \code{"rna"}, or \code{"protein"}.  
    In the original ClustalW implementation, this parameter is also called
    \code{-type}; \code{"auto"} is also possible in the original
    ClustalW, but, in this package, \code{"auto"} is deactivated.
    The \code{type} argument is mandatory if \code{inputSeqs} is
    a character vector or the file name of a FASTA file (see above).
    If \code{inputSeqs} is an object of class
    \code{\linkS4class{AAStringSet}}, \code{\linkS4class{DNAStringSet}},
    or \code{\linkS4class{RNAStringSet}}, the type of sequences is
    determined by the class of \code{inputSeqs} and the \code{type}
    parameter is not necessary. If it is nevertheless specified and the
    type does not match the class of \code{inputSeqs}, the function
    stops with an error.}
  \item{order}{how the sequences should be ordered in the output object;
    if \code{"aligned"} is chosen, the sequences are ordered in the way
    the multiple sequence alignment algorithm orders them. If
    \code{"input"} is chosen, the sequences in the output object are
    ordered in the same way as the input sequences. For MUSCLE, the
    choice \code{"input"} is not available for sequence data that is
    read directly from a FASTA file. Even if sequences are supplied
    directly via R, the sequences must have unique names, otherwise
    the input order cannot be recovered. If the sequences do not have
    names or if the names are not unique, the \code{\link{msaMuscle}}
    function assignes generic unique names \code{"Seq1"}-\code{Seqn}
    to the sequences and issues a warning.}
  \item{verbose}{if \code{TRUE}, the algorithm displays detailed
    information and progress messages.}
  \item{help}{if \code{TRUE}, information about algorithm-specific
    parameters is displayed. In this case, no multiple sequence
    alignment is performed and the function quits after displaying
    the additional help information.}
  \item{...}{all other parameters are passed on to the multiple
    sequence algorithm, i.e. to one of the functions
    \code{\link{msaClustalW}}, \code{\link{msaClustalOmega}}, or
    \code{\link{msaMuscle}}. An overview of parameters that are
    available for the chosen method 
    is shown when calling \code{msa} with \code{help=TRUE}.
    For more details, see also the documentation of chosen
    multiple sequence alignment algorithm.}
}
\details{
  \code{msa} is a simple wrapper function that unifies the interfaces of
  the three functions \code{\link{msaClustalW}},
  \code{\link{msaClustalOmega}}, and \code{\link{msaMuscle}}. Which
  function is called, is controlled by the \code{method} argument.

  Note that the input sequences may be reordered by the multiple
  sequence alignment algorithms in order to group together similar
  sequences (see also description of argument \code{order} above).
  So, if the input order should be preserved or if the input order
  should be recovered later, we strongly recommend to always assign
  unique names to the input sequences. As noted in the description
  of the \code{inputSeqs} argument above, all functions, \code{msa()},
  \code{\link{msaClustalW}}, \code{\link{msaClustalOmega}}, and
  \code{\link{msaMuscle}}, also allow
  for direct reading from FASTA files. This is mainly for the reason of
  memory efficiency if the sequence data set is very large. Otherwise,
  we want to encourage users to first read the sequences into the R
  workspace. If sequences are read from a FASTA file
  directly, the order of output sequences is completely under
  the control of the respective
  algorithm and does not allow for checking whether the sequences are
  named uniquely in the FASTA file. The preservation of the input order
  works also for sequence data read from a FASTA file, but only for
  ClustalW and ClustalOmega; MUSCLE does not support this (see also
  argument \code{order} above and \code{\link{msaMuscle}}).
}
\value{
   Depending on the type of sequences for which it was called,
   \code{msa} returns a \code{\linkS4class{MsaAAMultipleAlignment}}, 
   \code{\linkS4class{MsaDNAMultipleAlignment}}, or
   \code{\linkS4class{MsaRNAMultipleAlignment}} object. 
   If called with \code{help=TRUE}, \code{msa} returns
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

  \url{http://www.clustal.org/omega/README}
  
  \url{http://www.drive5.com/muscle/muscle.html}
  
  Thompson, J. D., Higgins, D. G., and Gibson, T. J. (1994)
  CLUSTAL W: improving the sensitivity of progressive multiple sequence
  alignment through sequence weighting, position-specific gap penalties
  and weight matrix choice.
  \emph{Nucleic Acids Res.} \bold{22}(22):4673-4680. DOI:
  \doi{10.1093/nar/22.22.4673}.

  Sievers, F., Wilm, A., Dineen, D., Gibson, T. J., Karplus, K., Li, W.,
  Lopez, R., McWilliam, H., Remmert, M., Soeding, J., Thompson, J. D.,
  and Higgins, D. G. (2011) Fast, scalable generation of high-quality
  protein multiple sequence alignments using Clustal Omega.
  \emph{Mol. Syst. Biol.} \bold{7}:539. DOI:
  \doi{10.1038/msb.2011.75}.

  Edgar, R. C. (2004) MUSCLE: multiple sequence alignment with high 
  accuracy and high throughput.
  \emph{Nucleic Acids Res.} \bold{32}(5):1792-1797. DOI:
  \doi{10.1093/nar/gkh340}.

  Edgar, R. C. (2004) MUSCLE: a multiple sequence alignment method 
  with reduced time and space complexity.
  \emph{BMC Bioinformatics} \bold{5}:113. DOI:
  \doi{10.1186/1471-2105-5-113}.
}  
\seealso{\code{\link{msaClustalW}},
  \code{\link{msaClustalOmega}}, \code{\link{msaMuscle}},
  \code{\link{msaPrettyPrint}}, \code{\linkS4class{MsaAAMultipleAlignment}},
  \code{\linkS4class{MsaDNAMultipleAlignment}},
  \code{\linkS4class{MsaRNAMultipleAlignment}},
  \code{\linkS4class{MsaMetaData}}
}
\examples{
## read sequences
filepath <- system.file("examples", "exampleAA.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)

## call unified interface msa() for default method (ClustalW) and
## default parameters
msa(mySeqs)

## call ClustalOmega through unified interface
msa(mySeqs, method="ClustalOmega")

## call MUSCLE through unified interface with some custom parameters
msa(mySeqs, method="Muscle", gapOpening=12, gapExtension=3, maxiters=16,
    cluster="upgmamax", SUEFF=0.4, brenner=FALSE,
    order="input", verbose=FALSE)
}
\keyword{manip}

