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
        verbose=FALSE, help=FALSE, ...)
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
    \code{\link{msaClustalOmega}}, or \code{\link{msaMuscle}});
    In order to standardize interfaces,
    all algorithms consistently use negative values for the gap opening
    penalty.}
  \item{gapExtension}{gap extension penalty; the defaults are
    specific to the algorithm (see \code{\link{msaClustalW}},
    \code{\link{msaClustalOmega}}, or \code{\link{msaMuscle}});
    In order to standardize interfaces,
    all algorithms consistently use negative values for the gap extension
    penalty.}
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
  sequences. Therefore, it is necessary to assign names to input
  sequences in case the user wants to preserve a one-to-one assignment
  between input sequences and the sequences in the multiple sequence
  alignment. That is why all algorithms check input sequences for
  unique names and issue warnings if the input sequences are unnamed or
  if the names are not unique.
  As noted in the description of the \code{inputSeqs} argument above,
  all functions, \code{msa()}, \code{\link{msaClustalW}},
  \code{\link{msaClustalOmega}}, and \code{\link{msaMuscle}}, also allow
  for direct reading from FASTA files. This is mainly for the reason of
  memory efficiency if the sequence data set is very large. Otherwise,
  we want to encourage users to first read the sequences into the R
  workspace. In any case, if sequences are read from a FASTA file
  directly, this is completely under the control of the respective
  algorithm and does not allow for checking whether the sequences are
  named uniquely in the FASTA file.
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
  
  \url{http://www.clustal.org/download/clustalw_help.txt}

  \url{http://www.clustal.org/omega/README}
  
  \url{http://www.drive5.com/muscle/muscle.html}
  
  Thompson, J. D., Higgins, D. G., and Gibson, T. J. (1994)
  CLUSTAL W: improving the sensitivity of progressive multiple sequence
  alignment through sequence weighting, position-specific gap penalties
  and weight matrix choice.
  \emph{Nucleic Acids Res.} \bold{22}(22):4673-4680. DOI:
  \href{http://dx.doi.org/10.1093/nar/22.22.4673}{10.1093/nar/22.22.4673}.

  Sievers, F., Wilm, A., Dineen, D., Gibson, T. J., Karplus, K., Li, W.,
  Lopez, R., McWilliam, H., Remmert, M., Soeding, J., Thompson, J. D.,
  and Higgins, D. G. (2011) Fast, scalable generation of high-quality
  protein multiple sequence alignments using Clustal Omega.
  \emph{Mol. Syst. Biol.} \bold{7}:539. DOI:
  \href{http://dx.doi.org/10.1038/msb.2011.75}{10.1038/msb.2011.75}.

  Edgar, R. C. (2004) MUSCLE: multiple sequence alignment with high 
  accuracy and high throughput.
  \emph{Nucleic Acids Res.} \bold{32}(5):1792-1797. DOI:
  \href{http://dx.doi.org/10.1093/nar/gkh340}{10.1093/nar/gkh340}.

  Edgar, R. C. (2004) MUSCLE: a multiple sequence alignment method 
  with reduced time and space complexity.
  \emph{BMC Bioinformatics} \bold{5}:113. DOI:
  \href{http://dx.doi.org/10.1186/1471-2105-5-113}{10.1186/1471-2105-5-113}.
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
msa(mySeqs, method="ClustalOmega", cluster=120, gapOpening=-3,
    gapExtension=-1, maxGuidetreeIterations=20, outFmt="clustal",
    outputOrder="input-order", substitutionMatrix="blosum30",
    type="protein", verbose=FALSE)
}
\keyword{manip}

