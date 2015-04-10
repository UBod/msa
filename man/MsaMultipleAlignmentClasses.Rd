\name{MsaMultipleAnlignmentClasses}
\docType{class}
\alias{class:MsaAAMultipleAlignment}
\alias{MsaAAMultipleAlignment-class}
\alias{MsaAAMultipleAlignment}
\alias{class:MsaDNAMultipleAlignment}
\alias{MsaDNAMultipleAlignment-class}
\alias{MsaDNAMultipleAlignment}
\alias{class:MsaRNAMultipleAlignment}
\alias{MsaRNAMultipleAlignment-class}
\alias{MsaRNAMultipleAlignment}
\alias{show,MsaAAMultipleAlignment-method}
\alias{show,MsaDNAMultipleAlignment-method}
\alias{show,MsaRNAMultipleAlignment-method}
\alias{params,MsaAAMultipleAlignment-method}
\alias{params,MsaDNAMultipleAlignment-method}
\alias{params,MsaRNAMultipleAlignment-method}
\alias{print,MsaAAMultipleAlignment-method}
\alias{print,MsaDNAMultipleAlignment-method}
\alias{print,MsaRNAMultipleAlignment-method}
\alias{version,MsaAAMultipleAlignment-method}
\alias{version,MsaDNAMultipleAlignment-method}
\alias{version,MsaRNAMultipleAlignment-method}
\title{Classes \code{MsaAAMultipleAlignment},
  \code{MsaDNAMultipleAlignment}, and \code{MsaRNAMultipleAlignment}}
\description{S4 classes for storing multiple alignments of amino
  acid, DNA, and RNA sequences along with algorithm metadata}
\section{Objects}{
  Objects of these classes are returned by the multiple sequence
  alignment algorithms  \code{\link{msaClustalW}},
  \code{\link{msaClustalOmega}}, \code{\link{msaMuscle}}, and the
  wrapper function \code{\link{msa}}, all of which are
  provided by the \pkg{msa} package.}
\section{Details}{
  The class \code{MsaAAMultipleAlignment} extends the
  \code{\linkS4class{AAMultipleAlignment}} class, the class
  \code{MsaDNAMultipleAlignment} extends the
  \code{\linkS4class{DNAMultipleAlignment}} class, and the class
  \code{MsaRNAMultipleAlignment} extends the
  \code{\linkS4class{RNAMultipleAlignment}} class. All three classes
  extend their parent classes by the slots contained in the
  \code{\linkS4class{MsaMetaData}}, i.e. all three classes are class
  unions of the aforementioned parent classes and the class
  \code{\linkS4class{MsaMetaData}}.
}
\section{Methods}{
  \describe{
    \item{\code{print(x, show=c("alignment", "version", "call"))}:}{
      prints information about the object \code{x}; the \code{show}
      argument allows for determining what should be printed.
      The \code{show} must be a character vector and may contain any
      combination of the following five strings:
      if \code{show} contains \code{"alignment"}, the multiple
      sequence alignment is printed using the corresponding method
      from the \pkg{Biostrings} package.
      If \code{show} contains \code{"version"},
      the \code{version} slot is shown. If \code{show} contains
      \code{"call"}, the \code{call} slot is shown.
      If \code{show} contains \code{"standardParams"}, the
      settings of the parameters that are common to all three
      multiple sequence alignment algorithms are shown. If \code{show}
      contains \code{"algParams"}, the
      algorithm-specific parameters are shown. 
      The order in which the strings are placed in the \code{show}
      argument does not have an effect on the order in which
      data are printed. The default is
      \code{show=c("alignment", "version", "call")}, i.e. by default,
      the multiple sequence alignment is shown along with version and
      call information.
    }
    \item{\code{show(object)}:}{displays the alignment along with
      metadata; synonymous to calling \code{print} with default
      \code{show} argument \code{c("alignment", "version", "call")}}
    \item{\code{version(object)}:}{displays the algorithm with which
      the multiple alignment has been computed along with its
      version number (see also \code{\linkS4class{MsaMetaData}}).}
    \item{\code{params(x)}:}{accessor to the \code{params} slot (see
      also \code{\linkS4class{MsaMetaData}})}
  }
}
\author{Enrico Bonatesta and Christoph Horejs-Kainrath
  <msa@bioinf.jku.at>
}
\seealso{\code{\link{msa}}, \code{\link{msaClustalW}},
  \code{\link{msaClustalOmega}}, \code{\link{msaMuscle}},
  \code{\linkS4class{MsaMetaData}}
}
\examples{
## read sequences
filepath <- system.file("examples", "exampleAA.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)

## simple call with default values
myAlignment <- msaClustalOmega(mySeqs)

## show the algorithm version with which the results were created
version(myAlignment)

## show the results
show(myAlignment)

## print the results
print(myAlignment, show="alignment")
print(myAlignment, show=c("alignment", "version"))
print(myAlignment, show="standardParams")
print(myAlignment, show="algParams")
print(myAlignment, show=c("call", "version"))

## show the params
params(myAlignment)
}
\keyword{class}