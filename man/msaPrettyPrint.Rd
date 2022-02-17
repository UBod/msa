\name{msaPrettyPrint}
\alias{msaPrettyPrint}
\title{Pretty-Printing of Multiple Sequence Alignments}
\description{
  The \code{msaPrettyPrint} function provides an R interface to the
  powerful LaTeX package \pkg{texshade.sty} which allows for a
  highly customizable plots of multiple sequence alignments.}
\usage{
    msaPrettyPrint(x, y, output=c("pdf", "tex", "dvi", "asis"),
                   subset=NULL, file=NULL, alFile=NULL,
                   askForOverwrite=TRUE,  psFonts=FALSE, code=NA,
                   paperWidth=11, paperHeight=8.5, margins=c(0.1, 0.3),
                   shadingMode=c("identical", "similar", "functional"),
                   shadingModeArg=NA,
                   shadingColors=c("blues", "reds", "greens", "grays",
                                   "black"),
                   showConsensus=c("bottom", "top", "none"),
                   consensusColors=c("ColdHot", "HotCold", "BlueRed",
                                     "RedBlue", "GreenRed",
                                     "RedGreen", "Gray"),
                   consensusThreshold=50,
                   showLogo=c("top", "bottom", "none"),
                   logoColors=c("chemical", "rasmol", "hydropathy",
                                "structure", "standard area",
                                "accessible area"),
                   showLogoScale=c("none", "leftright", "left",
                                   "right"),
                   showNames=c("left", "right", "none"),
                   showNumbering=c("right", "left", "none"),
                   showLegend=TRUE, furtherCode=NA, verbose=FALSE)
}
\arguments{
  \item{x}{an object of class \code{\linkS4class{MultipleAlignment}},
    which includes the classes \code{\linkS4class{MsaAAMultipleAlignment}}, 
    \code{\linkS4class{MsaDNAMultipleAlignment}}, and
    \code{\linkS4class{MsaRNAMultipleAlignment}}.}
  \item{y}{argument for restricting the output to a subset of columns;
    can be a numeric vector of length 2 with a lower and an upper bound
    or an object of class \code{\linkS4class{IRanges}}. If missing,
    the entire multiple alignment is printed.}
  \item{output}{type of output to be generated (see details below)}
  \item{subset}{can be used to specify a subset of sequences in the
    multiple alignment \code{x} if not all sequences should be printed.}
  \item{file}{name of output file; if no name 
    is given, the name of the output file defaults to name of the
    object provided as argument \code{x} along with the proper suffix
    which depends on the type of output specified with the \code{output}
    argument. Note that this might lead to
    invalid file names if not the name of an object, but an R expression
    is passed as argument \code{x}.}
  \item{alFile}{name of alignment file to be created;
    \code{msaPrettyPrint} first writes the multiple alignment \code{x}
    to a \code{.fasta} file. The name of this file can be
    determined with the \code{alFile} argument. If no name 
    is given, the name of the output file defaults to name of the
    object provided as argument \code{x} along with the suffix
    \code{.fasta}. Note that this might lead to
    invalid file names if not the name of an object, but an R expression
    is passed as argument \code{x}.}
  \item{askForOverwrite}{if \code{TRUE} (default), 
    \code{msaPrettyPrint} asks whether existing files should be 
    overwritten or not. If \code{askForOverwrite} is set to \code{FALSE}, 
    files are overwritten without further notice.}
  \item{psFonts}{if \code{TRUE}, \code{msaPrettyPrint} produces LaTeX
    code that includes the LaTeX package \pkg{times.sty}; if \code{FALSE}, 
    \code{msaPrettyPrint} produces LaTeX
    code based on the standard LaTeX fonts (default). Ignored for
    \code{output="asis"}.}
  \item{code}{this argument can be used to specify the entire LaTeX code
    in the \code{texshade} environment. This overrides all arguments
    that customize the appearance of the output. Instead, all
    customizations must be done as LaTeX commands provided by the
    package \pkg{texshade.sty} directly. This option should only be used
    by expert users and for special applciations in which the
    possibilities of the customizations of the \code{msaPrettyPrint}
    function turn out to be insufficient.}
  \item{paperWidth,paperHeight}{paper format to be used in the resulting
    document; defaults to 11in x 8.5in (US letter in landscape
    orientation). Ignored for \code{output="asis"}.}
  \item{margins}{a numeric vector of length 2 with the horizontal and
    vertical margins, respectively; the default is 0.1in for the
    horizontal and 0.3in for the vertical margin.}
  \item{shadingMode}{shading mode; currently the shading modes
    \code{"identical"}, \code{"similar"}, and \code{"functional"}
    are supported (see documentation of
    \pkg{texshade.sty} for details).}
  \item{shadingModeArg}{for shading modes \code{"identical"} and
    \code{"similar"}, \code{shadingModeArg} must 
    be a single numeric threshold between 0 and 100 or two thresholds 
    between 0 and 100 in increasing order. For shading mode
    \code{"functional"}, valid \code{shadingModeArg} arguments are
    \code{"charge"}, \code{"hydropathy"}, \code{"structure"}, 
    \code{"chemical"}, \code{"rasmol"}, 
    \code{"standard area"}, and \code{"accessible area"} (see documentation of
    \pkg{texshade.sty} for details).}
  \item{shadingColors}{color scheme for shading; valid
    \code{"shadingColors"} arguments are \code{"blues"}, \code{"reds"},
    \code{"greens"}, \code{"grays"}, and \code{"black"} (see documentation of
    \pkg{texshade.sty} for details).}
  \item{showConsensus}{where to show the consensus sequence; 
     possible values are \code{"bottom"}, \code{"top"}, and
     \code{"none"} (the latter option suppresses printing of the
     consensus sequence).}
   \item{consensusColors}{color scheme for printing the consensus
     sequence; the following choices are possible: \code{"ColdHot"},
     \code{"HotCold"}, \code{"BlueRed"}, \code{"RedBlue"}, \code{"GreenRed"},
     \code{"RedGreen"}, and \code{"Gray"} (see documentation of
     \pkg{texshade.sty} for details).}
   \item{consensusThreshold}{one or two numbers between 0 and 100, where
     the second one is optional and must be larger than the
     first one (see documentation of \pkg{texshade.sty} for details)}
   \item{showLogo}{where to show a sequence logo; 
     possible values are \code{"top"}, \code{"bottom"}, or 
     \code{"none"} (the latter option suppresses printing of the
     consensus sequence). If a sequence logo and a consensus
     sequence should be shown together, they can only be located at
     opposite sides.}
   \item{logoColors}{color scheme for printing the sequence logo;
     the following choices are possible: \code{"chemical"}, 
     \code{"rasmol"}, \code{"hydropathy"}, \code{"structure"},
     \code{"standard area"}, and \code{"accessible area"} (see
     documentation of \pkg{texshade.sty} for details).}
   \item{showLogoScale}{where to plot the vertical axis of the
     sequence logo; possible values are \code{"left"}, \code{"right"},
     \code{"leftright"}, and \code{"none"} (the latter option suppresses
     that the axis is displayed).}
    \item{showNames}{where to print sequence names;
      possible values are \code{"left"}, \code{"right"}, and
      \code{"none"} (the latter option suppresses that names are displayed).}
    \item{showNumbering}{where to print sequence numbers;
      possible values are \code{"left"}, \code{"right"}, and
      \code{"none"} (the latter option suppresses that numbers are displayed).
      If sequence names and numbers should be shown together,
      they can only be located at opposite sides.}
    \item{showLegend}{if \code{TRUE} (default), a legend is 
        printed at the end of the alignment.}
      \item{furtherCode}{additional LaTeX code to be included in the
	\code{texshade} environment; all text passed as
	\code{furtherCode} is placed between the commands created by
	\code{msaPrettyPrint} and the end of the \code{texshade} environment.
        Note the difference to the \code{code} argument: while the
	\code{code} argument replaces all LaTeX code in the
	\code{texshade} environment, the code passed as
	\code{furtherCode} argument is added to the LaTeX code in the
	\code{texshade} environment.}
      \item{verbose}{if \code{TRUE} (default), progress messages are
	printed and also the output of running (PDF)LaTeX (if
	applicable) is printed to the R session.}
}
\details{The \code{msaPrettyPrint} function writes a multiple alignment
  to a \code{.fasta} file and creates LaTeX code for
  pretty-printing the multiple alignment on the basis of the
  LaTeX package \pkg{texshade.sty}. If \code{output="asis"},
  \code{msaPrettyPrint} prints a LaTeX fragment consisting of the
  \code{texshade} environment to the console. The parameters described
  above can be used to customize the way the multiple alignment is
  formatted. If \code{output="tex"}, a complete LaTeX file including
  preamble is created. For \code{output="dvi"} and \code{output="pdf"},
  the same kind of LaTeX file is created, but processed using (PDF)LaTeX
  to produce a final DVI or PDF file, respectively. The \code{file}
  argument be used to determine the file name of the final output file
  (except for the \code{output="asis"} which does not create an output
  file).

  The choice \code{output="asis"} is particularly useful for Sweave or
  knitr documents. If \code{msaPrettyPrint} is called with
  \code{output="asis"} in a code chunk with \code{results="tex"}
  (Sweave) or \code{results="asis"} (knitr), then the resulting LaTeX
  fragment consisting of the \code{texshade} environment is directly
  included in the LaTeX document that is created from the Sweave/knitr
  document.

  As noted above, if they are not specified explicitly, output file
  names are determined automatically. It is important to point out that
  all file names need to be LaTeX-compliant, i.e. no special characters
  and spaces (!) are allowed. If a file name would be invalid,
  \code{msaPrettyPrint} makes a default choice.

  Moreover, if sequence names are to be printed, there might be names
  that are not LaTeX-compliant and lead to LaTeX errors. In order
  to check that in advance, the function \code{\link{msaCheckNames}}
  is available.

  Note that \code{\link{texi2dvi}} and \code{\link{texi2pdf}} always
  save the resulting DVI/PDF files to the current working directory,
  even if the LaTeX source file is in a different directory.
  That is also the reason why the temporary file is created in the
  current working directory in the example below.
}
\value{\code{msaPrettyPrint} returns an invisible character vector consisting
  of the LaTeX fragment with the \code{texshade} environment.
}
\author{Ulrich Bodenhofer, Enrico Bonatesta, and Christoph Horejs-Kainrath 
  <msa@bioinf.jku.at>
}
\references{
  \url{http://www.bioinf.jku.at/software/msa}
  
  U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter
  (2015). msa: an R package for multiple sequence alignment. 
  \emph{Bioinformatics} \bold{31}(24):3997-3999. DOI:
  \doi{10.1093/bioinformatics/btv494}.

  \url{https://www.ctan.org/pkg/texshade}

  Beitz, E. (2000) TeXshade: shading and labeling of multiple
  sequence alignments using LaTeX2e
  \emph{Bioinformatics} \bold{16}(2):135-139. DOI:
  \doi{10.1093/bioinformatics/16.2.135}.
}
\seealso{\code{\link{msaCheckNames}}
}
\examples{
## read sequences
filepath <- system.file("examples", "exampleAA.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)

## call unified interface msa() for default method (ClustalW) and
## default parameters
myAlignment <- msa(mySeqs)

## show resulting LaTeX code with default settings
msaPrettyPrint(myAlignment, output="asis", askForOverwrite=FALSE)

## create PDF file according to some custom settings
\dontrun{
tmpFile <- tempfile(pattern="msa", tmpdir=".", fileext=".pdf")
tmpFile
msaPrettyPrint(myAlignment, file=tmpFile, output="pdf",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="rasmol",
               verbose=FALSE, askForOverwrite=FALSE)

library(Biobase)
openPDF(tmpFile)}
}
\keyword{graphs}