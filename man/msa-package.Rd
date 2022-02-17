\name{msa-package}
\alias{msaPackage}
\alias{msa-package}
\docType{package}
\title{Multiple Sequence Alignment}
\description{
  The \pkg{msa} package provides a unified R/Bioconductor interface to
  different multiple sequence alignment algorithms.
  Currently, \sQuote{ClustalW}, \sQuote{ClustalOmega}, and
  \sQuote{MUSCLE} are supported. All algorithms are usable without
  additional software packages and on all major platforms.
  The multiple sequence algorithms are complemented by an R interface to
  the powerful LaTeX package \pkg{texshade.sty} which allows for a
  highly customizable plots of multiple sequence alignments.
}
\details{
  \tabular{ll}{
    Package: \tab msa\cr
    Type: \tab Package\cr
    Version: \tab 1.1.2\cr
    Date: \tab 2015-09-29\cr
    License: \tab GPL-2\cr
  }
}
\author{Enrico Bonatesta, Christoph Horejs-Kainrath, and Ulrich Bodenhofer
  <msa@bioinf.jku.at>
}
\references{
  \url{http://www.bioinf.jku.at/software/msa}
  
  U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter
  (2015). msa: an R package for multiple sequence alignment. 
  \emph{Bioinformatics} \bold{31}(24):3997-3999. DOI:
  \doi{10.1093/bioinformatics/btv494}.

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

  Beitz, E. (2000) TeXshade: shading and labeling of multiple
  sequence alignments using LaTeX2e
  \emph{Bioinformatics} \bold{16}(2):135-139. DOI:
  \doi{10.1093/bioinformatics/16.2.135}.
}  
\seealso{\code{\link{msa}},  \code{\link{msaClustalW}},
  \code{\link{msaClustalOmega}}, \code{\link{msaMuscle}},
  \code{\link{msaPrettyPrint}}
}
\examples{
## read sequences
filepath <- system.file("examples", "exampleAA.fasta", package="msa")
mySeqs <- readAAStringSet(filepath)

## call unified interface msa() for default method (ClustalW) and
## default parameters
msa(mySeqs)
}
\keyword{package}
