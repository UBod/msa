# msa: An R Package for Multiple Sequence Alignment
The 'msa' package provides a unified R/Bioconductor interface to
the multiple sequence alignment algorithms ClustalW, ClustalOmega,
and Muscle. All three algorithms are integrated in the package,
therefore, they do not depend on any external software tools
and are available for all major platforms. The multiple sequence
alignment algorithms are complemented by a function for
pretty-printing multiple sequence alignments using the LaTeX
package TeXshade.

Although the package is maintained by Ulrich Bodenhofer, the package itself
has been implemented mainly by Enrico Bonatesta and Christoph Kainrath
(formerly Christoph Horejs-Kainrath).

## Installation

The package can be installed from
[Bioconductor](https://bioconductor.org/). Therefore, the the simplest way to install the package is to enter
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("msa")
```
into your R session. If, for what reason ever, you prefer to install the package manually, follow the instructions in the [user manual](https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf).

## User support

If you encounter any issues or if you have any question that might be of interest also for other users, before writing a private message to the package developers/maintainers, please create an issue in this repository and also consider posting on [Bioconductor Support](https://support.bioconductor.org/) or on [StackOverflow](https://stackoverflow.com/). For other matters regarding the package, please contact the package author.


## Citing this package

If you use this package for research that is published later, you are kindly asked to cite it as follows:

- U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter (2015). msa: an R package for multiple sequence alignment. *Bioinformatics* **31**(24):3997-3999. DOI: [10.1093/bioinformatics/btv494](http://doi.org/10.1093/bioinformatics/btv494).

Moreover, we insist that, any time you use/cite the package, you also cite the original paper in which the algorithm/method/package that you have been using has been introduced:

**ClustalW:**

- J. D. Thompson, D. G. Higgins, and T. J. Gibson (1994). CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice. *Nucleic Acids Res.* **22**(22):4673 4680. DOI: [10.1093/nar/22.22.4673](http://doi.org/10.1093/nar/22.22.4673).

**ClustalOmega:**

- F. Sievers, A. Wilm, D. Dineen, T. J. Gibson, K. Karplus, W. Li, R. Lopez, H. McWilliam, M. Remmert, J. Söding, J. D. Thompson, and D. G. Higgins (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. *Mol. Syst. Biol.* **7**:539. DOI: [10.1038/msb.2011.75](http://doi.org/10.1038/msb.2011.75).

**MUSCLE:**

- R. C. Edgar (2004). MUSCLE: a multiple sequence alignment method with reduced time and space complexity. *BMC Bioinformatics* **5**(5):113. DOI: [10.1186/1471-2105-5-113](http://doi.org/10.1186/1471-2105-5-113).
- R. C. Edgar (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. *Nucleic Acids Res.* **32**(5):1792 1797. DOI: [10.1093/nar/gkh340](http://doi.org/10.1093/nar/gkh340).

**TeXshade:**

- E. Beitz (2000). TeXshade: shading and labeling of multiple sequence alignments using LaTeX2e. *Bioinformatics* **16**(2):135-139. DOI: [10.1093/bioinformatics/16.2.135](http://doi.org/10.1093/bioinformatics/16.2.135).

