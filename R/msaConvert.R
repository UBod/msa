msaConvert <- function(x, type=c("seqinr::alignment",
                                 "bios2mds::align",
                                 "ape::AAbin",
                                 "ape::DNAbin",
                                 "phangorn::phyDat",
                                 "bio3d::fasta"))
{
    type <- match.arg(type)

    if (!is(x, "MultipleAlignment"))
        stop("x must be a 'MultipleAlignment' object")

    xn <- as.character(unmasked(x))

    if (type == "seqinr::alignment")
    {
        out <- list(nb=length(xn),
                    nam=names(xn),
                    seq=unname(xn),
                    com=NA)

        class(out) <- "alignment"
    }
    else if (type == "bios2mds::align")
    {
        out <- .Call("SplitCharVector2List", xn)
        names(out) <- names(xn)
        class(out) <- "align"
    }
    else if (type == "ape::AAbin")
    {
        if (!is(x, "AAMultipleAlignment"))
            stop("conversion to 'ape::AAbin' only supported for ",
                 "amino acid sequences")

        out <- .Call("SplitCharVector2Matrix", xn, "X")
        rownames(out) <- names(xn)
        class(out) <- "AAbin"
    }
    else if (type == "ape::DNAbin")
    {
        if (!is(x, "DNAMultipleAlignment"))
            stop("conversion to 'ape::DNAbin' only supported for ",
                 "DNA sequences")

        out <- .Call("SplitCharVector2Matrix", tolower(xn), "n")
        rownames(out) <- names(xn)
        class(out) <- "DNAbin"
    }
    else if (type == "phangorn::phyDat")
    {
        if (!is(x, "DNAMultipleAlignment"))
            stop("conversion to 'phangorn::phyDat' only supported for ",
                 "DNA sequences")

        if (requireNamespace("phangorn", quietly=TRUE))
            out <- phangorn::as.phyDat(x)
        else
            stop("conversion to 'phyDat' requires package 'phangorn'")
    }
    else if (type == "bio3d::fasta")
    {
        out <- list(id=names(xn),
                    ali=.Call("SplitCharVector2Matrix", xn, "-"),
                    call=x@call)
        rownames(out$ali) <- out$id
        class(out) <- "fasta"
    }

    out
}
