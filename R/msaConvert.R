msaConvert <- function(x, type=c("seqinr::alignment",
                                 "bios2mds::align"))
{
    type <- match.arg(type)

    if (!is(x, "MultipleAlignment"))
        stop("x must be a 'MultipleAlignment' object")

    x <- as.character(unmasked(x))

    if (type == "seqinr::alignment")
    {
        out <- list(nb=length(x),
                    nam=names(x),
                    seq=unname(x),
                    com=NA)

        class(out) <- "alignment"
    }
    else if (type == "bios2mds::align")
    {
        out <- .Call("SplitCharVector2List", x)
        names(out) <- names(x)
        class(out) <- "align"
    }

    out
}
