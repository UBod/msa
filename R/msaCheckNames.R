msaCheckNames <- function(x, replacement=" ", verbose=TRUE)
{
    if (!is(x, "MultipleAlignment"))
        stop("x must be a multiple alignment object")

    out <- x

    pattern <- "[^a-zA-Z0-9,;:.?!/\\-\\(\\)\\'\" ]"

    if (length(grep(pattern, rownames(x), perl=TRUE)) > 0)
    {
        if (verbose)
            message("sequence names contain invalid characters")

        rownames(out) <- gsub(pattern, replacement, rownames(x),
                              perl=TRUE)
    }

    invisible(out)
}
