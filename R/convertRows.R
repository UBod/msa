convertAlnRows <- function(rows, type)
{
    version <- rows[1]

    if (length(rows) < 3 ||
        ##!identical(grep("^CLUSTAL", rows[1L]), 1L) ||
        !identical(sub("^\\s+$", "", rows[2:3]), c("", "")))
        stop("There is an invalid aln file!")

    rows <- tail(rows, -3)
    rows <- sub("^(\\S+\\s+\\S+)\\s*\\d*$", "\\1", rows)

    markupPattern <- "^(\\s|\\*|:|\\.)*$"

    markupLines <- grep(markupPattern, rows, perl=TRUE)
    alnLines <- gaps(as(markupLines, "IRanges"), start=1, end=length(rows))
    nseq <- unique(width(alnLines))

    if (length(nseq) != 1)
        stop("There are missing alignment rows!")

    rows <- extractROWS(rows, alnLines)
    spaces <- regexpr("\\s+", rows)
    ids <- substr(rows, 1L, spaces - 1L)
    nsplits <- length(rows) %/% nseq

    if (!identical(ids, rep.int(head(ids, nseq), nsplits)))
        stop("The alignment rows are out of order!")

    alns <- substr(rows, spaces + attr(spaces, "match.length"), nchar(rows))

    chrs <- structure(do.call(paste,
                              c(split(alns, rep(seq_len(nsplits),
                                                each=nseq)), sep="")),
                      names=head(ids, nseq))

    type <- switch(type, dna="DNA", rna="RNA", protein="AA")

    out <- new(paste0("Msa", type, "MultipleAlignment"),
               unmasked=do.call(paste0(type, "StringSet"), list(chrs)),
               rowmask=as(IRanges(), "NormalIRanges"),
               colmask=as(IRanges(), "NormalIRanges"),
               version=version)
}
