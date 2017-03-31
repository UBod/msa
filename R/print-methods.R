print.MsaMetaData <- function(x, show)
{
    sep <- FALSE

    if ("version" %in% show)
    {
        cat( x@version, "\n")
        sep <- TRUE
    }

    stdParams <- c("gapOpening", "gapExtension", "maxiters",
                   "verbose")
    algParams <- setdiff(names(x@params), stdParams)

    if ("standardParams" %in% show)
    {
        if (sep) cat("\n")

        cat("Standard params:\n")
        cat(paste0("   ", stdParams, " = ",
            as.character(sapply(stdParams, function(p) x@params[[p]]))
            ),
            sep="\n", fill=FALSE)

        sep <- TRUE
    }

    if ("algParams" %in% show)
    {
        if (sep) cat("\n")

        cat("Options specific to ", x@version, ":\n", sep="")
        cat(paste0("   ", algParams, " = ",
                  as.character(sapply(algParams, function(p) x@params[[p]]))),
            sep="\n", fill=FALSE)

        sep <- TRUE
    }

    if ("call" %in% show)
    {
        if (sep) cat("\n")

        cat("Call:\n   ", x@call, "\n", sep="")
    }
}

print.MsaMultipleAlignmentChunk <- function(str, names=NULL, halfNrow=9, pos="",
                                            addOne=FALSE)
{
    lx <- length(str)
    iW <- nchar(as.character(lx - if (addOne) 1 else 0)) + 2

    cat(format("", width=iW), sep = "")

    if (length(names) > 0)
        cat(format(paste(" aln", pos), width=nchar(str)[1]), " names\n")
    else
        cat(paste(" aln", pos), "\n", sep="")

    if (addOne)
        str <- paste0(format(c(paste("[", 1:(lx - 1), "]", sep=""), "Con"),
                             width=iW, justify="right"), " ",
                      str, if (is.null(names)) "" else paste0(" ", names))
    else
        str <- paste0(format(paste("[", 1:lx, "]", sep=""), width=iW,
                             justify="right"), " ",
                      str, if (is.null(names)) "" else paste0(" ", names))

    if (lx <= 2 * halfNrow + 1)
        cat(paste(str, collapse="\n"), "\n")
    else
    {
        cat(paste(str[1:halfNrow], collapse="\n"), "\n")
        cat(format("...", width = iW, justify = "right"), "...\n")
        cat(paste(str[(lx - halfNrow + 1):lx], collapse="\n"), "\n")
    }
}

print.MsaMultipleAlignment <- function(x, show=c("alignment", "version", "call"),
                                       showNames=TRUE, showConsensus=TRUE,
                                       halfNrow=9, nameWidth=20, ...)
{
    show <- match.arg(show,
                      choices=c("alignment", "complete", "version", "call",
                                "standardParams", "algParams", "all"),
                      several.ok=TRUE)

    if ("all" %in% show)
        show <- c("complete", "version", "call", "standardParams", "algParams")

    print.MsaMetaData(x, show=show)

    if (any(c("alignment", "complete") %in% show))
    {
        nr <- nrow(x)
        nc <- ncol(x)

        if (identical(halfNrow, NA) || identical(halfNrow, -1))
            halfNrow <- nr

        if (!is.numeric(halfNrow) || length(halfNrow) != 1 ||
             round(halfNrow) != halfNrow || halfNrow < 1)
            stop("halfNrow must be a single whole number or NA")

        if (!is.numeric(nameWidth) || length(nameWidth) != 1 ||
             round(nameWidth) != nameWidth || nameWidth < 5)
            stop("nameWidth must be a single whole number at least as large as 5")

        if (nameWidth > getOption("width") - 20)
        {
            nameWidth <- getOption("width") - 20
            warning("nameWidth must be at least width - 20")
        }

        cat("\n", class(x), " with ", nr,
            ifelse(nr == 1, " row and ", " rows and "),
            nc, ifelse(nc == 1, " column\n", " columns\n"), sep = "")

        if (nr > 0)
        {
            strings <- unmasked(x)
            mdim <- maskeddim(x)

            if (sum(mdim) > 0)
            {
                if (mdim[1] > 0)
                {
                    strings <- BStringSet(strings)

                    maskStrings <- rep(BStringSet(paste(rep.int("#", nc),
                                                        collapse = "")),
                                       mdim[1])

                    i <- as.integer(rowmask(x))

                    if (!is.null(rownames(x)))
                        names(maskStrings) <- rownames(x)[i]
                    strings[i] <- maskStrings
                }

                if (mdim[2] > 0)
                {
                    strings <- as.matrix(strings)
                    strings[, as.integer(colmask(x))] <- "#"
                    strings <- BStringSet(apply(strings, 1, paste,
                                                collapse = ""))
                }
            }

            strings <- as.character(strings)

            iw <- nchar(as.character(length(strings))) + 2

            if (showConsensus)
            {
                cons <- msaConsensusSequence(x, ...)
                strings <- c(strings, cons)

                if (length(names(strings)) > 0)
                    names(strings)[length(strings)] <- "Consensus"

                addOne <- TRUE
            }
            else
                addOne <- FALSE

            names <- names(strings)

            if (showNames && length(names) > 0)
            {
                names <- names(strings)

                names <- ifelse(nchar(names) <= nameWidth,
                                names,
                                paste0(substr(names, 1, nameWidth - 3), "..."))

                chunkSize <- getOption("width") - iw - nameWidth - 3
            }
            else
            {
                names <- NULL
                chunkSize <- getOption("width") - iw - 2
            }

            seqLen <- nchar(strings)[1]

            if (seqLen < 7)
                strings <- format(strings, width=7, justify="left")

            if (nchar(strings)[1] <= chunkSize)
                print.MsaMultipleAlignmentChunk(strings, names, halfNrow=halfNrow,
                                                addOne=addOne)
            else if ("complete" %in% show)
            {
                starts <- seq(from=1, to=seqLen, by=chunkSize)
                stops <- pmin(starts + chunkSize - 1, seqLen)

                n <- length(starts)

                for (i in 1:n)
                {
                    aln <- substr(strings, starts[i], stops[i])

                    pos <- paste0("(", starts[i], "..", stops[i], ")")

                    if (nchar(pos) + 4 > chunkSize)
                        pos <- paste0(substr(pos, 1, chunkSize - 7), "...")

                    if (nchar(aln)[1] < nchar(pos) + 4)
                        aln <- format(aln, width=nchar(pos) + 4, justify="left")

                    print.MsaMultipleAlignmentChunk(aln, names, halfNrow=halfNrow,
                                                    pos=pos, addOne=addOne)

                    if (i < n)
                        cat("\n")
                }
            }
           else
            {
                w1 <- (chunkSize - 2) %/% 2
                w2 <- (chunkSize - 3) %/% 2

                strings <- paste0(substr(strings, start=1, stop=w1),
                                  "...",
                                  substr(strings, start=seqLen - w2 + 1, stop=seqLen))

                print.MsaMultipleAlignmentChunk(strings, names, halfNrow=halfNrow,
                                                addOne=addOne)
            }
        }
    }
}

setMethod("print", signature("MsaDNAMultipleAlignment"),
          print.MsaMultipleAlignment)
setMethod("print", signature("MsaRNAMultipleAlignment"),
          print.MsaMultipleAlignment)
setMethod("print", signature("MsaAAMultipleAlignment"),
          print.MsaMultipleAlignment)
