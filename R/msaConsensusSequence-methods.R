msaConsensusSequence.matrix <- function(x, type=c("Biostrings", "upperlower"),
                                        thresh=c(80, 20), ignoreGaps=FALSE,
                                        ...)
{
    type <- match.arg(type)

    if (is.null(rownames(x)) ||
        any(!(rownames(x) %in% c(LETTERS, "-", "+", ".", "*"))))
        stop("consensus matrix 'x' is not in proper format")

    sel <- match(c("+", "."), rownames(x))
    sel <- sel[which(!is.na(sel))]
    if (length(sel) > 0)
        x <- x[-sel, ]

    if (!is.numeric(thresh) || length(thresh) !=2 ||
        thresh[1] > 100 || thresh[2] < 0 || thresh[1] < thresh[2])
        stop("'thresh' must be a decreasing sequence of two numbers ",
             "between 0 and 100")

    thresh <- thresh / 100

    if (type == "Biostrings")
    {
        cs <- colSums(x)

        if (any(is.na(cs)))
        {
            res <- rep.int("#", ncol(x))

            sel <- which(!is.na(cs))

            if (length(sel) > 0)
            {
                sstr <- consensusString(x[, sel, drop=FALSE], ...)
                res[sel] <- unlist(strsplit(sstr, ""))
            }

            out <- paste(res, collapse="")
        }
        else
            out <- consensusString(x, ...)
    }
    else
    {
        if (ignoreGaps)
        {
            sel <- match("-", rownames(x))

            if (!is.na(sel))
                sel <- (1:nrow(x))[-sel]
            else
                sel <- (1:nrow(x))

            perColFunc <- function(y)
            {
                if (any(is.na(y)))
                    char <- "#"
                else
                {
                    y <- y[sel] / length(sel)

                    maxw <-which.max(y[sel])
                    maxi <- y[sel[maxw]]

                    char <- rownames(x)[sel[maxw]]

                    if (maxi < thresh[1])
                    {
                        if (maxi >= thresh[2])
                            char <- tolower(char)
                        else
                            char <- "."
                    }
                }

                char
            }
        }
        else
        {
            perColFunc <- function(y)
            {
                if (any(is.na(y)))
                    char <- "#"
                else
                {
                    y <- y / length(y)

                    maxw <- which.max(y)
                    maxi <- y[maxw]

                    char <- rownames(x)[maxw]

                    if (maxi < thresh[1])
                    {
                        if (maxi >= thresh[2])
                            char <- tolower(char)
                        else
                            char <- "."
                    }
                }

                char
            }
        }

        out <- paste(apply(x, 2, perColFunc), collapse="")
    }

    out
}

setMethod("msaConsensusSequence", signature("matrix"), msaConsensusSequence.matrix)


setMethod("msaConsensusSequence", signature("MultipleAlignment"),
          function(x, ...)
          {
              mat <- consensusMatrix(x)
              msaConsensusSequence(mat, ...)
          })
