msaConservationScore.matrix <- function(x, substitutionMatrix, gapVsGap=NULL,
                                        ...)
{
    if (!is.matrix(substitutionMatrix) ||
        nrow(substitutionMatrix) != ncol(substitutionMatrix) ||
        any(!(rownames(substitutionMatrix) %in% c(LETTERS, "-", "*"))) ||
        any(!(colnames(substitutionMatrix) %in% c(LETTERS, "-", "*"))) ||
        any(rownames(substitutionMatrix) != colnames(substitutionMatrix)))
        stop("substitution matrix is not in proper format")

    if (is.null(rownames(x)) ||
        any(!(rownames(x) %in% c(LETTERS, "-", "+", ".", "*"))))
        stop("consensus matrix 'x' is not in proper format")

    sel <- match(c("+", ".", "*"), rownames(x))
    sel <- sel[which(!is.na(sel))]
    if (length(sel) > 0)
        x <- x[-sel, ]

    sel <- which(rowSums(x, na.rm=TRUE) <= 0)

    if (length(sel) > 0)
        x <- x[-sel, ]

    sel <- match(c("*", "-"), rownames(substitutionMatrix))

    if (is.na(sel[2]) && !is.na(sel[1]))
    {
        rownames(substitutionMatrix)[sel[1]] <- "-"
        colnames(substitutionMatrix)[sel[1]] <- "-"
    }

    sel <- match(rownames(x), rownames(substitutionMatrix))

    if (any(is.na(sel)))
        stop("some letters occurring in alignment 'x' are missing ",
             "in substitution matrix")

    substitutionMatrix <- substitutionMatrix[sel, sel]

    if (!is.null(gapVsGap))
    {
        if (!is.numeric(gapVsGap) || length(gapVsGap) != 1)
            stop("'gapVsGap' must be NULL or a single number")

        substitutionMatrix["-", "-"] <- gapVsGap
    }

    out <- drop(apply(x, 2, function(y) crossprod(y, substitutionMatrix %*% y)))

    names(out) <- unlist(strsplit(msaConsensusSequence(x, ...), ""))

    out
}

setMethod("msaConservationScore", signature("matrix"),
          msaConservationScore.matrix)


setMethod("msaConservationScore", signature("MultipleAlignment"),
          function(x, ...)
          {
              mat <- consensusMatrix(x)
              msaConservationScore.matrix(mat, ...)
          })
