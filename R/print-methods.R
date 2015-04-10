print.MsaMetaData <- function(x, show=c("version", "standardParams",
                                        "algParams", "call"))
{
    show <- match.arg(show, several.ok=TRUE)

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

setMethod("print", signature("MsaMetaData"), print.MsaMetaData)


print.MsaMultipleAlignment <- function(x, show=c("version", "call"))
{
    print(as(x, "MsaMetaData"), show=show)
    cat("\n")
    print(as(x, substr(class(x), 4, nchar(class(x)))))
}

setMethod("print", signature("MsaDNAMultipleAlignment"),
          print.MsaMultipleAlignment)
setMethod("print", signature("MsaRNAMultipleAlignment"),
          print.MsaMultipleAlignment)
setMethod("print", signature("MsaAAMultipleAlignment"),
          print.MsaMultipleAlignment)
