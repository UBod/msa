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


print.MsaMultipleAlignment <- function(x, show=c("alignment", "version",
                                                 "call"))
{
    show <- match.arg(show,
                      choices=c("alignment", "version", "call",
                                "standardParams", "algParams"),
                      several.ok=TRUE)

    print.MsaMetaData(x, show=show)

    if ("alignment" %in% show)
    {
        cat("\n")
        print(as(x, substr(class(x), 4, nchar(class(x)))))
    }
}

setMethod("print", signature("MsaDNAMultipleAlignment"),
          print.MsaMultipleAlignment)
setMethod("print", signature("MsaRNAMultipleAlignment"),
          print.MsaMultipleAlignment)
setMethod("print", signature("MsaAAMultipleAlignment"),
          print.MsaMultipleAlignment)
