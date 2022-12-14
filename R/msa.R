msa <- function(inputSeqs,
                method=c("ClustalW", "ClustalOmega", "Muscle"),
                cluster="default",
                gapOpening="default",
                gapExtension="default",
                maxiters="default",
                substitutionMatrix="default",
                type="default",
                order=c("aligned", "input"),
                verbose=FALSE,
                help=FALSE,
                ...)
{
    if (help && identical(method, c("ClustalW", "ClustalOmega", "Muscle")))
    {
        file <- system.file("extdata/msaHelpPrefix.txt", package="msa")
        cat(readLines(file), sep="\n")
    }

    method <- match.arg(method)

    msaFun <- get(paste0("msa", method), envir=environment(msa))

    out <- msaFun(inputSeqs=inputSeqs,
                  cluster=cluster,
                  gapOpening=gapOpening,
                  gapExtension=gapExtension,
                  maxiters=maxiters,
                  substitutionMatrix=substitutionMatrix,
                  type=type,
                  order=order,
                  verbose=verbose,
                  help=help,
                  ...)

    if (is(out, "MsaMetaData"))
        out@call <- deparse(sys.call())

    if (is.null(out))
        invisible(out)
    else
        out
}
