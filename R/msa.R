msa <- function(inputSeqs,
                method=c("ClustalW", "ClustalOmega", "Muscle"),
                cluster="default",
                gapOpening="default",
                gapExtension="default",
                maxiters="default",
                substitutionMatrix="default",
                type="default",
                verbose=FALSE,
                help=FALSE,
                ...)
{
    method <- match.arg(method)

    msaArgs <- as.list(match.call(expand.dots=TRUE)[-1])
    msaArgs[["method"]] <- NULL

    out <- do.call(paste0("msa", method), msaArgs)

    if (is(out, "MsaMetaData"))
        out@call <- deparse(sys.call())

    out
}
