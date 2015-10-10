msaClustalOmega <- function(inputSeqs,
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
    ##########
    ## help ##
    ##########
    if (help)
    {
        printHelp("ClustalOmega")
        return(invisible(NULL))
    }

    if (!checkFunctionAvailable("ClustalOmega"))
        stop("ClustalOmega is not available via msa!")

    params <- list(...)
    ##create a copy of the parameter list which is only used to
    ##avoid params which are not checked:
    ##after every check of a parameter, the parameter is deleted in the copy(!)
    ##if all parameter were checked, the copy finally should be empty...
    ##if not, there are additional parameters
    paramsCopy <- lapply(params, function(x) TRUE)

    ########################
    ########################
    ##  Common Parameters ##
    ########################
    ########################

    #############
    # inputSeqs #
    #############
    ##check if input of inputSeqs has a valid form, stop if not;
    ##if inputSeq is a file name set Flag TRUE else FALSE
    params[["inputSeqIsFileFlag"]] <- checkInputSeq(inputSeqs)

    ########
    # type #
    ########
    ##Sequence type

    ##validation of type
    type <- checkType(type, inputSeqs, "msaClustalOmega")

    #############
    # inputSeqs #
    #############
    inputSeqs <- transformInputSeq(inputSeqs)

    #############
    # order     #
    #############
    order <- match.arg(order)
    params[["outputOrder"]] <- switch(order,
                                      aligned="tree-order",
                                      input="input-order")

    ###########
    # cluster #
    ###########
    ##soft maximum of sequences in sub-clusters
    ##default: cluster-size = 100

    ##set default values
    if(identical(cluster, "default") || is.null(cluster)) {
        cluster <- 100
    }

    ##check, if clusterSize is a positive integer
    if (length(cluster) != 1) {
        stop("The parameter cluster should be a single positive integer!")
    }

    if (!is.integer(cluster)) {
        if (is.numeric(cluster)) {
            ##stop if usage of floats
            if (cluster - round(cluster) != 0) {
                stop("The parameter cluster should be a positive integer!")
            }
            ##stop if using negative integers
            if (cluster < 0) {
                stop("The parameter cluster should be a positive integer!")
            }
            ##stop if using numbers bigger than .Machine$integer.max
            if (cluster > .Machine$integer.max) {
                stop("The parameter cluster is bigger than an integer!")
            }
            ##typecast
            cluster <- as.integer(cluster)
        } else {
            stop("The parameter cluster should be a positive integer!")
        }
    }

    ######################
    # substitutionMatrix #
    ######################
    ##default-value: Gonnet

    ##check whether value is BlosumN or Gonnet
    if (is.null(substitutionMatrix) ||
            identical(substitutionMatrix, "default")) {
            substitutionMatrix <- NULL
    } else {
            possibleValues <- c("BLOSUM30", "BLOSUM40", "BLOSUM50",
                                "BLOSUM65", "BLOSUM80", "Gonnet")
            if (!is.character(substitutionMatrix) ||
                 !(substitutionMatrix %in% possibleValues)){
                ##create a string with all possible Values named text
                text <- ""
                text <- paste(possibleValues, collapse=", ")
                stop("The parameter substitutionMatrix",
                     "only can have the values: \n", text)
        }
    }

    ##############
    # gapOpening #
    ##############
    if (!identical(gapOpening, "default"))
        warning("msaClustalOmega currently does not support to set\n",
                "gapOpening to a non-default value!\n")

    ################
    # gapExtension #
    ################
    if (!identical(gapExtension, "default"))
        warning("msaClustalOmega currently does not support to set\n",
                "gapExtension to a non-default value!\n")

    ############
    # maxiters #
    ############
    ##Number of (combined guide-tree/HMM) iterations
    ##default: maxiter (former: "iterations") = 0

    maxiters <- checkMaxiters(maxiters, 0, "msaClustalOmega")

    ###########
    # verbose #
    ###########
    ##Write parameter settings and progress messages to log file

    verbose <- checkVerbose(FALSE, verbose)

    ######################################
    ######################################
    ######################################
    ###  ALGORITHM SPECIFIC PARAMETERS ###
    ######################################
    ######################################
    ######################################

    ########
    # auto #
    ########
    ##set options automatically (might overwrite some of your options)
    ##default: auto=FALSE

    params[["auto"]] <- checkLogicalParams("auto", params, FALSE)

    ##delete param in copy
    paramsCopy[["auto"]] <- NULL

    #################
    # clusteringOut #
    #################
    ##Clustering Output File
    ##default: clusteringOut=NULL

    if(!is.null(params[["clusteringOut"]])) {
        tempList <- checkOutFile("clusteringOut", params)
        if (tempList$existingFile) {
            interactiveCheck("clusteringOut", params)
            params[["force"]] <- TRUE
        }
        params[["clusteringOut"]] <- tempList$param
    }

    ##delete param in copy
    paramsCopy[["clusteringOut"]] <- NULL

    ###########
    # dealign #
    ###########
    ##dealign input sequences
    ##default: dealign=FALSE

    params[["dealign"]] <- checkLogicalParams("dealign", params, FALSE)

    ##delete param in copy
    paramsCopy[["dealign"]] <- NULL

    #############
    # distMatIn #
    #############
    ##Pairwise distance matrix input file
    ##default: distMatIn=NULL

    if (!is.null(params[["distMatIn"]])) {
        checkInFile("distMatIn", params)
    }

    ##delete param in copy
    paramsCopy[["distMatIn"]] <- NULL

    ##############
    # distMatOut #
    ##############
    ##Pairwise distance matrix output file
    ##default: distMatOut=NULL

    if(!is.null(params[["distMatOut"]])) {
        tempList <- checkOutFile("distMatOut", params)
        if (tempList$existingFile) {
            interactiveCheck("distMatOut", params)
            params[["force"]] <- TRUE
        }
        params[["distMatOut"]] <- tempList$param
    }

    ##delete param in copy
    paramsCopy[["distMatOut"]] <- NULL

    #########
    # force #
    #########
    ##force file overwriting
    ##default: force=FALSE

    params[["force"]] <- checkLogicalParams("force", params, FALSE)

    ##delete param in copy
    paramsCopy[["force"]] <- NULL

    ########
    # full #
    ########
    ##Use full distance matrix for guide-tree calculation
    ##(might be slow; mBed is default)
    ##default: full=FALSE

    params[["full"]] <- checkLogicalParams("full", params, FALSE)

    ##delete param in copy
    paramsCopy[["full"]] <- NULL

    #############
    # fullIter #
    #############
    ##Use full distance matrix for guide-tree calculation during iteration
    ##(might be slowish; mBed is default)
    ##default: fullIter=FALSE

    params[["fullIter"]] <- checkLogicalParams("fullIter", params, FALSE)

    ##delete param in copy
    paramsCopy[["fullIter"]] <- NULL

    ###############
    # guideTreeIn #
    ###############
    ##Guide tree input file (skips distance computation and
    ##guide tree clustering step
    ##default: distMatIn=NULL

    if (!is.null(params[["guideTreeIn"]])) {
        checkInFile("guideTreeIn", params)
    }

    ##delete param in copy
    paramsCopy[["guideTreeIn"]] <- NULL


    ################
    # guideTreeOut #
    ################
    ##Guide tree output file
    ##default: log=NULL

    if(!is.null(params[["guideTreeOut"]])) {
        tempList <- checkOutFile("guideTreeOut", params)
        if (tempList$existingFile) {
            interactiveCheck("guideTreeOut", params)
            params[["force"]] <- TRUE
        }
        params[["guideTreeOut"]] <- tempList$param
    }

    ##delete param in copy
    paramsCopy[["guideTreeOut"]] <- NULL

    #########
    # hmmIn #
    #########
    ##HMM input files
    ##default: hmmIn=NULL

    if (!is.null(params[["hmmIn"]])) {
        checkInFile("hmmIn", params)
    }

    ##delete param in copy
    paramsCopy[["hmmIn"]] <- NULL

    #########
    # inFmt #
    #########
    ##Forced sequence input file format
    ##default:infmt=auto

    ##possible Values
    posVal <- c("auto", "fa", "fasta", "clu", "clustal", "msf", "phy",
                "phylip", "selex", "st", "stockholm", "vie", "vienna")
    ##params[["inFmt"]] <- checkSingleValParams("inFmt", params,
                                                     ##"auto", posVal)
    params[["inFmt"]] <- checkSingleValParamsNew("inFmt", params, posVal)

    ##delete param in copy
    paramsCopy[["inFmt"]] <- NULL

    ##############
    # isProfile #
    ##############
    ##disable check if profile, force profile
    ##default: isProfile=FALSE

    params[["isProfile"]] <- checkLogicalParams("isProfile", params, FALSE)

    ##delete param in copy
    paramsCopy[["isProfile"]] <- NULL

    #######
    # log #
    #######
    ##Log all non-essential output to this file
    ##DEACTIVATED: All log-messages are print to R console

    ##if(!is.null(params[["log"]])) {
    ##    tempList <- checkOutFile("log", params)
    ##    if (tempList$existingFile) {
    ##        interactiveCheck("log", params)
    ##        params[["force"]] <- TRUE
    ##    }
    ##    params[["log"]] <- tempList$param
    ##}

    ##delete param in copy
    ##paramsCopy[["log"]] <- NULL

    ###############
    # longVersion #
    ###############
    ##print long version information and exit
    ##default: longVersion=FALSE

    params[["longVersion"]] <- checkLogicalParams("longVersion", params, FALSE)

    ##delete param in copy
    paramsCopy[["longVersion"]] <- NULL

    ##########
    # macRam #
    ##########
    ##maximum number guidetree iterations
    ##default: macRam = 2048

    ##params[["macRam"]] <- checkIntegerParams("macRam", params, 2048)
    params[["macRam"]] <- checkIntegerParamsNew("macRam", params)

    ##delete param in copy
    paramsCopy[["macRam"]] <- NULL

    ##########################
    # maxGuidetreeIterations #
    ##########################
    ##maximum number guidetree iterations

    ##params[["maxGuidetreeIterations"]] <- checkIntegerParams(
    ##    "maxGuidetreeIterations", params, .Machine$integer.max)
    params[["maxGuidetreeIterations"]] <- checkIntegerParamsNew(
                                    "maxGuidetreeIterations", params)

    ##delete param in copy
    paramsCopy[["maxGuidetreeIterations"]] <- NULL

    ####################
    # maxHmmIterations #
    ####################
    ##maximum number of HMM iterations

    ##params[["maxHmmIterations"]] <- checkIntegerParams(
    ##    "maxHmmIterations", params, 0) ## 0 = no Hmm is performed

    params[["maxHmmIterations"]] <- checkIntegerParamsNew(
                                    "maxHmmIterations", params)

    ##delete param in copy
    paramsCopy[["maxHmmIterations"]] <- NULL

    #############
    # maxNumSeq #
    #############
    ##maximum allowed number of sequences

    ##params[["maxNumSeq"]] <- checkIntegerParams(
    ##    "maxNumSeq", params, .Machine$integer.max)

    params[["maxNumSeq"]] <- checkIntegerParamsNew(
        "maxNumSeq", params)

    ##delete param in copy
    paramsCopy[["maxNumSeq"]] <- NULL

    #############
    # maxSeqLen #
    #############
    ##maximum allowed sequence length

    ##params[["maxSeqLen"]] <- checkIntegerParams(
    ##    "maxSeqLen", params, .Machine$integer.max)
    params[["maxSeqLen"]] <- checkIntegerParamsNew(
                                "maxSeqLen", params)

    ##delete param in copy
    paramsCopy[["maxSeqLen"]] <- NULL

    ###########
    # outfile #
    ###########
    ##Multiple Sequence Alignment Output File
    ##default: outfile=NULL

    if(!is.null(params[["outfile"]])) {
        tempList <- checkOutFile("outfile", params)
        if (tempList$existingFile) {
            interactiveCheck("outfile", params)
            params[["force"]] <- TRUE
        }
        params[["outfile"]] <- tempList$param
    }

    ##delete param in copy
    paramsCopy[["outfile"]] <- NULL

    #########
    # outFmt #
    ##########
    ##MSA output file format
    ##default in original:outfmt=fasta
    ##default (and only possibility) in this version: outfmt=clustal

    ##possible Values
    posVal <- c("auto", "fa", "fasta", "clu", "clustal", "msf", "phy",
                "phylip", "selex", "st", "stockholm", "vie", "vienna")
    ##params[["outFmt"]] <- checkSingleValParams("outFmt", params,
    ##                                           "fasta", posVal)
    params[["outFmt"]] <- checkSingleValParamsNew("outFmt", params, posVal)

    ##FIXME TODO (for later version)
    ##until now, only the clustal-format is implemented
    if (!is.null(params[["outFmt"]])) {
        if (!identical (params[["outFmt"]], "clustal") &&
            !identical (params[["outFmt"]], "clu")) {
            stop("Until now, the parameter outFmt is only implemented ",
                 "for value \"clustal\" \n",
                 "the other formats will be ",
                 "realized in a later version.")
        }
    }

    ##delete param in copy
    paramsCopy[["outFmt"]] <- NULL

    #############
    # percentId #
    #############
    ##convert distances into percent identities
    ##default: percentId=FALSE

    params[["percentId"]] <- checkLogicalParams("percentId", params, FALSE)

    ##delete param in copy
    paramsCopy[["percentId"]] <- NULL

    ############
    # profile1 #
    ############
    ##pre-aligned multiple sequence file (aligned columns will be kept fix)
    if (!is.null(params[["profile1"]])) {
        checkInFile("profile1", params)
    }

    ##delete param in copy
    paramsCopy[["profile1"]] <- NULL

    ############
    # profile2 #
    ############
    ##pre-aligned multiple sequence file (aligned columns will be kept fix)

    if (!is.null(params[["profile2"]])) {
        if (is.null(params[["profile1"]])){
            stop("The parameter profile1 is NULL, \n",
                 "so the parameter profile2 can't have a value! \n",
                 "Please insert file for parameter profile1 \n",
                 "or change the parameters profile1 and profile2!")
        }
        checkInFile("profile2", params)
    }

    ##delete param in copy
    paramsCopy[["profile2"]] <- NULL

    #################
    # residueNumber #
    #################
    ##in Clustal format print residue numbers
    ##default: residuenumber=FALSE

    params[["residueNumber"]] <- checkLogicalParams("residueNumber",
                                                    params, FALSE)

    ##delete param in copy
    paramsCopy[["residueNumber"]] <- NULL

    ###########
    # threads #
    ###########
    ##threads=<N>
    ##number of processors to use
    ##default: threads=1

    ##params[["threads"]] <- checkIntegerParams("threads", params, 1)
    params[["threads"]] <- checkIntegerParamsNew("threads", params)
    params[["threads"]] <- checkPositiveParams("threads", params)

    ##delete param in copy
    paramsCopy[["threads"]] <- NULL

    #############
    # useKimura #
    #############
    ##use Kimura distance correction for aligned sequences
    ##default: use-kimura=FALSE

    params[["useKimura"]] <- checkLogicalParams("useKimura", params, FALSE)

    ##useKimura==TRUE AND percentId==TRUE is not allowed!!!
    if(params[["useKimura"]] & params[["percentId"]]){
        stop("Percentage Identity cannot be calcuted \n",
             "if Kimura Distances are used!",
             "You have to set either the parameter percentID or \n",
             "the parameter useKimura to FALSE!")
    }

    ##delete param in copy
    paramsCopy[["useKimura"]] <- NULL

    ###########
    # version #
    ###########
    ##Print version information and exit
    ##default: version=FALSE

    params[["version"]] <- checkLogicalParams("version", params, FALSE)

    ##delete param in copy
    paramsCopy[["version"]] <- NULL

    ########
    # wrap #
    ########
    ##wrap=<N>
    ##number of residues before line-wrap in output
    ##default: wrap=60

    ##params[["wrap"]] <- checkIntegerParams("wrap", params, 60)
    params[["wrap"]] <- checkIntegerParamsNew("wrap", params)
    params[["wrap"]] <- checkPositiveParams("wrap", params)

    ##delete param in copy
    paramsCopy[["wrap"]] <- NULL


    #################
    #################
    ## FINAL CHECK ##
    #################
    #################
    ##check, whether there are additional parameters
    ##see line 23

    if (length(paramsCopy) != 0){
        stop("The following parameters are not known \n",
             "(or have been specified",
             "more often than once):\n    ",
             paste(names(paramsCopy), collapse=", ", sep=""))
    }

    inputSeqNames <- names(inputSeqs)

    names(inputSeqs) <- paste0("Seq", 1:length(inputSeqs))

    result <- .Call("RClustalOmega", inputSeqs, cluster, 6,
                    1, maxiters, substitutionMatrix, type,
                    verbose, params, PACKAGE="msa");

    out <- convertAlnRows(result$msa, type)

    if (length(inputSeqNames) > 0)
    {
        perm <- match(names(out@unmasked), names(inputSeqs))
        names(out@unmasked) <- inputSeqNames[perm]
    }
    else
        names(out@unmasked) <- NULL

    standardParams <- list(gapOpening=gapOpening,
                           gapExtension=gapExtension,
                           maxiters=maxiters,
                           verbose=verbose)

    out@params <- c(standardParams, params)
    out@call <- deparse(sys.call())
    out
}
