msaClustalW <- function(inputSeqs,
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
        printHelp("ClustalW")
        return(invisible(NULL))
    }

    if (!checkFunctionAvailable("ClustalW"))
        stop("ClustalW is not available via msa!")

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
    params[["inputSeqIsFileFlag"]] <-  checkInputSeq(inputSeqs)

    ########
    # type #
    ########
    ##Sequence type, in clustalW also known as seqType

    ##validation of type
    type <- checkType(type, inputSeqs, "msaClustalW")

    #############
    # inputSeqs #
    #############
    inputSeqs <- transformInputSeq(inputSeqs)

    #############
    # order     #
    #############
    order <- match.arg(order)

    params[["outorder"]] <- order

    ###########
    # cluster #
    ###########

    ##set default values
    if (is.null(cluster) || identical(cluster, "default")){
        cluster <- "nj"
    }

    ##more than one value
    if (length(cluster) != 1) {
        stop("The parameter cluster can only have one value. \n" ,
               "Possible values are \"nj\" or \"upgma\"!")
    }
    cluster <- tolower(cluster)
    ##valid values
    if (!(cluster %in% c("nj", "upgma"))) {
        stop("The parameter cluster can only have ",
             "the values \"nj\" or \"upgma\"!")
    }

    ##FIXME TODO: check substitutionMatrix!!!
    ######################
    # substitutionMatrix #
    ######################
    ##For proteins: GONNET
    ##These matrices were derived using almost the same procedure as the
    ##PAM but are much more up to date and are based on a far larger
    ##data set. They appear to be more sensitive than the Dayhoff series.
    ##We use the GONNET 80, 120, 160, 250 and 350 matrices.
    ##This series is the default for Clustal W version 1.8.

    ##For DNA:
    ##For DNA, a single matrix (not a series) is used.
    ##Two hard-coded matrices are available:
    ##1) IUB. This is the default scoring matrix used by BESTFIT for the
    ##comparison of nucleic acid sequences. X's and N's are treated as matches
    ##to any IUB ambiguity symbol. All matches score 1.9;
    ##all mismatches for IUB symbols score 0.
    ##2) CLUSTALW(1.6). The previous system used by Clustal W, in which matches
    ##score 1.0 and mismatches score 0.
    ##All matches for IUB symbols also score 0.

    ##set both flags FALSE
    params[["substitutionMatrixIsDefaultFlag"]] <- FALSE
    params[["substitutionMatrixIsStringFlag"]] <- FALSE

    ##default-value
    if (is.null(substitutionMatrix) ||
            identical(substitutionMatrix, "default")) {
                params[["substitutionMatrixIsDefaultFlag"]] <- TRUE
    } else if (is.character(substitutionMatrix) &&
               !is.matrix(substitutionMatrix) &&
               grepl("\\.", substitutionMatrix, perl=TRUE)) {
        if (length(substitutionMatrix) != 1) {
            stop("You are using more than one file for substitutionMatrix. \n",
                 "It should only be a single character string!")
        }
        if (!file.exists(substitutionMatrix)){
            stop("The file for parameter substitutionMatrix does not exist!")
        }
        params[["substitutionMatrixIsStringFlag"]] <- TRUE
    } else if (is.character(substitutionMatrix) &&
            ##name of a matrix that should be used
        !is.matrix(substitutionMatrix)) {
        ##check whether value is BLOSUM, PAM, GONNET, or ID;
        if (type == "protein")
        {
            possibleValues <- c("blosum", "pam", "gonnet", "id")
            if (!(substitutionMatrix %in% possibleValues)){
                ##create a string with all possible Values named text
                text <- ""
                text <- paste(possibleValues, collapse=", ")
                stop("The parameter substitutionMatrix ",
                     "only can have the values: \n", text)
            }
            params[["substitutionMatrixIsStringFlag"]] <- TRUE
        }
        else
        {
            possibleValues <- c("iub", "clustalw")
            if (!(substitutionMatrix %in% possibleValues)){
                ##create a string with all possible Values named text
                text <- ""
                text <- paste(possibleValues, collapse=", ")
                stop("The parameter substitutionMatrix ",
                     "only can have the values: \n", text)
            }

            params[["substitutionMatrixIsStringFlag"]] <- FALSE
            params[["substitutionMatrixIsDefaultFlag"]] <- TRUE
            params[["pwdnamatrix"]] <- substitutionMatrix
            substitutionMatrix <- "default"
        }
    } else {
        ##real matrix
        reqNames <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
                      "B", "Z", "X", "*")

        if (type == "protein")
        {
            rowPerm <- match(reqNames, rownames(substitutionMatrix))
            if (any(is.na(rowPerm)))
                stop("substitutionMatrix does not contain all necessary rows")

            colPerm <- match(reqNames, colnames(substitutionMatrix))
            if (any(is.na(colPerm)))
                stop("substitutionMatrix does not contain all necessary columns")

            substitutionMatrix <- substitutionMatrix[rowPerm, colPerm]

            if (!isSymmetric(substitutionMatrix))
                stop("substitutionMatrix should be a symmetric matrix!")
        }
        else
        {
            reqNuc <- if (type == "dna") c("A", "G", "C", "T")
                      else c("A", "G", "C", "U")

            if (any(is.na(match(reqNuc, rownames(substitutionMatrix)))))
                    stop("substitutionMatrix does not contain all necessary rows")

            if (any(is.na(match(reqNuc, colnames(substitutionMatrix)))))
                stop("substitutionMatrix does not contain all necessary columns")

            rowSel <- which(rownames(substitutionMatrix) %in% reqNames)
            colSel <- which(colnames(substitutionMatrix) %in% reqNames)

            substitutionMatrix <- substitutionMatrix[rowSel, colSel]

            fakeAAmat <- matrix(0, length(reqNames), length(reqNames))
            rownames(fakeAAmat) <- reqNames
            colnames(fakeAAmat) <- reqNames
            fakeAAmat[rownames(substitutionMatrix), colnames(substitutionMatrix)] <-
                substitutionMatrix

            substitutionMatrix <- fakeAAmat

            params[["dnamatrix"]] <- NULL
        }
    }

    ##############
    # gapOpening #
    ##############
    ##The gap open score. Must be negative

    gapOpening <- checkGapOpening(gapOpening, type, substitutionMatrix,
                                  defaultDNAValue=15.0, defaultAAValue=10.0)

    ################
    # gapExtension #
    ################
    ##The gap extend score. Must be negative

    gapExtension <- checkGapExtension(gapExtension, type, substitutionMatrix,
                                      defaultDNAValue=6.66,
                                      defaultAAValue=0.2)

    ############
    # maxiters #
    ############
    ##Maximum number of iterations
    ##Muscle: Maxiters == ClustalW: Numiters

    maxiters <- checkMaxiters(maxiters, 3, "msaClustalW")

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

    ############################
    ############################
    ##  ***verbs(do things*** ##
    ############################
    ############################

    ###########
    # options #
    ###########
    ##list the command line parameters
    params[["options"]] <- checkLogicalParams("options", params, FALSE)

    ##delete param in copy
    paramsCopy[["options"]] <- NULL

    #########
    # check #
    #########
    ##outline the command line parameters
    params[["check"]] <- checkLogicalParams("check", params, FALSE)

    ##delete param in copy
    paramsCopy[["check"]] <- NULL

    ############
    # fullhelp #
    ############
    ##output full help content
    params[["fullhelp"]] <- checkLogicalParams("fullhelp", params, FALSE)

    ##delete param in copy
    paramsCopy[["fullhelp"]] <- NULL

    #########
    # align #
    #########
    ##do full multiple alignment
    params[["align"]] <- checkLogicalParams("align", params, FALSE)

    ##delete param in copy
    paramsCopy[["align"]] <- NULL

    ########
    # tree #
    ########
    ##calculate NJ tree
    ## params[["tree"]] <- checkLogicalParams("tree", params, FALSE)

    ##delete param in copy
    ## paramsCopy[["tree"]] <- NULL

    #######
    # pim #
    #######
    ##output percent identity matrix (while calculating the tree)
    params[["pim"]] <- checkLogicalParams("pim", params, FALSE)

    ##delete param in copy
    paramsCopy[["pim"]] <- NULL

    #############
    # bootstrap #
    #############
    ##bootstrap a NJ tree
    ##DEACTIVATED, will be realized in later version

    ##params[["bootstrap"]] <- checkLogicalParams("bootstrap", params,  FALSE)

    ##delete param in copy
    ##paramsCopy[["bootstrap"]] <- NULL

    ###############
    # bootstrapNo #
    ###############
    ##bootstrap a NJ tree
    ##DEACTIVATED, will be realized in later version

    ##if (!is.null(params[["bootstrapNo"]] ) && params[["bootstrap"]]!= TRUE){
	##	stop("If you use bootstrapNo, please set bootstrap TRUE.")
	##}
    ##params[["bootstrapNo"]] <- checkIntegerParamsNew("bootstrapNo", params)
    ##params[["bootstrapNo"]] <- checkPositiveParams("bootstrapNo", params)

    ##delete param in copy
    ##paramsCopy[["bootstrapNo"]] <- NULL

    ###########
    # convert #
    ###########    ##output the input sequences in a different file format.
    params[["convert"]] <- checkLogicalParams("convert", params, FALSE)

    ##delete param in copy
    paramsCopy[["convert"]] <- NULL

    #############################
    #############################
    ##  ***General settings*** ##
    #############################
    #############################

    #############
    # quicktree #
    #############
    ##use FAST algorithm for the alignment guide tree
    params[["quicktree"]] <- checkLogicalParams("quicktree", params, FALSE)

    ##delete param in copy
    paramsCopy[["quicktree"]] <- NULL

    ############
    # negative #
    ############
    ##protein alignment with negative values in matrix

    params[["negative"]] <- checkLogicalParams("negative", params, FALSE)

    ##delete param in copy
    paramsCopy[["negative"]] <- NULL

    ###########
    # outfile #
    ###########
    ##sequence alignment file name

    ##FIXME TODO in later version (==ILV):
    ##activate param outfile

    ##if(!is.null(params[["outfile"]])) {
    ##    tempList <- checkOutFile("outfile", params)
    ##    if (tempList$existingFile) {
    ##        interactiveCheck("outfile", params)
    ##    }
    ##    params[["outfile"]] <- tempList$param
    ##}

    ##delete param in copy
    ##paramsCopy[["outfile"]] <- NULL

    ##FIXME TODO ILV:
    ##Activate the parameter output!
    ##Until now, the only possible value for output is clustal, which is also
    ##default. A more sophisticated implementation should be available in higher
    ##versions.

    ##########
    # output #
    ##########

    ##possible Values
    ##FIXME TODO (ILV) activate the following line and
    ##deactivate the second line
    ##posVal <- c("gcg", "gde", "pir", "phylip", "nexus", "fasta", "clustal")
    posVal <- "clustal"
    ##FIXME TODO (ILV) remove this if-clause
    if (!is.null(params[["output"]]) &&
        !identical(params[["output"]], posVal)) {
        stop("Until now, the only value for parameter \n",
             "output is \"clustal\", which is default. \n",
             "A more sophisticated implementation should \n",
             "be available in higher versions of the package.")
    }
    params[["output"]] <- checkSingleValParamsNew("output", params, posVal)

    ##delete param in copy
    paramsCopy[["output"]] <- NULL

    ########
    # case #
    ########
    ##for GDE output only

    ##possible Values
    posVal <- c("lower", "upper")
    params[["case"]] <- checkSingleValParamsNew("case", params, posVal)
    ##delete param in copy
    paramsCopy[["case"]] <- NULL

    ##########
    # seqnos #
    ##########
    ##for Clustal output only

    ##Attention: param <seqnos> can still be NULL after check

    if (is.null(params[["seqnos"]])){
        params[["seqnosFlag"]] <- TRUE
        } else {
        params[["seqnosFlag"]] <- FALSE
        if (length(params[["seqnos"]]) != 1) {
            stop("The parameter seqnos should be a single string! \n",
                 "Possible values are \"on\", or \"off\"!")
        }
        if (!is.character(params[["seqnos"]])) {
            stop("The parameter <seqnos> should be a string! \n",
                 "Possible values are \"on\", or \"off\"!")
        }
        ##possible Values
        posVal <- c("on", "off")
        params[["seqnos"]] <- checkIsValue("seqnos", params, posVal)
    }

    ##delete param in copy
    paramsCopy[["seqnos"]] <- NULL

    ###############
    # seqno_range #
    ###############
    ##NEW: for all output formats

    ##possible Values
    posVal <- c("off", "on")
    params[["seqno_range"]] <- checkSingleValParamsNew("seqno_range", params,
                                                    posVal)

    ##delete param in copy
    paramsCopy[["seqno_range"]] <- NULL

    #########
    # range #
    #########
    ##RANGE=m,n
    ##sequence range to write starting m to m+n

    #default value c(-1,-1), means range not set!
    if (!is.null(params[["range"]])){
        if (length(params[["range"]])!=2){
            stop("The parameter range needs a vector of length 2! \n",
                 "Both values should be positive integers!")
        }
        if (!is.vector(params[["range"]])) {
            stop("The parameter range should be a vector ",
                 "with 2 positive integers!")
        }
        if (any(is.na(params[["range"]])) || any(is.nan(params[["range"]]))) {
            stop("The parameter range should consist of 2 positive ",
                   "integers, \n",
                   "not with NAs or NaNs!")
        }
        if (!is.integer(params[["range"]])) {
            ##stop if usage of floats
            if (params[["range"]][[1]] - round(params[["range"]][[1]]) != 0 |
                params[["range"]][[2]] - round(params[["range"]][[2]]) != 0) {
                stop("The parameter range should consist of integers, \n",
                     "not numeric values!")
            }
            ##typecast if possible
            if (params[["range"]][[1]] <= .Machine$integer.max &
                params[["range"]][[2]] <= .Machine$integer.max) {
                params[["range"]][[1]] = as.integer(params[["range"]][[1]])
                params[["range"]][[2]] = as.integer(params[["range"]][[2]])
            } else {
                stop("The values in parameter range ",
                     " are bigger than integer!")
            }
        }
        if (params[["range"]][[1]] < 0 | params[["range"]][[2]] < 0) {
            stop("The parameter range needs positive integer values!")
        }
    }

    ##delete param in copy
    paramsCopy[["range"]] <- NULL

    #############
    # maxseqlen #
    #############
    ##maximum allowed input sequence length

    ##deactivated due to the fact, that the input sequence can also be a file.
    ##If this is the case, and maxseqlen is setted, problems occur in msa().
    ##A possible catch in R is only possible with loss of performance, so we
    ##decided to renounce this unnecessary param

    ##params[["maxseqlen"]] <- checkIntegerParamsNew("maxseqlen", params)
    ##delete param in copy
    ##paramsCopy[["maxseqlen"]] <- NULL

    #########
    # stats #
    #########
    ##log some alignment statistics to file

    if(!is.null(params[["stats"]])) {
        tempList <- checkOutFile("stats", params)
        if (tempList$existingFile) {
            interactiveCheck("stats", params)
        }
        params[["stats"]] <- tempList$param
    }

    ##delete param in copy
    paramsCopy[["stats"]] <- NULL

    ####################################
    ####################################
    ##  ***Fast Pairwise Alignment*** ##
    ####################################
    ####################################

    ##########
    # ktuple #
    ##########
    ##Fast pairwise alignment word size used to find matches between
    ##the sequences.Decrease for sensitivity; increase for speed.

    params[["ktuple"]] <- checkIntegerParamsNew("ktuple", params)
    ##constraint: if type=="protein", only ktuple <= 2 allowed
    if (!is.null(params[["ktuple"]])) {
        if (type=="protein") {
            if (params[["ktuple"]] > 2){
                stop("If you are using proteins, ktuple should be <=2!")
            }
        }
    }

    ##delete param in copy
    paramsCopy[["ktuple"]] <- NULL

    ############
    # topdiags #
    ############
    ##Fast pairwise alignment number of match regions are used to create
    ##the pairwise alignment. Decrease for speed; increase for sensitivity.

    params[["topdiags"]] <- checkIntegerParamsNew("topdiags", params)


    ##delete param in copy
    paramsCopy[["topdiags"]] <- NULL

    ##########
    # window #
    ##########
    ##Fast pairwise alignment window size for joining word matches.
    ##Decrease for speed; increase for sensitivity.


    params[["window"]] <- checkIntegerParamsNew("window", params)


    ##delete param in copy
    paramsCopy[["window"]] <- NULL

    ###########
    # pairgap #
    ###########
    ##Fast pairwise alignment gap penalty for each gap created.


    params[["pairgap"]]  <- checkIntegerParamsNew("pairgap", params)


    ##delete param in copy
    paramsCopy[["pairgap"]] <- NULL

    #########
    # score #
    #########
    ##Fast pairwise alignment score type to output.

    posVal <- c("percent", "absolute")
    params[["score"]] <- checkSingleValParamsNew("score", params, posVal)

    ##delete param in copy
    paramsCopy[["score"]] <- NULL

    ####################################
    ####################################
    ##  ***Slow Pairwise Alignment*** ##
    ####################################
    ####################################

    ############
    # pwmatrix #
    ############
    ##Slow pairwise alignment protein sequence comparison matrix series
    ##used to score alignment.

    ##if filename (seperated with ".") check and use file;
    ##else check whether value is BLOSUM, PAM, GONNET, or ID;
    ##if nothing is given, use GONNET
    if (!is.null(params[["pwmatrix"]]) &&
            grepl("\\.", params[["pwmatrix"]], perl=TRUE)) {
        checkInFile("pwmatrix", params)
    } else {
        posVal <- c("blosum", "pam", "gonnet", "id")
        params[["pwmatrix"]] <- checkSingleValParamsNew(
                                "pwmatrix", params, posVal)
    }

    ##delete param in copy
    paramsCopy[["pwmatrix"]] <- NULL

    ###############
    # pwdnamatrix #
    ###############
    ##Slow pairwise alignment nucleotide sequence comparison matrix
    ##used to score alignment.

    ##if filename (seperated with ".") check and use file;
    ##else check whether value is iub or clustalw;
    ##if nothing is given, use iub
    if (!is.null(params[["pwdnamatrix"]]) &&
            grepl("\\.", params[["pwdnamatrix"]], perl=TRUE)) {
        checkInFile("pwdnamatrix", params)
    } else {
        posVal <- c("iub", "clustalw")
        params[["pwdnamatrix"]] <- checkSingleValParamsNew("pwdnamatrix",
                                                           params, posVal)
    }

    ##delete param in copy
    paramsCopy[["pwdnamatrix"]] <- NULL

    #############
    # pwgapopen #
    #############
    ##Slow pairwise alignment score for the first residue in a gap.

    params[["pwgapopen"]] <- checkNumericParamsNew("pwgapopen", params)

    if (is.numeric(params[["pwgapopen"]]))
        params[["pwgapopen"]] <- abs(params[["pwgapopen"]])

    ##delete param in copy
    paramsCopy[["pwgapopen"]] <- NULL

    ############
    # pwgapext #
    ############
    ##Slow pairwise alignment score for each additional residue in a gap.

    params[["pwgapext"]] <- checkNumericParamsNew("pwgapext", params)

    if (is.numeric(params[["pwgapext"]]))
        params[["pwgapext"]] <- abs(params[["pwgapext"]])

    ##delete param in copy
    paramsCopy[["pwgapext"]] <- NULL

    ################################
    ################################
    ##  ***Multiple Alignments*** ##
    ################################
    ################################

    ###########
    # newtree #
    ###########
    ##file for new guide tree

    ##FIXME TODO (ILV):
    ##activate the parameter newtree

    ##if(!is.null(params[["newtree"]])) {
    ##    tempList <- checkOutFile("newtree", params)
    ##    if (tempList$existingFile) {
    ##        interactiveCheck("newtree", params)
    ##    }
    ##    params[["newtree"]] <- tempList$param
    ##}

    ##delete param in copy
    ##paramsCopy[["newtree"]] <- NULL

    ###########
    # usetree #
    ###########
    ##file for old guide tree

    if (!is.null(params[["usetree"]])) {
        checkInFile("usetree", params)
    }

    ##delete param in copy
    paramsCopy[["usetree"]] <- NULL

    #############
    # dnamatrix #
    #############
    ##DNA weight matrix=IUB, CLUSTALW or filename

    ##if filename (separated with ".") check and use file;
    ##else check whether value is iub or clustalw;
    ##if nothing is given, use iub
    if (!is.null(params[["dnamatrix"]]) &&
            grepl("\\.", params[["dnamatrix"]], perl=TRUE)) {
        checkInFile("dnamatrix", params)
    } else if (is.null(params[["pwdnamatrix"]])) {
        posVal <- c("iub", "clustalw")
        params[["pwdnamatrix"]] <- checkSingleValParamsNew("dnamatrix",
                                                           params, posVal)
    }

    ##delete param in copy
    paramsCopy[["dnamatrix"]] <- NULL


    ###########
    # endgaps #
    ###########
    ##no end gap separation penalty

    params[["endgaps"]] <- checkLogicalParams("endgaps", params, FALSE)

    ##delete param in copy
    paramsCopy[["endgaps"]] <- NULL

    ###########
    # gapdist #
    ###########
    ##gap separation pen. range

    params[["gapdist"]] <- checkIntegerParamsNew("gapdist", params)
    ##delete param in copy
    paramsCopy[["gapdist"]] <- NULL

    ##########
    # nopgap #
    ##########
    ##residue-specific gaps off

    params[["nopgap"]] <- checkLogicalParams("nopgap", params, FALSE)

    ##delete param in copy
    paramsCopy[["nopgap"]] <- NULL

    ##########
    # nohgap #
    ##########
    ##hydrophilic gaps off

    params[["nohgap"]] <- checkLogicalParams("nohgap", params, FALSE)

    ##delete param in copy
    paramsCopy[["nohgap"]] <- NULL

    ##########
    # novgap #
    ##########
    ##
    if (!is.null(params[["novgap"]])) {
        params[["novgap"]] <- checkLogicalParams("novgap", params, TRUE)
    }

    ##delete param in copy
    paramsCopy[["novgap"]] <- NULL

    ################
    # hgapresidues #
    ################
    ##list hydrophilic res.

    if (!is.null(params[["hgapresidues"]])){
        ##check if hgapresidues is a string
        if (!is.character(params[["hgapresidues"]])){
            stop("The parameter hgapresidues should be a string!")
        }
    }

    ##delete param in copy
    paramsCopy[["hgapresidues"]] <- NULL

    ##########
    # maxdiv #
    ##########
    ##% ident. for delay
    params[["maxdiv"]] <- checkIntegerParamsNew("maxdiv", params)
    ##delete param in copy
    paramsCopy[["maxdiv"]] <- NULL

    ###############
    # transweight #
    ###############
    ##transitions weighting
    params[["transweight"]] <- checkNumericParamsNew("transweight", params)
    ##delete param in copy
    paramsCopy[["transweight"]] <- NULL


    #############
    # iteration #
    #############

    ##possible Values
    posVal <- c("tree", "alignment", "none")
    params[["iteration"]] <- checkSingleValParamsNew("iteration",
                                                     params,  posVal)

    ##delete param in copy
    paramsCopy[["iteration"]] <- NULL

    #############
    # noweights #
    #############
    ##disable sequence weighting

    params[["noweights"]] <- checkLogicalParams("noweights", params, FALSE)

    ##delete param in copy
    paramsCopy[["noweights"]] <- NULL

    ###############################
    ###############################
    ##  ***Profile Alignments*** ##
    ###############################
    ###############################

    ###########
    # profile #
    ###########
    ##Merge two alignments by profile alignment

    params[["profile"]] <- checkLogicalParams("profile", params, FALSE)

    ##delete param in copy
    paramsCopy[["profile"]] <- NULL

    ############
    # profile1 #
    ############

    if (!is.null(params[["profile1"]])) {
        checkInFile("profile1", params)
    }

    ##delete param in copy
    paramsCopy[["profile1"]] <- NULL

    ############
    # profile2 #
    ############

    if (!is.null(params[["profile2"]])) {
        checkInFile("profile2", params)
    }

    ##delete param in copy
    paramsCopy[["profile2"]] <- NULL

    ############
    # newtree1 #
    ############
    ##file for new guide tree for profile1

    ##FIXME TODO ILV
    ##activate param newtree1

    ##if(!is.null(params[["newtree1"]])) {
    ##    tempList <- checkOutFile("newtree1", params)
    ##    if (tempList$existingFile) {
    ##        interactiveCheck("newtree1", params)
    ##    }
    ##   params[["newtree1"]] <- tempList$param
    ##}

    ##delete param in copy
    ##paramsCopy[["newtree1"]] <- NULL

    ############
    # newtree2 #
    ############
    ##file for new guide tree for profile2

    ##FIXME TODO ILV
    ##activate param newtree1

    ##if(!is.null(params[["newtree2"]])) {
    ##    tempList <- checkOutFile("newtree2", params)
    ##    if (tempList$existingFile) {
    ##        interactiveCheck("newtree2", params)
    ##    }
    ##    params[["newtree2"]] <- tempList$param
    ##}

    ##delete param in copy
    ##paramsCopy[["newtree2"]] <- NULL

    ############
    # usetree1 #
    ############
    ##file for old guide tree for profile1


    if (!is.null(params[["usetree1"]])) {
        checkInFile("usetree1", params)
    }

    ##delete param in copy
    paramsCopy[["usetree1"]] <- NULL

    ############
    # usetree2 #
    ############
    ##file for old guide tree for profile2


    if (!is.null(params[["usetree2"]])) {
        checkInFile("usetree2", params)
    }

    ##delete param in copy
    paramsCopy[["usetree2"]] <- NULL

    ###########################################
    ###########################################
    ##  ***Sequence to Profile Alignments*** ##
    ###########################################
    ###########################################

    #############
    # sequences #
    #############
    ##Sequentially add profile2 sequences to profile1 alignment

    params[["sequences"]] <- checkLogicalParams("sequences", params, FALSE)

    ##delete param in copy
    paramsCopy[["sequences"]] <- NULL

    ############
    # newtree #
    ############
    ##file for new guide tree

    ##already implemented in Multiple Alignments

    ###########
    # usetree #
    ###########
    ##file for old guide tree

    ##already implemented in Multiple Alignments

    ##################################
    ##################################
    ##  ***Structural Alignments*** ##
    ##################################
    ##################################

    #############
    # nosecstr1 #
    #############
    ##do not use secondary structure-gap penalty mask for profile 1

    params[["nosecstr1"]] <- checkLogicalParams("nosecstr1", params, FALSE)

    ##delete param in copy
    paramsCopy[["nosecstr1"]] <- NULL

    #############
    # nosecstr2 #
    #############
    ##do not use secondary structure-gap penalty mask for profile 2

    params[["nosecstr2"]] <- checkLogicalParams("nosecstr2", params, FALSE)

    ##delete param in copy
    paramsCopy[["nosecstr2"]] <- NULL

    #############
    # secstrout #
    #############
    ##output in alignment file

    ##possible Values
    posVal <- c("structure", "mask", "both", "none")
    params[["secstrout"]] <- checkSingleValParamsNew("secstrout",
                                                     params, posVal)

    ##delete param in copy
    paramsCopy[["secstrout"]] <- NULL

    ############
    # helixgap #
    ############
    ##gap penalty for helix core residues

    params[["helixgap"]] <- checkIntegerParamsNew("helixgap", params)
    ##delete param in copy
    paramsCopy[["helixgap"]] <- NULL

    #############
    # strandgap #
    #############
    ##gap penalty for strand core residues

    params[["strandgap"]] <- checkIntegerParamsNew("strandgap", params)
    ##delete param in copy
    paramsCopy[["strandgap"]] <- NULL

    ###########
    # loopgap #
    ###########
    ##gap penalty for loop regions

    params[["loopgap"]] <- checkIntegerParamsNew("loopgap", params)
    ##delete param in copy
    paramsCopy[["loopgap"]] <- NULL

    ###############
    # terminalgap #
    ###############
    ##gap penalty for structure termini

    params[["terminalgap"]] <- checkIntegerParamsNew("terminalgap", params)
    ##delete param in copy
    paramsCopy[["terminalgap"]] <- NULL

    ##############
    # helixendin #
    ##############
    ##number of residues inside helix to be treated as terminal

    params[["helixendin"]] <- checkIntegerParamsNew("helixendin", params)
    ##delete param in copy
    paramsCopy[["helixendin"]] <- NULL


    ###############
    # helixendout #
    ###############
    ##number of residues outside helix to be treated as terminal

    params[["helixendout"]] <- checkIntegerParamsNew("helixendout", params)
    ##delete param in copy
    paramsCopy[["helixendout"]] <- NULL

    ###############
    # strandendin #
    ###############
    ##number of residues inside strand to be treated as terminal

    params[["strandendin"]] <- checkIntegerParamsNew("strandendin", params)
    ##delete param in copy
    paramsCopy[["strandendin"]] <- NULL

    ################
    # strandendout #
    ################
    ##number of residues outside strand to be treated as terminal

    params[["strandendout"]] <- checkIntegerParamsNew("strandendout", params)
    ##delete param in copy
    paramsCopy[["strandendout"]] <- NULL

    #################
    #################
    ## ***Trees*** ##
    #################
    #################

    ##############
    # outputtree #
    ##############
    ##

    ##possible Values
    posVal <- c("nj", "phylip", "dist", "nexus")
    params[["outputtree"]] <- checkSingleValParamsNew("outputtree",
                                                      params, posVal)

    ##delete param in copy
    paramsCopy[["outputtree"]] <- NULL

    ########
    # seed #
    ########
    ##seed number of bootstraps

    params[["seed"]] <- checkIntegerParamsNew("seed", params)

    ##delete param in copy
    paramsCopy[["seed"]] <- NULL

    ##########
    # kimura #
    ##########
    ##use Kimura's correction

    params[["kimura"]] <- checkLogicalParams("kimura", params, FALSE)

    ##delete param in copy
    paramsCopy[["kimura"]] <- NULL

    ############
    # tossgaps #
    ############
    ##ignore positions with gaps

    params[["tossgaps"]] <- checkLogicalParams("tossgaps", params, FALSE)

    ##delete param in copy
    paramsCopy[["tossgaps"]] <- NULL

    ##############
    # bootlabels #
    ##############
    ##position of bootstrap values in tree display

    ##possible Values
    posVal <- c("node", "branch")
    params[["bootlabels"]] <- checkSingleValParamsNew("bootlabels",
                                                   params, posVal)

    ##delete param in copy
    paramsCopy[["bootlabels"]] <- NULL

    #################
    #################
    ## FINAL CHECK ##
    #################
    #################
    ##check, whether there are additional parameters
    ##see line 25

    if (length(paramsCopy) != 0){
        stop("The following parameters are not known \n",
             "(or have been specified",
             "more often than once):\n    ",
             paste(names(paramsCopy), collapse=", ", sep=""))
    }

    inputSeqNames <- names(inputSeqs)

    names(inputSeqs) <- paste0("Seq", 1:length(inputSeqs))

    result <- .Call("RClustalW", inputSeqs, cluster, abs(gapOpening),
                    abs(gapExtension), maxiters, substitutionMatrix,
                    type, verbose, params, PACKAGE="msa")

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
