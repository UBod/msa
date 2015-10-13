msaMuscle <- function(inputSeqs,
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
        printHelp("Muscle")
        return(invisible(NULL))
    }

    if (!checkFunctionAvailable("Muscle"))
        stop("Muscle is not available via msa!")

    params <- list(...)
    ##create a copy of the parameter list which is only used to
    ##avoid params which are not checked:
    ##after every check of a parameter, the parameter is deleted in the copy(!)
    ##if all parameter were checked, the copy finally should be empty...
    ##if not, there are additional parameters
    paramsCopy <- lapply(params, function(x) TRUE)

    ##set names according to Biostrings

    ##set default values according to muscle_userguide3.8
    ##and check if the parameters are consistnt

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
    type <- checkType(type, inputSeqs, "msaMuscle")

    ##check the profile score
    ##necessary, because default-values for gapOpening and gapExtension
    ##depending on used profile score;
    ##center, smoothScoreCeiling, minBestColScore and minSmoothScore too
    temporaryHelp <- checkProfileScore(type, params)
    params[["le"]] <- temporaryHelp[["le"]]
    params[["sp"]] <- temporaryHelp[["sp"]]
    params[["sv"]] <- temporaryHelp[["sv"]]
    params[["spn"]] <- temporaryHelp[["spn"]]

    #############
    # inputSeqs #
    #############
    ##transform the input Sequences to a string vector
    inputSeqs <- transformInputSeq(inputSeqs)

    #############
    # order     #
    #############
    order <- match.arg(order)

    if (order == "input")
    {
        if (params[["inputSeqIsFileFlag"]])
            stop("msaMuscle does not support order=\"input\" for reading\n",
                 "sequences directly from a FASTA file.")
        else if (is.null(names(inputSeqs)) ||
                 length(unique(names(inputSeqs))) != length(inputSeqs))
        {
            warning("order=\"input\" requires input sequences to be named\n",
                    "uniquely! Assigning default names 'Seq1'..'Seqn'\n",
                    "to sequences.")
            names(inputSeqs) <- paste0("Seq", 1:length(inputSeqs))
        }
    }

    ###########
    # cluster #
    ###########
	##Perform fast clustering of input sequences.
	##Use the -tree1 option to save the tree.

    ##default-value
    if (identical(cluster, "default") || is.null(cluster)) {
        cluster <- "upgma"
    } else {
        possibleValues <- c("upgma", "upgmamax", "upgmamin",
                            "upgmb", "neighborjoining")

        ##check if cluster contains only one single value
        if(length(cluster) != 1) {
            stop("The parameter cluster contains more than one value!")
        }

        ##check if cluster is of type character
        if(!is.character(cluster)) {
            stop("The parameter cluster should be a string!")
        }

        ##make cluster more robust
        cluster <- tolower(cluster)

        ##check, if input of parameter is valid
        if (!(cluster %in% possibleValues)){
            ##create a string with all possible Values named text
            text <- ""
            text <- paste(possibleValues, collapse=", ")
            stop("The parameter cluster only can have the values: \n", text)
        }
    }


    ######################
    # substitutionMatrix #
    ######################
    ## Substitution matrix in NCBI or WU-BLAST format.
    ## If you specify your own matrix, you should also specify:
    ## -gapOpening <g>, -gapExtension <e> -center 0.0
    ## Note that <g> and <e> MUST be negative.
    ##default:
    ##protein & le => vtml_la
    ##protein & sp => PAM200
    ##protein & sv => vtml_sp
    ##dna|rna => nuc_sp

    ##default-value
    if (is.null(substitutionMatrix) ||
        identical(substitutionMatrix, "default")) {
        substitutionMatrix <- NULL;
    }

    ##check if substitutionMatrix is a matrix
    if ((!is.null(substitutionMatrix) && !is.matrix(substitutionMatrix)) ||
        identical(mode(substitutionMatrix), "list"))
        stop("The parameter substitutionMatrix should be a matrix!")

    if (!is.null(substitutionMatrix))
    {
        headerNames <- c("A", "C", "D", "E", "F",
                         "G", "H", "I", "K", "L",
                         "M", "N", "P", "Q", "R",
                         "S", "T", "V", "W", "Y")

        if (type == "protein")
            reqNames <- headerNames
        else if (type == "dna")
            reqNames <- c("A", "C", "G", "T")
        else
            reqNames <- c("A", "C", "G", "U")

        rowPerm <- match(reqNames, rownames(substitutionMatrix))
        if (any(is.na(rowPerm)))
            stop("substitutionMatrix does not contain all necessary rows")

        colPerm <- match(reqNames, colnames(substitutionMatrix))
        if (any(is.na(colPerm)))
            stop("substitutionMatrix does not contain all necessary columns")

        substitutionMatrix <- substitutionMatrix[rowPerm, colPerm]

        if (type == "rna")
            reqNames <- c("A", "C", "G", "T")

        auxMat <- matrix(0, length(headerNames), length(headerNames))
        rownames(auxMat) <- headerNames
        colnames(auxMat) <- headerNames
        auxMat[reqNames, reqNames] <- substitutionMatrix
        substitutionMatrix <- auxMat

        if (!isSymmetric(substitutionMatrix))
            stop("substitutionMatrix should be a symmetric matrix!")

        if (any(is.na(substitutionMatrix)) || any(is.na(substitutionMatrix)) ||
            any(is.infinite(substitutionMatrix)))
            stop("substitutionMatrix contains invalid values!")

        params[["le"]] <- FALSE
        params[["sv"]] <- FALSE

        if (type == "protein")
        {
            params[["sp"]] <- TRUE
            params[["spn"]] <- FALSE
        }
        else
        {
            params[["sp"]] <- FALSE
            params[["spn"]] <- TRUE
        }

        paramsCopy[["le"]] <- NULL
        paramsCopy[["sv"]] <- NULL
        paramsCopy[["sp"]] <- NULL
        paramsCopy[["spn"]] <- NULL
     }

    ##############
    # gapOpening #
    ##############
    ##The gap open score. Must be negative
    ##defaultValues:
    ##le: -2.9
    ##sp: -1439
    ##sv: -300
    ##spn_dna:400
    ##spn_rna:420

    if (params$le) {
        gapOpening <- checkGapOpening2(gapOpening, substitutionMatrix, 2.9)
    }
    else if (params$sp) {
        gapOpening <- checkGapOpening2(gapOpening, substitutionMatrix, 1439)
    }
    else if (params$sv) {
        gapOpening <- checkGapOpening2(gapOpening, substitutionMatrix, 300)
    }
    else if (params$spn) {
        if (identical(type,"dna")) {
        gapOpening <- checkGapOpening2(gapOpening, substitutionMatrix, 400)
        }
        if (identical(type,"rna")) {
            gapOpening <- checkGapOpening2(gapOpening, substitutionMatrix, 420)
        }
        if (identical(type,"protein")) {
           stop("If you use sequences of type \"protein\", \n",
                "you can't use the parameter \"spn\"!")
        }
    }

    ##FIXME TODO: check default-Value for type=protein
    ################
    # gapExtension #
    ################
    ##The gap extend score. Must be negative
    ##GapExtension only used if type=dna/rna
    ##default-Value:
    ##type= "dna" => gapExtension=0
    ##type= "rna" => gapExtension=0
    ##type= "protein" 0> gapExtension=???

    gapExtension <- checkGapExtension(gapExtension,
                                      type, substitutionMatrix, 0, 0)


    ############
    # maxiters #
    ############
    ##Maximum number of iterations

    maxiters <- checkMaxiters(maxiters, 16, "msaMuscle")


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

    ####################
    ####################
    ##  VALUE-OPTIONS ##
    ####################
    ####################

    #################
    # anchorspacing #
    #################
    ##Minimum spacing between anchor columns

    params[["anchorspacing"]] <- checkIntegerParamsNew("anchorspacing", params)

    ##delete param in copy
    paramsCopy[["anchorspacing"]] <- NULL


    ##########
    # center #
    ##########
    ##Center parameter. Should be 0 or negative
    ##default-Values:
    ##le=TRUE: -0.52
    ##spn=TRUE, type="rna": -300
    ##sp=TRUE or SV=TRUE or spn=TRUE and type="dna": 0

        params[["center"]] <- checkNumericParamsNew("center", params)
        params[["center"]] <- checkNegativeParams("center", params)


    ##delete param in copy
    paramsCopy[["center"]] <- NULL

    ############
    # cluster1 #
    ############
    ##Clustering method. cluster1 is used in iteration 1 and 2,
    ##cluster2 in later iterations.

    posVal <- c("upgma", "upgmamax", "upgmamin", "upgmb", "neighborjoining")
    params[["cluster1"]] <- checkSingleValParamsNew("cluster1", params, posVal)

    ##delete param in copy
    paramsCopy[["cluster1"]] <- NULL

    ############
    # cluster2 #
    ############
    ##Clustering method. cluster1 is used in iteration 1 and 2,
    ##cluster2 in later iterations.

    posVal <- c("upgma", "upgmb", "upgmamax", "upgmamin", "neighborjoining")
    params[["cluster2"]] <- checkSingleValParamsNew("cluster2", params, posVal)

    ##delete param in copy
    paramsCopy[["cluster2"]] <- NULL

    #############
    # diagbreak #
    #############
    ##Maximum distance between two diagonals that allows
    ##them to merge into one diagonal.

    params[["diagbreak"]] <- checkIntegerParamsNew("diagbreak", params)
    params[["diagbreak"]] <- checkPositiveParams("diagbreak", params)

    ##delete param in copy
    paramsCopy[["diagbreak"]] <- NULL

    ##############
    # diaglength #
    ##############
    ##Minimum length of diagonal

    params[["diaglength"]] <- checkIntegerParamsNew("diaglength", params)
    params[["diaglength"]] <- checkPositiveParams("diaglength", params)

    ##delete param in copy
    paramsCopy[["diaglength"]] <- NULL

    ##############
    # diagmargin #
    ##############
    ## Discard this many positions at ends of diagonal

    params[["diagmargin"]] <- checkIntegerParamsNew("diagmargin", params)
    params[["diagmargin"]] <- checkPositiveParams("diagmargin", params)

    ##delete param in copy
    paramsCopy[["diagmargin"]] <- NULL

    #############
    # distance1 #
    #############
    ##Distance measure for iteration 1.

    ##possibleValues
    if (type == "protein") {
        posVal <- c("kmer6_6", "kmer20_3", "kmer20_4", "kbit20_3")
    } else {
        posVal <- c("kmer6_6", "kmer20_3", "kmer20_4", "kbit20_3", "kmer4_6")
    }
    ##check, if input of distance1 is valid
    if (!is.null(params[["distance1"]])) {
        params[["distance1"]] <- checkIsValue("distance1", params, posVal)
    }

    ##delete param in copy
    paramsCopy[["distance1"]] <- NULL

    #############
    # distance2 #
    #############
    ##Distance measure for iterations 2, 3 ...

    params[["distance2"]] <- checkSingleValParamsNew(
                      "distance2", params, c("pctidkimura", "pctidlog"))

    ##delete param in copy
    paramsCopy[["distance2"]] <- NULL


    #########
    # hydro #
    #########
    ##Window size for determining whether a region is hydrophobic

    params[["hydro"]] <- checkIntegerParamsNew("hydro", params)

    ##delete param in copy
    paramsCopy[["hydro"]] <- NULL

    ###############
    # hydrofactor #
    ###############
    ##Multiplier for gap open/close penalties in hydrophobic regions

    params[["hydrofactor"]] <- checkNumericParamsNew("hydrofactor", params)

    ##delete param in copy
    paramsCopy[["hydrofactor"]] <- NULL

    #######
    # in1 #
    #######
    if (!is.null(params[["in1"]])) {
        checkInFile("in1", params)
    }

    ##delete param in copy
    paramsCopy[["in1"]] <- NULL

    #######
    # in2 #
    #######
    if (!is.null(params[["in2"]])) {
        checkInFile("in2", params)
    }

    ##delete param in copy
    paramsCopy[["in2"]] <- NULL

    ############
    # maxhours #
    ############
    ##Maximum time to run in hours. The actual time may exceed
    ##the requested limit by a few minutes. Decimals are allowed,
    ##so 1.5 means one hour and 30 minutes.

    params[["maxhours"]] <- checkNumericParamsNew("maxhours", params)
    params[["maxhours"]] <- checkPositiveParams("maxhours", params)

    ##delete param in copy
    paramsCopy[["maxhours"]] <- NULL

    ############
    # maxtrees #
    ############
    ##Maximum number of new trees to build in iteration 2

    params[["maxtrees"]] <- checkIntegerParamsNew("maxtrees", params)
    params[["maxtrees"]] <- checkPositiveParams("maxtrees", params)

    ##delete param in copy
    paramsCopy[["maxtrees"]] <- NULL


    ###################
    # minbestcolscore #
    ###################
    ##Minimum score a column must have to be an anchor
    ##default-Values
    ##le=TRUE=2.0
    ##sp=TRUE: 300.0
    ##sv=TRUE: 130.0
    ##spn=TRUE: 90.0

    params[["minbestcolscore"]] <- checkNumericParamsNew("minbestcolscore",
                                                         params)

    ##delete param in copy
    paramsCopy[["minbestcolscore"]] <- NULL

    ##################
    # minsmoothscore #
    ##################
    ##Minimum smoothed score a column must have to be an anchor
    ##default-Values
    ##le=TRUE=1.0
    ##sp=TRUE: 125.0
    ##sv=TRUE: 40.0
    ##spn=TRUE: 90.0

    params[["minsmoothscore"]] <- checkNumericParamsNew("minsmoothscore",
                                                        params)

    ##delete param in copy
    paramsCopy[["minsmoothscore"]] <- NULL

    ############
    # objscore #
    ############
    ##Objective score used by tree dependent refinement.
    ##sp=sum-of-pairs score.
    ##spf=sum-of-pairs score (dimer approximation)
    ##spm=sp for < 100 seqs, otherwise spf
    ##dp=dynamic programming score
    ##ps=average profile-sequence score
    ##xp=cross profile score

    posVal <- c("dp", "ps", "sp", "spf", "spm", "xp")
    params[["objscore"]] <- checkSingleValParamsNew("objscore", params, posVal)


    #delete param in copy
    paramsCopy[["objscore"]] <- NULL


    ################
    # refinewindow #
    ################
    ##Length of window for -refinew

    params[["refinewindow"]] <- checkIntegerParamsNew("refinewindow", params)
    params[["refinewindow"]] <- checkPositiveParams("refinewindow", params)

    ##delete param in copy
    paramsCopy[["refinewindow"]] <- NULL


    #########
    # root1 #
    #########
    ##Method used to root tree; root1 is used in iteration 1 and 2,
    ##root2 in later iterations

    posVal <- c("pseudo", "midlongestspan", "minavgleafdist")
    params[["root1"]] <- checkSingleValParamsNew("root1", params, posVal)

    ##delete param in copy
    paramsCopy[["root1"]] <- NULL

    #########
    # root2 #
    #########
    ##Method used to root tree; root1 is used in iteration 1 and 2,
    ##root2 in later iterations

    posVal <-  c("pseudo", "midlongestspan", "minavgleafdist")
    params[["root2"]] <- checkSingleValParamsNew("root2", params, posVal)

    ##delete param in copy
    paramsCopy[["root2"]] <- NULL

    ###################
    # smoothscoreceil #
    ###################
    ##Maximum value of column score for smoothing purposes
    ##default-Values
    ##le=TRUE=3.0
    ##sp=TRUE: 200.0
    ##sv=TRUE: 90.0
    ##spn=TRUE: 999.0

    checkNumericParamsNew("smoothscoreceil", params)
    ##delete param in copy
    paramsCopy[["smoothscoreceil"]] <- NULL


    ################
    # smoothwindow #
    ################
    ##Window used for anchor column smoothing

    params[["smoothwindow"]] <- checkIntegerParamsNew("smoothwindow", params)
    params[["smoothwindow"]] <- checkPositiveParams("smoothwindow", params)

    if (!is.null(params[["smoothwindow"]]) &&
        params[["smoothwindow"]] %% 2 == 0) {
            stop("The parameter smoothwindow must be odd!")
    }

    ##delete param in copy
    paramsCopy[["smoothwindow"]] <- NULL

    #########
    # SUEFF #
    #########
    ##Constant used in UPGMB clustering. Determines the relative
    ##fraction of average linkage (SUEFF) vs. nearest-neighbor linkage
    ##(1 - SUEFF).

    params[["SUEFF"]] <- checkNumericParamsNew("SUEFF", params)

    ##check, if input is between 0 and 1
    params[["SUEFF"]] <- checkIntervalParamsNew("SUEFF", params, 0, 1)

    ##delete param in copy
    paramsCopy[["SUEFF"]] <- NULL

    #########
    # tree1 #
    #########
    ##Save tree produced in first or second iteration to given
    ##file in Newick (Phylip-compatible) format

    ##not implemented yet

    #########
    # tree2 #
    #########
    ##Save tree produced in first or second iteration to given
    ##file in Newick (Phylip-compatible) format

    ##not implemented yet

    ###########
    # usetree #
    ###########
    ##Use given tree as guide tree.
    ##Must by in Newick (Phyip-compatible) format

    ##not implemented yet

    ###########
    # weight1 #
    ###########
    ##Sequence weighting scheme. weight1 is used in iterations 1 and 2.
    ##weight2 is used for tree-dependent refinement.
    ##none=all sequences have equal weight.
    ##henikoff=Henikoff & Henikoff weighting scheme.
    ##henikoffpb=Modified Henikoff scheme as used in PSI-BLAST.
    ##clustalw=CLUSTALW method.
    ##threeway=Gotoh three-way method.

    posVal <- c("none",
                "henikoff",
                "henikoffpb",
                "gsc",
                "clustalw",
                "threeway")

    params[["weight1"]] <- checkSingleValParamsNew("weight1", params, posVal)

    ##delete param in copy
    paramsCopy[["weight1"]] <- NULL


    ###########
    # weight2 #
    ###########
    ##Sequence weighting scheme. weight1 is used in iterations 1 and 2.
    ##weight2 is used for tree-dependent refinement.
    ##none=all sequences have equal weight.
    ##henikoff=Henikoff & Henikoff weighting scheme.
    ##henikoffpb=Modified Henikoff scheme as used in PSI-BLAST.
    ##clustalw=CLUSTALW method.
    ##threeway=Gotoh three-way method.

    posVal <- c("none",
            "henikoff",
            "henikoffpb",
            "gsc",
            "clustalw",
            "threeway")
    params[["weight2"]] <- checkSingleValParamsNew("weight2", params, posVal)

    ##delete param in copy
    paramsCopy[["weight2"]] <- NULL


    ###################
    ###################
    ##  FLAG-OPTIONS ##
    ###################
    ###################

    ###########
    # anchors #
    ###########
    ##Use anchor optimization in tree dependent refinement iterations

    if (!is.null(params[["anchors"]])) {
        params[["anchors"]] <- checkLogicalParams("anchors", params, TRUE)
    }

    ##delete param in copy
    paramsCopy[["anchors"]] <- NULL

    ###########
    # brenner #
    ###########
    ##Use Steven Brenner's method for computing the root alignment

    params[["brenner"]] <- checkLogicalParams("brenner", params, FALSE)

    ##delete param in copy
    paramsCopy[["brenner"]] <- NULL

    #######
    # clw #
    #######
    ##Write output in CLUSTALW format

    ##not activated

    ##params[["clw"]] <- checkLogicalParams("clw", params, FALSE)

    ##delete param in copy
    ##paramsCopy[["clw"]] <- NULL

    #############
    # clwstrict #
    #############
    ##Write output in CLUSTALW format with the "CLUSTAL W (1.81)"
    ##header rather than the MUSCLE version.
    ##This is useful when a post-processing step is
    ##picky about the file header

    ##set by default!!!!
    ##not activated

    ##params[["clwstrict"]] <- checkLogicalParams("clwstrict", params, FALSE)

    ##delete param in copy
    ##paramsCopy[["clwstrict"]] <- NULL

    ########
    # core #
    ########
    ##Do not catch exceptions

    if (!is.null(params[["core"]])) {
        params[["core"]] <- checkLogicalParams("core", params, TRUE)
    }

    ##delete param in copy
    paramsCopy[["core"]] <- NULL

    #########
    # diags #
    #########
    ##Use diagonal optimizations. Faster, especially for closely
    ##related sequences, but may be less accurate.

    params[["diags"]] <- checkLogicalParams("diags", params, FALSE)

    ##delete param in copy
    paramsCopy[["diags"]] <- NULL

    ##########
    # diags1 #
    ##########
    ##Use diagonal optimizations in first iteration

    params[["diags1"]] <- checkLogicalParams("diags1", params, FALSE)

    ##delete param in copy
    paramsCopy[["diags1"]] <- NULL

    ##########
    # diags2 #
    ##########
    ##Use diagonal optimizations in second iteration

    params[["diags2"]] <- checkLogicalParams("diags2", params, FALSE)

    ##delete param in copy
    paramsCopy[["diags2"]] <- NULL

    #########
    # dimer #
    #########
    ##Use dimer approximation for the SP score
    ##(faster, slightly less accurate)

    params[["dimer"]] <- checkLogicalParams("dimer", params, FALSE)

    ##delete param in copy
    paramsCopy[["dimer"]] <- NULL

    #########
    # group #
    #########
    ##Group similar sequences together in the output.
    ##This is the default. See also -stable.

    ##if (!is.null(params[["group"]])){
    ##    params[["group"]] <- checkLogicalParams("group", params, TRUE)
    ##}

    ##delete param in copy
    ##paramsCopy[["group"]] <- NULL

    ########
    # html #
    ########
    ##Write output in HTML format

    ##not implemented yet
    ##not activated

    ##params[["html"]] <- checkLogicalParams("html", params, FALSE)

    ##delete param in copy
    ##paramsCopy[["html"]] <- NULL

    ######
    # le #
    ######
    ##Use log-expectation profile score (VTML240).
    ##Alternatives are to use -sp or -sv.
    ##This is the default for amino acid sequences.

    ##param already checked within function checkProfileScore

    ##delete param in copy
    paramsCopy[["le"]] <- NULL


    #######
    # msf #
    #######
    ##Write output in MSF format (default is FASTA).
    ##Designed to be compatible with the GCG package.

    ##not implemented yet
    ##not activated

    ##params[["msf"]] <- checkLogicalParams("msf", params, FALSE)

    ##delete param in copy
    ##paramsCopy[["msf"]] <- NULL

    #############
    # noanchors #
    #############
    ##Disable anchor optimization. Default is -anchors

    params[["noanchors"]] <- checkLogicalParams("noanchors", params, FALSE)

    ##check anchors <-> noanchors
    if (!is.null(params[["anchors"]])){
        ##both positive
        if (params[["anchors"]] && params[["noanchors"]]){
            stop("The parameters anchors and noanchors \n",
                 "can't be positive at the same time!")
        }
        ##both negative
        if (!params[["anchors"]] && !params[["noanchors"]]){
            stop("The parameters anchors and noanchors \n",
                 "can't be negative at the same time!")
        }
    }
    ##delete param in copy
    paramsCopy[["noanchors"]] <- NULL

    ##########
    # nocore #
    ##########
    ##Catch exceptions and give an error message if possible

    params[["nocore"]] <- checkLogicalParams("nocore", params, FALSE)
    ##check core <-> nocore
    if (!is.null(params[["core"]])){
        ##both positive
        if (params[["core"]] && params[["nocore"]]){
            stop("The parameters core and nocore \n",
                 "can't be positive at the same time!")
        }
        ##both negative
        if (!params[["core"]] && !params[["nocore"]]){
            stop("The parameters core and nocore \n",
                 "can't be negative at the same time!")
        }
    }

    ##delete param in copy
    paramsCopy[["nocore"]] <- NULL

    ########
    # phyi #
    ########
    ##Write output in Phylip interleaved format

    ##not implemented yet
    ##not activated

    ##params[["phyi"]] <- checkLogicalParams("phyi", params, FALSE)

    ##delete param in copy
    ##paramsCopy[["phyi"]] <- NULL

    ########
    # phys #
    ########
    ##Write output in Phylip sequential format

    ##not implemented yet
    ##not activated

    ##params[["phys"]] <- checkLogicalParams("phys", params, FALSE)

    ##delete param in copy
    ##paramsCopy[["phys"]] <- NULL

    ###########
    # profile #
    ###########
    ##compute profile-profile alignment. Input alignments must be given
    ##using -in1 and -in2 options.

    params[["profile"]] <- checkLogicalParams("profile", params, FALSE)

    if (params[["profile"]]) {
        if (is.null(params[["in1"]]) || is.null(params[["in2"]])) {
            stop("The parameter profile needs the following parameters: \n",
                  "in1 and in2!")
        }
    }

    ##delete param in copy
    paramsCopy[["profile"]] <- NULL

    ##########
    # refine #
    ##########
    ##Input file is already aligned, skip first two iterations
    ##and begin tree dependent refinement.

    params[["refine"]] <- checkLogicalParams("refine", params, FALSE)

    ##delete param in copy
    paramsCopy[["refine"]] <- NULL

    ###########
    # refinew #
    ###########
    ##Refine an alignment by dividing it into non-overlapping windows
    ##and re-aligning each window.
    ##Typically used for whole-genome nucleotide alignments.

    params[["refinew"]] <- checkLogicalParams("refinew", params, FALSE)

    ##delete param in copy
    paramsCopy[["refinew"]] <- NULL

    ######
    # sp #
    ######
    ##Use sum-of-pairs protein profile score (PAM200). Default is -le

    ##param already checked within function checkProfileScore

    ##delete param in copy
    paramsCopy[["sp"]] <- NULL

    #######
    # spn #
    #######
    ##Use sum-of-pairs nucleotide profile score. This is the only option
    ##for nucleotides, and is therefore the default. The substitution scores
    ##and gap penalty scores are "borrowed" from BLASTZ.

    ##param already checked within function checkProfileScore

    ##delete param in copy
    paramsCopy[["spn"]] <- NULL

    ###########
    # spscore #
    ###########
    ##Compute alignment score of profile-profile alignment.
    ##Input alignments must be given using -in1 and -in2 options.
    ##These must be pre-aligned with gapped columns as needed,
    ##i.e. must be of the same length (have same number of columns).

    params[["spscore"]] <- checkLogicalParams("spscore", params, FALSE)

    ##delete param in copy
    paramsCopy[["spscore"]] <- NULL

    ##########
    # stable #
    ##########
    ##Preserve input order of sequences in output file.
    ##Default is to group sequences by similarity (-group).
    ##WARNING THIS OPTION WAS BUGGY AND IS NOT SUPPORTED IN v3.8

    ##The termgapshalflonger-parameter is not fully supported in version 3.8
    ##not activated
    ##params[["stable"]] <- checkLogicalParams("stable", params, FALSE)

    ##delete param in copy
    ##paramsCopy[["stable"]] <- NULL

    ######
    # sv #
    ######
    ##Use sum-of-pairs profile score (VTML240). Default is -le

    ##param already checked within function checkProfileScore

    ##delete param in copy
    paramsCopy[["sv"]] <- NULL

    #############
    # termgaps4 #
    #############
    ##Use 4-way test for treatment of terminal gaps.
    ##(Cannot be disabled in this version)

    ##The termgaps4-parameter is not fully supported in version 3.8
    ##not activated

    ##params[["termgaps4"]] <- checkLogicalParams("termgaps4", params, TRUE)

    ##delete param in copy
    ##paramsCopy[["termgaps4"]] <- NULL

    ################
    # termgapsfull #
    ################
    ##Terminal gaps penalized with full penalty

    ##The termgapsfull-parameter is not fully supported in version 3.8
    ##not activated

    ##params[["termgapsfull"]] <-
    ##checkLogicalParams("termgapsfull", params, FALSE)

    ##if (identical(params[["termgapsfull"]],TRUE)) {
    ##    warning("The parameter termgapsfull is not fully
    ##    supported in version 3.8")
    ##}

    ##delete param in copy
    ##paramsCopy[["termgapsfull"]] <- NULL

    ################
    # termgapshalf #
    ################
    ##Terminal gaps penalized with half penalty

    ##The termgapshalf-parameter is not fully supported in version 3.8
    ##not activated

    ##params[["termgapshalf"]] <-
    ## checkLogicalParams("termgapshalf", params, TRUE)

    ##delete param in copy
    ##paramsCopy[["termgapshalf"]] <- NULL

    ######################
    # termgapshalflonger #
    ######################
    ##Terminal gaps penalized with half penalty if gap relative
    ##to longer sequence, otherwise with full penalty.

    ##The termgapshalflonger-parameter is not fully supported in version 3.8
    ##not activated

    ##params[["termgapshalflonger"]] <-
    ##checkLogicalParams("termgapshalflonger", params, FALSE)

    ##if (identical(params[["termgapshalflonger"]],TRUE)) {
    ##    warning("The parameter termgapshalflonger is not fully supported
    ##    in version 3.8")
    ##}

    ##delete param in copy
    ##paramsCopy[["termgapshalflonger"]] <- NULL

    ###########
    # version #
    ###########
    ##Write version string to stdout and exit

    params[["version"]] <- checkLogicalParams("version", params, FALSE)

    ##delete param in copy
    paramsCopy[["version"]] <- NULL

    #################
    #################
    ## FINAL CHECK ##
    #################
    #################
    ##check, whether there are additional parameters
    ##see line 24

    if (length(paramsCopy) != 0){
        stop("The following parameters are not known  \n",
             "(or have been specified",
             "more often than once):\n    ",
             paste(names(paramsCopy), collapse=", ", sep=""))
    }

    inputSeqNames <- names(inputSeqs)

    names(inputSeqs) <- paste0("Seq", 1:length(inputSeqs))

    result <- .Call("RMuscle", inputSeqs, cluster, -abs(gapOpening),
                    -abs(gapExtension), maxiters, substitutionMatrix, type,
                    verbose, params, PACKAGE="msa")

    out <- convertAlnRows(result$msa, type)

    if (length(inputSeqNames) > 0)
    {
        if (order == "aligned")
        {
            perm <- match(names(out@unmasked), names(inputSeqs))
            names(out@unmasked) <- inputSeqNames[perm]
        }
        else
        {
            perm <- match(names(inputSeqs), names(out@unmasked))
            out@unmasked <- out@unmasked[perm]
            names(out@unmasked) <- inputSeqNames
        }
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
