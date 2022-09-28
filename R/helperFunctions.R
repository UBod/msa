## Function that returns the Input-Parameters as a character string
## in the form "param1Name = param1Value, param2Name = param2Value, ..."
printParams <-function(type,
                       cluster,
                       gapOpening,
                       gapExtension,
                       maxiters,
                       verbose,
                       params){
    result <-  capture.output(cat(
                    rbind(  "cluster", "=", cluster, ",",
                            "gapOpening", "=" , gapOpening, ", ",
                            "gapExtension", "=" , gapExtension,  ", ",
                            "maxiters", "=" , maxiters,  ", ",
                            "type", "=" , type,  ", ",
                            "verbose", "=" , verbose),
                    as.character(rbind(",", names(params), "=" ,params))))
    return (result)
}

##############################################################################

## Function that returns the Input-Parameters as a list
getParamsList <-function(type,
                         cluster,
                         gapOpening,
                         gapExtension,
                         maxiters,
                         verbose,
                         params){
    result <- list()
    result$type <- type
    result$cluster <- cluster
    result$gapOpening <- gapOpening
    result$gapExtension <- gapExtension
    result$maxiters <- maxiters
    result$verbose <- verbose
    result$params <- params #list of all other params
    return (result)
}

###############################################################################

##function, which tests a given input sequence
##whether it is a character string, or a XStringSet
transformInputSeq <- function(inputSeq) {
    if (!is(inputSeq, "character")) {
        if (is(inputSeq, "XStringSet")) {
            inputSeq <- as.character(inputSeq)
        } else {
            stop("The Parameter inputSeq is not valid. \n",
                 "Possible inputs are <character> or <XStringSet>!")
        }
    }

    return(toupper(inputSeq))
}

###############################################################################

##function which returns the type of a given input sequence
##returns "dna", "rna" or "protein" (depending of the input),
##if input sequence is NOT a BStringSet or a character string.
##The later two return NULL

getTypeOfInputSeq <- function(inputSeq) {
    if (!is(inputSeq, "character")) {
        if (is(inputSeq, "AAStringSet")) {
            type <- "protein"
        } else if (is(inputSeq, "DNAStringSet")) {
            type <- "dna"
        } else if (is(inputSeq, "RNAStringSet")) {
            type <- "rna"
        } else if (is(inputSeq, "BStringSet")) {
            type <- NULL
        } else {
            stop("The Parameter inputSeq is not valid.\n",
                 "Possible inputs are <character>, or <XStringSet>!")
        }
    } else {
        type <- NULL
    }
    return(type)
}

###############################################################################

##auxiliary function, that checks the type, whether it is NULL, or not one
##of the following possible values c("protein", "dna", "rna"), the later
##not for T-Coffee. In both cases, an exception is thrown.
##If no exception is thrown, it returns the type.
checkOneType <- function(type, msaName) {

    ##default-value
    if (is.null(type)) {
        stop("Your input sequence is a character string, a filename,\n",
             "or a BStringSet object. So please specify the sequence type.\n",
             "Possible values are: \"protein\", \"dna\", or \"rna\"!")
    }

    ##make it more robust
    type=tolower(type)

    ##check, if input of type is valid
    if (tolower(msaName) %in% c("msaclustalw",
                                "msaclustalomega",
                                "msamuscle")) {
        if (!type %in% c("protein", "dna", "rna")) {
            stop("The type parameter only can have the values ",
                   "\"protein\", \"dna\", or \"rna\"!")
        }
    }

    return(type)
}


###############################################################################

##auxiliary function, which tests a given type with the type of inputSeqs
##(only possible, if type is given AND the input Sequence is a
##XStringSet(without BStringSet)
checkDoubleGivenType <- function(type1, type2) {
    if (!(length(type1) == 0) & !(length(type2) == 0)) {
        if (!is.na(type1) & !is.na(type2)) {
            if (type1 != type2){
                stop("The param type is in conflict to the input Sequence!")
            }
        }
    }
}

###############################################################################


##function that converts a string of type "paramVal1,paremVal2,etc."
##(only seperated with commatas, no blanks in between)
##into a vector of type c("ParamVal1", "ParamVal2", "etc.")
convertStringToVector <- function(parameterName, params){
    helpVector <- c()
    ##default-value
    if (is.null(params[[parameterName]])) {
        helpVector <- NULL
    } else {
        ##make it robust
        if (!is.character(params[[parameterName]])) {
            stop("The parameter ", parameterName, " should be a string!")
        }
        ##do not allow empty strings
        if (params[[parameterName]] == "") {
            stop("The parameter ", parameterName,
                            " can't be an empty string!")
        }
        ##convert
        for (i in 1:length(strsplit(params[[parameterName]], split=",",
                        fixed=FALSE, perl=FALSE, useBytes=FALSE)[[1]])){
            helpVector <- c(helpVector, strsplit(params[[parameterName]],
                            split=",", fixed=FALSE, perl=FALSE,
                            useBytes=FALSE)[[1]][[i]])
        }
    }
    return(helpVector)
}

###############################################################################

##auxiliary function, that tests a param, which has exact one single character
##string, if it is of a set of possible values. If yes, the param returns.
checkIsValue <- function(parameterName, params, possibleValues){

    ##check if input is NULL
    if (is.null(params[[parameterName]])) {
        stop("The parameter ", parameterName,
                        " contains no value!")
    }

    ##check if params$parameterName contains only one single value
    if(length(params[[parameterName]]) != 1) {
        stop("The parameter ", parameterName,
                        " contains more than one value!")
    }

    ##check if is params$parameterName of type character
    if(!is.character(params[[parameterName]])) {
        stop("The parameter ", parameterName, " should be a string!")
    }

    ##check, if input of parameter is valid
    if (!(tolower(params[[parameterName]]) %in% possibleValues)){
        ##create a string with all possible Values named text
        text <- ""
        text <- paste(possibleValues, collapse=", ")
        stop("The parameter ", parameterName,
                        " only can have the values: \n", text)
    }

    return(tolower(params[[parameterName]]))
}

###############################################################################
##function for an interactive check, whether a output file should be
##overwritten or not. If not, the function stops.

interactiveCheck <- function(parameterName, params){
    message(c("File for param ", parameterName, " exists. Overwrite? (y/n)"))
    answer <- try(tolower(scan(what=character(), nmax=1,
                            quiet=TRUE)), silent=TRUE)
    if (nchar(answer) != 1 || substr(answer, 1, 1) != "y") {
        stop("The file for param ", parameterName, " exists already! \n",
             "You didn't allow overwriting. So remove/rename \n",
             "the file and try it again (or allow overwriting)!")
    }
}

###############################################################################

##auxiliary function that replaces backslashes in file names by slashes on
##Windows systems
stratifyFilenames <- function(x)
{
    if (identical(.Platform$OS.type, "windows"))
        gsub("\\", "/", x, fixed=TRUE)
    else
        x
}


###############################################################################

##auxiliary function that prints the help with all common issues.
##Furthermore, it adds a specific help (depending on algorithm)

printHelp <- function(method)
{
    file <- system.file("extdata/msaCommonHelp.txt", package="msa")
    cat(readLines(file), sep="\n")

    cat ("Press [enter] to continue")
    line <- readline()

    file <- system.file(paste0("extdata/msa", method, "Help.txt"),
                        package="msa")
    file.show(file)
}
