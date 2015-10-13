##function, that tests verbose, set the default value if necessary
##and returns verbose.
checkVerbose <- function(defaultValue, verbose){
    ##default-value
    if (is.null(verbose) || identical(verbose, "default")) {
        verbose <- defaultValue
    }

    if (!(identical(verbose,TRUE)|identical(verbose,FALSE))) {
        stop("The verbose parameter must be logical!")
    }
    return (verbose)
}

###############################################################################

##function, which tests a given input sequence
##whether it is a character string, a XStringSet,  or a
##file name.
##if file name, it returns a flag TRUE, else FALSE
checkInputSeq <- function(inputSeq) {
    result <- FALSE
    if (is(inputSeq, "character")) {
        ##check if Input is File
        if (length(inputSeq) == 1) {
            ##check whether file-ending is ".fa" or ".fasta"
            if(grepl("\\.fa", inputSeq, perl=TRUE) ||
                    grepl("\\.fasta", inputSeq, perl=TRUE)) {
                ##check whether file exists
                if (file.exists(inputSeq)) {
                    result <- TRUE
                } else {
                    stop("The file for inputSeq does not exist!")
                }
            } else {
                ##any other file
                if(grepl("\\.", inputSeq, perl=TRUE)) {
                    stop("For inputSeq, only \".fasta\", or \".fa\" -Files \n",
                                    "are allowed!")
                }
            }
        }
    } else {
        if (!is(inputSeq, "XStringSet")) {
            stop("The parameter inputSeq is not valid! \n",
                 "Possible inputs are <character>, <XStringSet>, or a file.")
        }
    }
    return(result)
}


###############################################################################

##function, that tests the type with the two auxiliary functions
##checkOneType() and checkDoubleGivenType()
checkType <- function(type, inputSeqs, msaName){
    type2 <- getTypeOfInputSeq(inputSeqs)

    if(is.null(type) || identical(type, "default")) {
        ##type <- type of inputSeqs
        type <- type2
    }

    ##validation of type
    type <- checkOneType(type, msaName)
    ##check, if type == type of inputSeqs
    checkDoubleGivenType(type, type2)
    return(type)
}

###############################################################################

##function, that tests the input of gapOpening.
##If the value is numeric, everything is ok and the function returns
##the gapOpening parameter. If the input is not numeric, an exception is thrown.
##Same for missing substitutionMatrix
checkGapOpening <- function(gapOpening, type, substitutionMatrix,
        defaultDNAValue,  defaultAAValue){
    if (is.null(gapOpening) || identical(gapOpening, "default")) {
        if (type  == "protein"){
            gapOpening <- defaultAAValue
        } else {
            gapOpening <- defaultDNAValue
        }
    }

    ##check, if input of gapOpening is valid
    if (is.numeric(gapOpening)) {
        if (is.matrix(gapOpening)) {
            stop("The parameter gapOpening should be \n",
                 "a numeric, not a matrix!")
        }
        if (length(gapOpening) != 1) {
            stop("The parameter gapOpening should be \n",
                 "a numeric, not a vector!")
        }
        if (is.nan(gapOpening)) {
            stop("The parameter gapOpening should be \n",
                 "a numeric, not a NaN!")
        }
    } else {
        stop("The parameter gapOpening should be a numeric!")
    }
    return(abs(gapOpening))
}

###############################################################################

##used in MUSCLE, analoguous to checkGapOpening, but only ONE defaut value
##function, that tests the input of gapOpening.
##If the value is numeric, everything is ok and the function returns
##the gapOpening parameter. If the input is not numeric, an exception is thrown.
checkGapOpening2 <- function(gapOpening, substitutionMatrix,
        defaultValue){
    ##set defaultValue
    if (is.null(gapOpening) || identical(gapOpening, "default")) {
        gapOpening <- defaultValue
    }


    ##check, if input of gapOpening is valid
    if (is.numeric(gapOpening)) {
        if (is.matrix(gapOpening)) {
            stop("The parameter gapOpening should be \n",
                 "a numeric, not a matrix!")
        }
        if (length(gapOpening) != 1) {
            stop("The parameter gapOpening should be \n",
                 "a numeric, not a vector!")
        }
        if (is.nan(gapOpening)) {
            stop("The parameter gapOpening should be \n",
                 "a numeric, not a NaN!")
        }
    } else {
        stop("The parameter gapOpening should be a numeric!")
    }
    return(abs(gapOpening))
}

###############################################################################

##function analoguous to checkGapOpening
checkGapExtension <- function(gapExtension, type, substitutionMatrix,
        defaultDNAValue, defaultAAValue){
    if (is.null(gapExtension) || identical(gapExtension, "default"))  {
        if (type  == "protein"){
            gapExtension <- defaultAAValue
        } else {
            gapExtension <- defaultDNAValue
        }
    }

    ##check, if input of gapExtension is valid
    if (is.numeric(gapExtension)) {
        if (is.matrix(gapExtension)) {
           stop("The parameter gapExtension should be \n",
                "a numeric, not a matrix!")
        }
        if (length(gapExtension) != 1) {
           stop("The parameter gapExtension should be \n",
                "a numeric, not a vector!")
        }
        if (is.nan(gapExtension)) {
            stop("The parameter gapExtension should be \n",
                 "a numeric, not a NaN!")
        }
    } else {
        stop("The parameter gapExtension should be a numeric!")
    }

    return(abs(gapExtension))
}

###############################################################################

##function, that tests the input of maxIters.
##set the default value, if necessary
##stops, if not using positive integers
checkMaxiters <- function(maxIters, defaultValue, algorithmName){
    ##default-value
    if(is.null(maxIters)|| identical(maxIters, "default")) {
        maxIters <- defaultValue
    }

    ##check, if input of maxiters is valid
    if (length(maxIters) != 1) {
        stop("The parameter maxiters should be a single positive integer!")
    }
    if (is.integer(maxIters)) {
        if (maxIters < 0) {
            stop("The parameter maxiters should be a positive integer!")
        }
        ##stop if using 0 in Muscle or ClustalW
        if (algorithmName %in% c("msaMuscle", "msaClustalW") &&
                maxIters == 0) {
            stop("The parameter maxiters should be a positive integer!")
        }
    } else {
        if (is.numeric(maxIters)) {
            if (is.matrix(maxIters)) {
                stop("The parameter maxiters should be a positive integer,\n",
                     "not a matrix!")
            }
            if (length(maxIters) != 1) {
                stop("The parameter maxiters should be a positive integer,\n",
                        "not a vector!")
            }
            if (is.nan(maxIters)) {
                stop("The parameter maxiters should be a negative numeric,\n",
                     "not a NaN!")
            }
            ##stop if usage of floats
            if (maxIters - round(maxIters) != 0) {
                stop("The parameter maxiters should be a positive integer!")
            }

            ##stop if using maxiters <= 0 in Muscle or ClustalW
            if (algorithmName %in% c("msaMuscle", "msaClustalW") &&
                    maxIters <= 0) {
                stop("The parameter maxiters should be a positive integer!")
            }
            ##stop if using maxiters < 0 in ClustalOmega
            if (identical(algorithmName, "msaClustalOmega") && maxIters < 0) {
                stop("The parameter maxiters should be a positive integer!")
            }
            ##typecast
            if (maxIters < .Machine$integer.max) {
                maxIters <- as.integer(maxIters)
            } else {
                stop("The parameter maxiters is bigger than an integer!")
            }
        } else {
            stop("The parameter maxiters should be a positive integer!")
        }
    }
    return(maxIters)
}
###############################################################################

##function, that tests a param whether it is logical or not and if
##the default value needs to be set. If it isn't logical,
##an exception is thrown, otherwise, the function returns the param
checkLogicalParams <- function(parameterName, params, defaultValue){
    ##default-value
    if (is.null(params[[parameterName]])) {
        params[[parameterName]] <- defaultValue
    }

    if (!(identical(params[[parameterName]],TRUE)|
                identical(params[[parameterName]],FALSE))) {
        stop("The parameter ", parameterName, " must be logical, \n",
               "NAs are not allowed.")
    }
    return(params[[parameterName]])
}

###############################################################################

##function, that tests a param, which has exact one single character string,
##if it is of a set of possible values. If yes, the param returns.
##Furthermore, a default-value is setted if necessary.
checkSingleValParams <- function(parameterName, params,
        defaultValue, possibleValues){
    ##default-value
    if (is.null(params[[parameterName]])) {
        params[[parameterName]] <- defaultValue
    } else {
        if (length(params[[parameterName]]) != 1) {
            stop("The parameter ", parameterName,
                            " only can have one value!")
        }
        params[[parameterName]] <- checkIsValue(parameterName,
                params, possibleValues)
        params[[parameterName]] <- tolower(params[[parameterName]])
    }
    return(params[[parameterName]])
}

###############################################################################

##function, that tests a param, which has exact one single character string,
##if it is of a set of possible values. If yes, the param returns.
##No default-value!!!
checkSingleValParamsNew <- function(parameterName,
        params,
        possibleValues){

    if (!is.null(params[[parameterName]])) {
        if (length(params[[parameterName]]) != 1) {
            stop("The parameter ", parameterName,
                            " only can have one value!")
        }
        params[[parameterName]] <- checkIsValue(parameterName,
                params, possibleValues)
        params[[parameterName]] <- tolower(params[[parameterName]])
    }
    return(params[[parameterName]])
}

###############################################################################

##function, that tests a param, if it has exact one single character string,

checkString <- function(parameterName, params){

    if (!is.null(params[[parameterName]])) {
        if (length(params[[parameterName]]) != 1) {
            stop("The parameter ", parameterName,
                            " demands a single string!")
        }
        if (!is.character(params[[parameterName]])) {
            stop("The parameter ", parameterName,
                            " demands a single string!")
        }
    }
    return(params[[parameterName]])
}

###############################################################################

##function, that tests a param, which has as input any vector of type
##c("","",...), or list("","",...), whether all inputs are of a set of
##possible values. If yes, the vector returns.
checkValueParams <- function(parameterName, params, possibleValues){
    for (i in 1: length(params[[parameterName]])) {
        ##check if is params$parameterName of type character
        if(!is.character(params[[parameterName]][[i]])) {
            stop("The parameter ", parameterName,
                            " should contain strings!")
        }

        ##check, if input of parameter is valid
        if (!(tolower(params[[parameterName]])[[i]] %in% possibleValues)){
            ##create a string with all possible Values named text
            text <- ""
            text <- paste(possibleValues, collapse=", ")
            stop("The parameter ", parameterName,
                   " only can have the values: \n", text,
                   "\n Check, whether there are blanks or typos in between!")
        }
    }
    return(tolower(params[[parameterName]]))
}

###############################################################################

##function, that tests, whether an input of a parameter is an Integer or not;
##sets default-value if necessary and returns the parameter
checkIntegerParams <- function(parameterName, params, defaultValue) {
    ##default-value
    if (is.null(params[[parameterName]])) {
        params[[parameterName]] <- as.integer(defaultValue)
    }

    ##check, if input of parameter is valid
    if (!is.integer(params[[parameterName]])) {
        if (is.numeric(params[[parameterName]])) {
            if (is.matrix(params[[parameterName]])) {
                stop("The parameter ", parameterName,
                                " should be an integer, not a matrix!")
            }
            if (length(params[[parameterName]]) != 1) {
                stop("The parameter ", parameterName,
                                " should be an integer, not a vector!")
            }
            ##stop if usage of floats
            if (params[[parameterName]] -
                round(params[[parameterName]]) != 0) {
                    stop("The parameter ", parameterName,
                           " should be an integer, not numeric!")
            }
            if (params[[parameterName]] <= .Machine$integer.max) {
                params[[parameterName]] <- as.integer(params[[parameterName]])
            } else {
                stop("The parameter ", parameterName,
                       " is bigger than an integer!")
            }
        } else {
            stop("The parameter ", parameterName,
                   " should be an integer or at least numeric!")
        }
    } else {
        if (is.matrix(params[[parameterName]])) {
            stop("The parameter ", parameterName,
                   " should be an integer, not a matrix!")
        }
        if (length(params[[parameterName]]) != 1) {
            stop("The parameter ", parameterName,
                   " should be an integer, not a vector!")
        }
    }
    return(params[[parameterName]])
}

###############################################################################

##function, that tests, whether an input of a parameter is an Integer or not;
##sets NO DEFAULT-VALUE and returns the parameter
checkIntegerParamsNew <- function(parameterName, params) {
    if (!is.null(params[[parameterName]])) {
        ##check, if input of parameter is valid
        if (!is.integer(params[[parameterName]])) {
            if (is.numeric(params[[parameterName]])) {
                if (is.matrix(params[[parameterName]])) {
                    stop("The parameter ", parameterName,
                                    " should be an integer, not a matrix!")
                }
                if (length(params[[parameterName]]) != 1) {
                    stop("The parameter ", parameterName,
                                    " should be an integer, not a vector!")
                }
                ##stop if usage of floats
                if (params[[parameterName]] -
                    round(params[[parameterName]]) != 0) {
                        stop("The parameter ", parameterName,
                               " should be an integer, not numeric!")
                }
                if (params[[parameterName]] <= .Machine$integer.max) {
                    params[[parameterName]] <- as.integer(
                                                  params[[parameterName]])
                } else {
                    stop("The parameter ", parameterName,
                           " is bigger than an integer!")
                }
            } else {
                stop("The parameter ", parameterName,
                       " should be an integer or at least numeric!")
            }
        } else {
            if (is.matrix(params[[parameterName]])) {
                stop("The parameter ", parameterName,
                       " should be an integer, not a matrix!")
            }
            if (length(params[[parameterName]]) != 1) {
                stop("The parameter ", parameterName,
                       " should be an integer, not a vector!")
            }
        }
    }
    return(params[[parameterName]])
}
###############################################################################

##function, that tests the param if it is positive
checkPositiveParams <- function(parameterName, params){
    if (!is.null(params[[parameterName]])) {
        if (is.numeric(params[[parameterName]])) {
            if (is.matrix(params[[parameterName]])) {
                stop("The parameter ", parameterName,
                       " should be a positive value, not a matrix!")
            }
            if (length(params[[parameterName]]) != 1) {
                stop("The parameter ", parameterName,
                       " should be a positive value, not a vector!")
            }
            if (params[[parameterName]] < 0) {
                stop("The parameter ", parameterName, " should be positive!")
            }
        } else {
            stop("The parameter ", parameterName,
                   " should be a positive numeric!")
        }
    }
    return(params[[parameterName]])
}

###############################################################################

##function, that tests the param if it is negative
checkNegativeParams <- function(parameterName, params){

    if (!is.null(params[[parameterName]])) {
        if (is.numeric(params[[parameterName]])) {
            if (is.matrix(params[[parameterName]])) {
                stop("The parameter ", parameterName,
                       " should be a negative value, not a matrix!")
            }
            if (length(params[[parameterName]]) != 1) {
                stop("The parameter ", parameterName,
                       " should be a negative value, not a vector!")
            }
            if (params[[parameterName]] > 0) {
                stop("The parameter ", parameterName, " should be negative!")
            }
        } else {
            stop("The parameter ", parameterName,
                   " should be a negative numeric!")
        }
    }
    return(params[[parameterName]])
}

###############################################################################

##function, that tests, whether an input of a parameter is numeric or not;
##sets default-value if necessary and returns the parameter
checkNumericParams <- function(parameterName, params, defaultValue) {
    ##default-value
    if (is.null(params[[parameterName]])) {
        params[[parameterName]] <- defaultValue
    }

    ##check, if input of parameter is valid
    if (is.numeric(params[[parameterName]])) {
        if (is.matrix(params[[parameterName]])) {
            stop("The parameter ", parameterName,
                   " should be numeric, not a matrix!")
        }
        if (length(params[[parameterName]]) != 1) {
            stop("The parameter ", parameterName,
                   " should be numeric, not a vector!")
        }
    } else {
        stop("The parameter ", parameterName, " should be numeric!")
    }
    return(params[[parameterName]])
}

###############################################################################

##function, that tests, whether an input of a parameter is numeric or not;
##sets NO DEFAULT-VALUE and returns the parameter
checkNumericParamsNew <- function(parameterName, params) {
    if (!is.null(params[[parameterName]])) {
        ##check, if input of parameter is valid
        if (is.numeric(params[[parameterName]])) {
            if (length(params[[parameterName]]) != 1) {
                stop("The parameter ", parameterName,
                                " should be numeric, not a vector!")
            }
            if (is.nan(params[[parameterName]])) {
                stop("The parameter ", parameterName,
                       " should be numeric, NaN is not allowed!")
            }
            if (is.matrix(params[[parameterName]])) {
                stop("The parameter ", parameterName,
                       " should be numeric, not a matrix!")
            }
        } else {
            stop("The parameter ", parameterName, " should be numeric!")
        }
    }
    return(params[[parameterName]])
}

###############################################################################

##function, which evaluates, whether a parameter value is in between an
##interval or not; sets default-value if necessary and returns the parameter
checkIntervalParams <- function(parameterName,
                                params,
                                defaultValue,
                                lowerB,
                                upperB){
    ##default-value
    if (is.null(params[[parameterName]])) {
        params[[parameterName]] <- defaultValue
    }

    ##check, if input of filter is valid
    if (is.numeric(params[[parameterName]])) {
        ##check if filter is negative
        if (params[[parameterName]] < lowerB |
                params[[parameterName]] > upperB) {
            stop("The parameter ", parameterName,
                   " should be in the interval [",
                   lowerB, ",", upperB, "]!")
        }
    } else {
        stop("The parameter ", parameterName,
               " should be a numeric in [",
               lowerB, ",", upperB, "]!")
    }
    return(params[[parameterName]])
}

###############################################################################

##function, which evaluates, whether a parameter value is in between an
##interval or not; sets NO DEFAULT-VALUE and returns the parameter
checkIntervalParamsNew <- function(parameterName,
                                   params,
                                   lowerB,
                                   upperB){

    if (!is.null(params[[parameterName]])) {
        ##check, if input of filter is valid
        if (is.numeric(params[[parameterName]])) {
            ##check if filter is negative
            if (params[[parameterName]] < lowerB |
                    params[[parameterName]] > upperB) {
                stop("The parameter ", parameterName,
                                " should be in the interval [",
                                lowerB, ",", upperB, "]!")
            }
        } else {
            stop("The parameter ", parameterName,
                   " should be a numeric in [",
                   lowerB, ",", upperB, "]!")
        }
    }
    return(params[[parameterName]])
}


###############################################################################
##function that checks, whether a input file exists or not
##-throw exception if not, wrong directory or empty file

checkInFile <- function(parameterName, params){
    if (!is.character(params[[parameterName]]) ||
            length(params[[parameterName]]) != 1) {
        stop("The parameter ", parameterName,
               " must be single character string!")
    }
    if (!file.exists(params[[parameterName]])){
        stop("The file for parameter ", parameterName ," does not exist!")
    }
    if (file.info(params[[parameterName]])$size == 0){
        stop("The file for parameter ", parameterName ," is empty!")
    }
}

###############################################################################
##function that checks, whether a output file exists or not
##-create directory if not
##-returns list with 2 params
## 1. checked file path
## 2. flag if file exists

checkOutFile <- function(parameterName, params){
    result <- list()
    if (!is.character(params[[parameterName]]) ||
        length(params[[parameterName]]) != 1) {
        stop("The parameter ", parameterName,
               " must be single character string")
    }
    if (file.exists(params[[parameterName]])){
        result[["existingFile"]] <- TRUE
        result[["param"]] <- params[[parameterName]]
    } else {
        result[["existingFile"]] <- FALSE
        result[["param"]] <- params[[parameterName]]
    }
    return(result)
}

###############################################################################
##function for a profile score check, whether le, sp, sv or spn are used
##returns list with 4 parameters, all boolean:
##result$le
##result$sp
##result$sv
##result$spn

checkProfileScore <- function(type, params){
    result <- list()
    ##defaultValues, if all 4 parameters (le, sp, sv, spn) are NULL
    if(is.null(params[["le"]]) && is.null(params[["sp"]]) &&
       is.null(params[["sv"]]) && is.null(params[["spn"]])){
        if (identical(type, "protein")) {
            result[["le"]] <- TRUE
            result[["sp"]] <- FALSE
            result[["sv"]] <- FALSE
            result[["spn"]] <- FALSE
        } else if (identical(type, "rna") || identical(type, "dna")) {
            result[["le"]] <- FALSE
            result[["sp"]] <- FALSE
            result[["sv"]] <- FALSE
            result[["spn"]] <- TRUE
        }
        ##check, if all are boolean
        ##if any of the parameters is NULL, the default-Value is set
    } else {
        if (identical(type, "protein")) {
            params[["sp"]] <- checkLogicalParams("sp", params, FALSE)
            ##if sp==TRUE set le=FALSE
            if (params[["sp"]]) {
                params[["le"]] <- FALSE
            }
            params[["sv"]] <- checkLogicalParams("sv", params, FALSE)
            ##if sv==TRUE set le=FALSE
            if (params[["sv"]]) {
                params[["le"]] <- FALSE
            }
            params[["le"]] <- checkLogicalParams("le", params, TRUE)
            params[["spn"]] <- checkLogicalParams("spn", params, FALSE)
        } else {
            params[["spn"]] <- checkLogicalParams("spn", params, TRUE)
            params[["le"]] <- checkLogicalParams("le", params, FALSE)
            params[["sp"]] <- checkLogicalParams("sp", params, FALSE)
            params[["sv"]] <- checkLogicalParams("sv", params, FALSE)
        }
        ##consistency check
        ##type==RNA|DNA =>only spn==TRUE, all others FALSE possible
        if (identical(type, "rna") || identical(type, "dna")) {
            if (!params[["spn"]] | params[["le"]] | params[["sp"]] |
                 params[["sv"]]){
                stop("The used profile score is inconsistent. \n",
                     "If you use nucleotides, ",
                     "the parameter spn should be TRUE! \n",
                     "All others (sp, sv, le) should be FALSE!")
            }
        }
        ##type==protein =>only spn==FALSE possible
        if (identical(type, "protein")){
            if (params[["spn"]]) {
                stop("The used profile score is inconsistent. \n",
                     "If you use proteins, ",
                     "the parameter spn should be FALSE!")
            }
            ##type==protein =>only 1 of the others (sp, sv, le) TRUE
            if ((params[["sv"]] && params[["le"]]) ||
                (params[["sv"]] && params[["sp"]]) ||
                (params[["sp"]] && params[["le"]]) ||
                (params[["sp"]] && params[["le"]] && params[["sv"]]))
            {
                stop("The used profile score is inconsistent. \n",
                     "Only one of the parameter sp, sv, le can be TRUE!")
            }
        }

        ##all 4 are negative
        if (!params[["spn"]] && !params[["sp"]] &&
            !params[["sv"]] && !params[["le"]]) {
            stop("The used profile score is inconsistent. \n",
                 "You are not allowed to set all 4 possibilities FALSE!")
        }
        result[["le"]] <- params[["le"]]
        result[["sp"]] <- params[["sp"]]
        result[["sv"]] <- params[["sv"]]
        result[["spn"]] <- params[["spn"]]
    }
    return(result)
}

###############################################################################
##consistency check for a profile score, whether le, sp, sv or spn are used
##stops, if any inconsistency appears

checkProfileScoreNew <- function(type, params){

    if (identical(type, "protein")) {
        ##type==protein =>only spn=FALSE possible
        if (params[["spn"]]) {
            stop("The used profile score is inconsistent. \n",
                 "If you use proteins, the prameter spn should be FALSE!")
        }
        ##type==protein =>only 1 of the others (sp, sv, le) TRUE
        if ((params[["sv"]] && params[["le"]]) ||
            (params[["sv"]] && params[["sp"]]) ||
            (params[["sp"]] && params[["le"]]) ||
            (params[["sp"]] && params[["le"]] && params[["sv"]])){
            stop("The used profile score is inconsistent. \n",
                 "Only one of the parameters sp, sv, le can be TRUE!")
        }
    } else {
        ##consistency check
        ##type==RNA|DNA =>only spn=TRUE, all others FALSE possible
        if (params[["le"]] | params[["sp"]] |params[["sv"]]){
            stop("The used profile score is inconsistent. \n",
                "If you use nucleotides, the parameter spn should be TRUE! \n",
                "All others (sp, sv, le) should be FALSE!")
        }
    }
}

###############################################################################
checkFunctionAvailable <- function(name) {
    #mylibs <- library.dynam()
    #hasFunction <- FALSE
    #for (i in 1:length(mylibs)) {
    #    cur <- mylibs[[i]]
    #    if (identical(name, cur[[1]])) {
    #        hasFunction <- TRUE
    #        return(hasFunction)
    #    }
    #}
    #return(hasFunction)
    return(TRUE)
}
