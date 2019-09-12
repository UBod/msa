msaPrettyPrint <- function(x, y, output=c("pdf", "tex", "dvi", "asis"),
                           subset=NULL, file=NULL, alFile=NULL,
                           askForOverwrite=TRUE, psFonts=FALSE, code=NA,
                           paperWidth=11, paperHeight=8.5, margins=c(0.1, 0.3),
                           shadingMode=c("identical", "similar", "functional"),
                           shadingModeArg=NA,
                           shadingColors=c("blues", "reds", "greens", "grays",
                                           "black"),
                           showConsensus=c("bottom", "top", "none"),
                           consensusColors=c("ColdHot", "HotCold", "BlueRed",
                                             "RedBlue", "GreenRed", "RedGreen",
                                             "Gray"),
                           consensusThreshold=50,
                           showLogo=c("top", "bottom", "none"),
                           logoColors=c("chemical", "rasmol", "hydropathy",
                                        "structure", "standard area",
                                        "accessible area"),
                           showLogoScale=c("none", "leftright",
                                           "left", "right"),
                           showNames=c("left", "right", "none"),
                           showNumbering=c("right", "left", "none"),
                           showLegend=TRUE, furtherCode=NA, verbose=FALSE)
{
    xname <- deparse(substitute(x))
    output <- match.arg(output)

    if (is.null(file))
    {
        if (length(grep("[^\\w]", xname, perl=TRUE)) > 0)
        {
            warning("Cannot use default file name '", xname, ".", output,
                    "' because it would contain invalid characters => ",
                    "resorting to 'msaPrettyPrintOutput.", output, "'!")
            xname <- "msaPrettyPrintOutput"
        }

        file <- paste(xname, output, sep=".")
    }

    if (is.null(alFile))
        alFile <- tempfile(pattern="seq", tmpdir=tempdir(), fileext=".fasta")
    else if (is.character(alFile) &&
             substr(alFile, nchar(alFile) - 5, nchar(alFile)) == ".fasta")
    {
        if (askForOverwrite && file.exists(alFile))
        {
            message("File ", alFile, " exists. Overwrite? (y/N)")

            answer <- try(tolower(scan(what=character(), nmax=1,
                                       quiet=TRUE)), silent=TRUE)

            if (nchar(answer) != 1 || substr(answer, 1, 1) != "y")
                return(invisible(NULL))
        }
    }
    else
        stop("The parameter alFile has an invalid argument!")

    if (!is(x, "MultipleAlignment"))
        stop("The parameter x has an invalid argument! \n",
             "x must be a multiple alignment object!")

    if (output != "asis")
    {
        if (!is.numeric(paperWidth) || length(paperWidth) != 1 ||
            paperWidth <= 0)
            stop("The parameter paperWidth must be ",
                 "single positive number (unit: inches)!")

        if (!is.numeric(paperHeight) || length(paperHeight) != 1 ||
            paperHeight <= 0)
            stop("The parameter paperHeight must be ",
                 "single positive number (unit: inches)!")

        if (!is.numeric(margins) || length(margins) != 2)
            stop("The parameter margins must be ",
                 "two positive numbers (unit: inches)!")
    }

    if (!identical(subset, NULL) && !identical(subset, NA))
    {
        if (is.numeric(subset))
        {
            if (max(subset) < .Machine$integer.max)
                subset <- as.integer(subset)
            else
               stop("One or more values for parameter subset ",
                    "are larger than integer!")
        }
        else if (!is.integer(subset))
            stop("The parameter subset has an invalid argument!")

        if (length(subset) < 2)
            stop("The parameter subset is expected to be \n",
                 " a vector with at least 2 entries!")

        if (!all(subset %in% 1:nrow(x)))
            stop("Some values in parameter subset are out of range!")
    }
    else if (length(rowmask(x)) > 0)
    {
        if (setdiff(IRanges(start=1, end=nrow(x)), rowmask(x))
            <= .Machine$integer.max)
            subset <- as.integer(setdiff(IRanges(start=1, end=nrow(x)),
                                         rowmask(x)))
        else
            stop("There is no typecast possible in parameter subset!")
    }
    else
        subset <- 1:nrow(x)

    shadingMode <- match.arg(shadingMode)
    shadingColors <- match.arg(shadingColors)
    showConsensus <- match.arg(showConsensus)
    consensusColors <- match.arg(consensusColors)
    showLogo <- match.arg(showLogo)
    logoColors <- match.arg(logoColors)
    showLogoScale <- match.arg(showLogoScale)
    showNames <- match.arg(showNames)
    showNumbering <- match.arg(showNumbering)

    if (!is.numeric(consensusThreshold) || length(consensusThreshold) < 1 ||
        length(consensusThreshold) > 2 ||
        any(consensusThreshold < 0) || any(consensusThreshold > 100))
        stop("The parameter consensusThreshold must be \n",
             "one or two numbers between 0 and 100 !")
    else if (length(consensusThreshold) == 2 &&
             consensusThreshold[1] >= consensusThreshold[2])
        stop("The second percentage in consensusThreshold must be \n",
             "at least as large as the first one!")

    if (shadingMode %in% c("identical", "similar"))
    {
        if (!identical(shadingModeArg, NA) &&
                (!is.numeric(shadingModeArg) ||
                 length(shadingModeArg) > 2 ||
                 length(shadingModeArg) < 1 ||
                 (length(shadingModeArg) == 2 &&
                     shadingModeArg[1] > shadingModeArg[2])||
                 shadingModeArg[1] < 0 ||
                 shadingModeArg[1] > 100 ||
                 (length(shadingModeArg) == 2 &&
                     (shadingModeArg[2] < 0 ||
                      shadingModeArg[2] > 100))))
           stop("If identical or similarity shading is used, shadingModeArg\n",
                 "must be a single numeric threshold between 0 and 100 or\n",
                 "two thresholds between 0 and 100 in increasing order!")
    }
    else if (identical(shadingMode, "functional"))
    {
        if (!identical(shadingModeArg, NA))
            shadingModeArg <- match.arg(shadingModeArg,
                                       c("charge", "hydropathy", "structure",
                                         "chemical", "rasmol", "standard area",
                                         "accessible area"))
        else
            stop("Missing shadingModeArg for functional shading mode. \n",
                 "Valid values are: \n",
                 "\"charge\", \n",
                 "\"hydropathy\", \n",
                 "\"structure\", \n",
                 "\"chemical\",\n",
                 " \"rasmol\",\n",
                 "\"standard area\",\n",
                 "\"accessible area\"!")
    }
    else if (!identical(shadingMode, NA))
        stop("The parameter shadingModeArg has an invalid argument!")

    if (showConsensus != "none" && showConsensus == showLogo)
        stop("Cannot display consensus sequence and sequence logo ",
             "on the same side!")

    if (showNames != "none" && showNames == showNumbering)
        stop("Cannot display sequence names and numbering on the same side!")

    if (!identical(code, NA) && !is.character(code))
        stop("The parameter code has an invalid argument!")

    if (!identical(furtherCode, NA) && !is.character(furtherCode))
        stop("The parameter furtherCode has an invalid argument!")

    if (missing(y))
        toShow <- IRanges(start=1, end=ncol(x))
    else if (is(y, "IRanges"))
    {
        if (all(start(y) >= 1) && all(end(y) <= ncol(x)))
            toShow <- reduce(y)
        else
            stop("The parameter y has invalid ranges: out of bounds!")
    }
    else if ((is.numeric(y) || is.integer(y)) && length(y) == 2 && y[1] >= 1 &&
             y[2] <= ncol(x) && y[1] < y[2])
        toShow <- IRanges(start=y[1], end=y[2])
    else
        stop("The parameter y has an invalid argument!")

    if (length(colmask(x)) > 0)
        toShow <- setdiff(toShow, colmask(x))

    if (sum(width(toShow)) == 0)
        stop("Sequences empty or everything masked: nothing to be shown!")

    jobname <- ""
    suffix <- ""

    if (output != "asis")
    {
        if (!is.character(file) || length(file) > 1)
            stop("The argument for parameter file must be \n",
                 "a single character string!")

        if (substr(file, nchar(file) - 2, nchar(file)) != output)
            stop("The file name suffix and output type do not match!")

        jobname <- substr(file, 1, nchar(file) - 4)

        if (length(grep("[^\\w\\-/\\\\:.]", jobname, perl=TRUE)) > 0)
        {
            warning("Cannot use file name '", file,
                    "' because it contains invalid characters => \n",
                    "resorting to 'msaPrettyPrintOutput.", output, "'!")
            jobname <- "msaPrettyPrintOutput"
            file <- paste0(jobname, output)
        }

        if (askForOverwrite && file.exists(file))
        {
            message("File ", file, " exists. Overwrite? (y/N)")

            answer <- try(tolower(scan(what=character(), nmax=1,
                                       quiet=TRUE)), silent=TRUE)

            if (nchar(answer) != 1 || substr(answer, 1, 1) != "y")
                return(invisible(NULL))
        }
    }

    writeXStringSet(as(unmasked(x), "XStringSet")[subset], filepath=alFile)

    if (verbose)
        message("Multiple alignment written to temporary file ", alFile)

    texOutput <- paste0("\\begin{texshade}{", stratifyFilenames(alFile), "}")

    if (is(x, "AAMultipleAlignment"))
        texOutput <- c(texOutput, "\\seqtype{P}")
    else
        texOutput <- c(texOutput, "\\seqtype{N}")

    if (length(toShow) == 1)
    {
        if (sum(width(toShow)) < ncol(x))
            texOutput <- c(texOutput, paste("\\setends{consensus}{",
                                            start(toShow), "..", end(toShow),
                                            "}", sep=""))
    }
    else
    {
        coList <- sapply(1:length(toShow),
                         function(i) paste(start(toShow)[i], "..",
                                           end(toShow)[i], sep=""))

        texOutput <- c(texOutput, paste("\\setdomain{consensus}{",
                                        paste(coList, collapse=","), "}",
                                        sep=""))
    }

    if (identical(code, NA))
    {
        if (identical(shadingModeArg, NA))
            texOutput <- c(texOutput,
                               paste("\\shadingmode{", shadingMode, "}",
                                     sep=""))
        else
            texOutput <- c(texOutput,
                               paste("\\shadingmode[",
                                     shadingModeArg, "]{",
                                     shadingMode, "}", sep=""))

        if (length(consensusThreshold) == 2)
            texOutput <- c(texOutput, paste("\\threshold[",
                                            consensusThreshold[2], "]{",
                                            consensusThreshold[1], "}",
                                            sep=""))
        else
            texOutput <- c(texOutput, paste("\\threshold{",
                                            consensusThreshold[1], "}",
                                            sep=""))

        if (showConsensus != "none")
        {
            texOutput <- c(texOutput,
                           paste("\\showconsensus[", consensusColors,
                                 "]{", showConsensus, "}", sep=""))
        }
        else
            texOutput <- c(texOutput, "\\hideconsensus")

        texOutput <- c(texOutput, paste("\\shadingcolors{",
                                        shadingColors, "}", sep=""))

        if (showLogo != "none")
            texOutput <- c(texOutput,
                           paste("\\showsequencelogo[", logoColors,
                                 "]{", showLogo, "}", sep=""))

        if (showLogoScale == "none")
            texOutput <- c(texOutput, "\\hidelogoscale")
        else
            texOutput <- c(texOutput,
                           paste("\\showlogoscale{", showLogoScale, "}",
                                 sep=""))

        if (showNames != "none")
        {
            seqNames <- rownames(x)[subset]
            pattern <- "[^a-zA-Z0-9,;:.?!/\\-\\(\\)\\'\" ]"
            seqNames <- gsub(pattern, " ", seqNames, perl=TRUE)

            texOutput <- c(texOutput,
                           paste("\\shownames{", showNames, "}", sep=""),
                           paste("\\nameseq{", 1:length(subset), "}{",
                                 seqNames, "}", sep=""))
        }
        else
            texOutput <- c(texOutput, "\\hidenames")

        if (showNumbering != "none")
            texOutput <- c(texOutput,
                           paste("\\shownumbering{", showNumbering, "}",
                                 sep=""))
        else
            texOutput <- c(texOutput, "\\hidenumbering")

        if (showLegend)
            texOutput <- c(texOutput, "\\showlegend")

        if (!identical(furtherCode, NA))
            texOutput <- c(texOutput, furtherCode)
    }
    else
        texOutput <- c(texOutput, code)

    texOutput <- c(texOutput, "\\end{texshade}")

    if (output == "asis")
        cat(texOutput, sep="\n")
    else
    {
        texHeader <- c("\\documentclass[10pt]{article}", "")

        if (psFonts)
            texHeader <- c(texHeader, "\\usepackage{times}")

        texHeader <- c(texHeader, "\\usepackage{texshade}")

        texHeader <- c(texHeader, "", "\\headheight=0pt", "\\headsep=0pt",
                       "\\hoffset=0pt", "\\voffset=0pt",
                       paste0("\\paperwidth=", paperWidth, "in"),
                       paste0("\\paperheight=", paperHeight, "in"),
                       "\\ifx\\pdfoutput\\undefined",
                       "\\relax",
                       "\\else",
                       "\\pdfpagewidth=\\paperwidth",
                       "\\pdfpageheight=\\paperheight",
                       "\\fi",
                       paste0("\\oddsidemargin=", margins[1] - 1, "in"),
                       paste0("\\topmargin=", margins[2] - 1, "in"),
                       paste0("\\textwidth=",
                              paperWidth - 2 * margins[1], "in"),
                       paste0("\\textheight=",
                              paperHeight - 2 * margins[2],"in"),
                       "", "\\pagestyle{empty}", "", "\\begin{document}")
        texFooter <- "\\end{document}"

        if (output == "tex")
            writeLines(c(texHeader, texOutput, texFooter), con=file)
        else
        {
            texfile <- paste(jobname, "tex", sep=".")

            if (askForOverwrite && file.exists(texfile))
            {
                message("File ", texfile, " exists. Overwrite? (y/N)")

                answer <- try(tolower(scan(what=character(), nmax=1,
                                           quiet=TRUE)), silent=TRUE)

                if (nchar(answer) != 1 || substr(answer, 1, 1) != "y")
                    return(invisible(NULL))
            }

            writeLines(c(texHeader, texOutput, texFooter), con=texfile)

            if (verbose)
                message("File ", texfile, " created")

            texi2dvi(texfile, quiet=!verbose, pdf=identical(output, "pdf"),
                     texinputs=system.file("tex", package="msa"),
                     clean=TRUE, index=FALSE)
        }

        if (verbose)
            message("Output file ", file, " created")
    }

    invisible(texOutput)
}
