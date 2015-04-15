/*********************************************************************
 * Clustal Omega - Multiple sequence alignment
 *
 * Copyright (C) 2010 University College Dublin
 *
 * Clustal-Omega is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This file is part of Clustal-Omega.
 *
 ********************************************************************/

/*
 *  RCS $Id: main.cpp 234 2011-04-13 05:26:16Z andreas $
 */

/*
 * We are using a mix of C and C++, which means that linking has to be
 * done with a C++ compiler. By using this "fake" main c++ function,
 * automake is convinced to use a C++ compiler for linking.
 *
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "RClustalOmega.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

extern "C" {
#include "mymain.h"
#include "clustal/util.h"
#include "squid/squid.h"
}

using namespace std;
using namespace Rcpp;

bool hasClustalOmegaEntry(Rcpp::List params, const char* entry) {
	Rcpp::CharacterVector namesCV = params.names();
	int n = namesCV.size();
	vector<string> namesVector = Rcpp::as< std::vector< std::string > >(namesCV);
	for (int i = 0; i < n; i++) {
		const char* vec = namesVector[i].c_str();
		if (strcmp(vec, entry) == 0) {
			return !Rf_isNull(params[entry]);
		}
	}
	return false;
}

void appendDoubleToString(SEXP value, const char *str, char *buffer) {
    double doubleValue = as<double>(value);
    sprintf(buffer, "%s=%f", str, doubleValue);
}

void appendIntToString(SEXP value, const char *str, char *buffer) {
    int intValue = as<int>(value);
    sprintf(buffer, "%s=%i", str, intValue);
}

void appendStringToString(SEXP value, const char *str, char *buffer) {
    string stringValue = as<string>(value);
    sprintf(buffer, "%s=%s", str, stringValue.c_str());
}

void appendString(char ***argv, int &count, const char *str) {
    count++;
    *argv = (char**)realloc(*argv, count * sizeof(char*));
    if(*argv == NULL)
        Rprintf("Error (re)allocating memory\n");
    else {
        (*argv)[count - 1] = (char*)malloc(strlen(str)*sizeof(char) + 1);
        strcpy((*argv)[count - 1], str);
    }
}

void appendIntValue(char ***argv, int &count, const char *str, SEXP intValue) {
    if (!intValue || Rf_isNull(intValue)) {
        return;
    }
    char buffer[128];
    appendIntToString(intValue, str, buffer);
    appendString(argv, count, buffer);
}

void appendDoubleValue(char ***argv, int &count, const char *str, SEXP doubleValue) {
    if (!doubleValue || Rf_isNull(doubleValue)) {
        return;
    }
    char buffer[128];
    appendDoubleToString(doubleValue, str, buffer);
    appendString(argv, count, buffer);
}

void appendStringValue(char ***argv, int &count, const char *str, SEXP stringValue) {
    if (!stringValue || Rf_isNull(stringValue)) {
        return;
    }
    char buffer[128];
    appendStringToString(stringValue, str, buffer);
    appendString(argv, count, buffer);
}

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str) {
    List lst(list);
    SEXP element = NULL;
    CharacterVector nam = lst.names();    
    for (int i = 0; i < nam.size(); i++) {
        string nameStr = as<string>(nam[i]);
        if(strcmp(nameStr.c_str(), str) == 0) {
           element = lst[str];
           break;
        }
    }
    return element;
}

char* getChar(string str) {
    char * writable = new char[str.size() + 1];
    std::copy(str.begin(), str.end(), writable);
    writable[str.size()] = '\0';
    return writable;
}

SEXP RClustalOmega(SEXP rInputSeqs,
                   SEXP rCluster,
                   SEXP rGapOpening,
                   SEXP rGapExtension,
                   SEXP rMaxiters,
                   SEXP rSubstitutionMatrix,
                   SEXP rType,
                   SEXP rVerbose,
                   SEXP rParams) {

	Rcpp::List retList;
    try {
        Rcpp::List params(rParams);

        bool inputSeqIsFileFlag = as<bool>(params["inputSeqIsFileFlag"]);

        ClustalOmegaInput msaInput;

        int argc = 0;
        char **argv = NULL;

        appendString(&argv, argc, "clustalo");
        R_len_t n = Rf_length(rInputSeqs);
        if (inputSeqIsFileFlag) {
        	appendString(&argv, argc, "-i");  //input
            string inputSeqFile = as<string>(rInputSeqs);
            appendString(&argv, argc, inputSeqFile.c_str());
        } else {
            vector<string> myInputSeqs = Rcpp::as<vector<string> >(rInputSeqs);
            CharacterVector myInputSeqsCV(rInputSeqs);
            vector<string> mySeqNames =
                    Rcpp::as<vector<string> >(myInputSeqsCV.attr("names"));

            R_len_t i;
            char **seqs = (char **)R_alloc(n, sizeof(char *));
            char **seqNames = (char **)R_alloc(n, sizeof(char *));
            int *seqLengths = (int *)R_alloc(n, sizeof(int));
            int seqLength = n;

            for (i = 0; i < n; i++) {
                seqs[i] = getChar(myInputSeqs[i]);
                seqNames[i] = getChar(mySeqNames[i]);
                seqLengths[i] = strlen(seqs[i]);
            }
            msaInput.inputSeqs = seqs;
            msaInput.seqLength = seqLength;
            msaInput.seqNames = seqNames;
            //hardcoded value, if sequence is an R sequence
            appendString(&argv, argc, "--R");
        }
        appendString(&argv, argc, "-o");
        //used in some cases ...
        appendString(&argv, argc, "tempClustalOmega.aln");
        //output format, we expect
        appendString(&argv, argc, "--outfmt=clustal");
        /*FIXME TODO in later version: activate other formats
        appendStringValue(&argv, argc, "--outfmt",
                getListElement(rParams, "outFmt"));*/
        appendStringValue(&argv, argc, "--infmt",
                getListElement(rParams, "inFmt"));
        appendStringValue(&argv, argc, "--seqtype", rType);
        //overwrite files ...
        appendString(&argv, argc, "--force");
        appendDoubleValue(&argv, argc, "--gapopen", rGapOpening);
        appendDoubleValue(&argv, argc, "--gapext", rGapExtension);
        appendIntValue(&argv, argc, "--cluster-size", rCluster);
        appendIntValue(&argv, argc, "--iter", rMaxiters);

        if (hasClustalOmegaEntry(rParams, "clusteringOut")) {
			appendStringValue(&argv, argc, "--clustering-out",
					getListElement(rParams, "clusteringOut"));
        }
        if (hasClustalOmegaEntry(rParams, "macRam")) {
			appendIntValue(&argv, argc, "--MAC-RAM",
					getListElement(rParams, "macRam"));
        }
        if (hasClustalOmegaEntry(rParams, "maxGuidetreeIterations")) {
			appendIntValue(&argv, argc, "--max-guidetree-iterations",
					getListElement(rParams, "maxGuidetreeIterations"));
        }
        if (hasClustalOmegaEntry(rParams, "maxHmmIterations")) {
			appendIntValue(&argv, argc, "--max-hmm-iterations",
					getListElement(rParams, "maxHmmIterations"));
        }
        if (hasClustalOmegaEntry(rParams, "maxNumSeq")) {
			appendIntValue(&argv, argc, "--maxnumseq",
					getListElement(rParams, "maxNumSeq"));
        }
        if (hasClustalOmegaEntry(rParams, "maxSeqLen")) {
			appendIntValue(&argv, argc, "--maxseqlen",
					getListElement(rParams, "maxSeqLen"));
        }
        if (hasClustalOmegaEntry(rParams, "threads")) {
			appendIntValue(&argv, argc, "--threads",
					getListElement(rParams, "threads"));
        }
        if (hasClustalOmegaEntry(rParams, "wrap")) {
			appendIntValue(&argv, argc, "--wrap",
					getListElement(rParams, "wrap"));
        }
        if (hasClustalOmegaEntry(rParams, "distMatIn")) {
			appendStringValue(&argv, argc, "--distmat-in",
					getListElement(rParams, "distMatIn"));
        }
        if (hasClustalOmegaEntry(rParams, "distMatOut")) {
			appendStringValue(&argv, argc, "--distmat-out",
					getListElement(rParams, "distMatOut"));
        }
        if (hasClustalOmegaEntry(rParams, "outputOrder")) {
			appendStringValue(&argv, argc, "--output-order",
					getListElement(rParams, "outputOrder"));
        }
        if (hasClustalOmegaEntry(rParams, "guideTreeIn")) {
			appendStringValue(&argv, argc, "--guidetree-in",
					getListElement(rParams, "guideTreeIn"));
        }
        if (hasClustalOmegaEntry(rParams, "guideTreeOut")) {
			appendStringValue(&argv, argc, "--guidetree-out",
					getListElement(rParams, "guideTreeOut"));
        }
        if (hasClustalOmegaEntry(rParams, "hmmIn")) {
			appendStringValue(&argv, argc, "--hmm-in",
					getListElement(rParams, "hmmIn"));
        }
        if (hasClustalOmegaEntry(rParams, "log")) {
			appendStringValue(&argv, argc, "--log",
					getListElement(rParams, "log"));
        }
        if (hasClustalOmegaEntry(rParams, "outfile")) {
			appendStringValue(&argv, argc, "--outfile",
					getListElement(rParams, "outfile"));
        }
        if (hasClustalOmegaEntry(rParams, "profile1")) {
			appendStringValue(&argv, argc, "--profile1",
					getListElement(rParams, "profile1"));
        }
        if (hasClustalOmegaEntry(rParams, "profile2")) {
			appendStringValue(&argv, argc, "--profile2",
					getListElement(rParams, "profile2"));
        }
        if (hasClustalOmegaEntry(rParams, "auto")) {
			bool autoFlag = as<bool>(params["auto"]);
			if (autoFlag) { //14
				appendString(&argv, argc, "--auto");
			}
        }
        if (hasClustalOmegaEntry(rParams, "dealign")) {
			bool dealignFlag = as<bool>(params["dealign"]);
			if (dealignFlag) { //15
				//Rprintf("DEALIGN");
				appendString(&argv, argc, "--dealign");
			}
        }
        if (hasClustalOmegaEntry(rParams, "force")) {
			bool forceFlag = as<bool>(params["force"]);
			if (forceFlag) { //16
				appendString(&argv, argc, "--force");
			}
        }
        if (hasClustalOmegaEntry(rParams, "full")) {
			bool fullFlag = as<bool>(params["full"]);
			if (fullFlag) { //17
				//Rprintf("FULL");
				appendString(&argv, argc, "--full");
			}
        }
        if (hasClustalOmegaEntry(rParams, "fullIter")) {
			bool fullIterFlag = as<bool>(params["fullIter"]);
			if (fullIterFlag) { //18
				//Rprintf("FULLITER");
				appendString(&argv, argc, "--full-iter");
			}
        }
        if (hasClustalOmegaEntry(rParams, "help")) {
			bool helpFlag = as<bool>(params["help"]);
			if (helpFlag) { //19
				//Rprintf("HELP");
				appendString(&argv, argc, "--help");
			}
        }
        if (hasClustalOmegaEntry(rParams, "isProfile")) {
			bool isProfileFlag = as<bool>(params["isProfile"]);
			if (isProfileFlag) { //20
				//Rprintf("ISPROFILE");
				appendString(&argv, argc, "--is-profile");
			}
        }
        if (hasClustalOmegaEntry(rParams, "longVersion")) {
			bool longVersionFlag = as<bool>(params["longVersion"]);
			if (longVersionFlag) { //21
				//Rprintf("longVersion");
				appendString(&argv, argc, "--long-version");
			}
        }
        if (hasClustalOmegaEntry(rParams, "percentId")) {
			bool percentIdFlag = as<bool>(params["percentId"]);
			if (percentIdFlag) { //22
				//Rprintf("percentId");
				appendString(&argv, argc, "--percent-id");
			}
        }
        if (hasClustalOmegaEntry(rParams, "residueNumber")) {
			bool residueNumberFlag = as<bool>(params["residueNumber"]);
			if (residueNumberFlag) { //23
				//Rprintf("residueNumber");
				appendString(&argv, argc, "--residuenumber");
			}
        }
        if (hasClustalOmegaEntry(rParams, "useKimura")) {
			bool useKimuraFlag = as<bool>(params["useKimura"]);
			if (useKimuraFlag) { //24
				//Rprintf("useKimura");
				appendString(&argv, argc, "--use-kimura");
			}
        }
        if (hasClustalOmegaEntry(rParams, "version")) {
			bool versionFlag = as<bool>(params["version"]);
			if (versionFlag) { //25
				//Rprintf("version");
				appendString(&argv, argc, "--version");
			}
        }

        bool verbose = as<bool>(rVerbose);
		if (verbose) {
			Rprintf("params:", argc);
			int cnt;
			for (cnt = 0; cnt < argc; cnt++) {
				Rprintf(" %s", argv[cnt]);
			}
			Rprintf("\n");
			appendString(&argv, argc, "-v");
		}

        if (!Rf_isNull(rSubstitutionMatrix)) {
            string substitutionMatrix = as<string>(rSubstitutionMatrix);
            const char* subMatrix = substitutionMatrix.c_str();
            if (strcmp("BLOSUM30", subMatrix) == 0) {
                msaInput.substitutionMatrix = 30;
                Rprintf("using BLOSUM30\n");
            } else if (strcmp("BLOSUM40", subMatrix) == 0) {
                msaInput.substitutionMatrix = 40;
                Rprintf("using BLOSUM40\n");
            } else if (strcmp("BLOSUM50", subMatrix) == 0) {
                msaInput.substitutionMatrix = 50;
                Rprintf("using BLOSUM50\n");
            } else if (strcmp("BLOSUM65", subMatrix) == 0) {
                msaInput.substitutionMatrix = 65;
                Rprintf("using BLOSUM65\n");
            } else if (strcmp("BLOSUM80", subMatrix) == 0) {
                msaInput.substitutionMatrix = 80;
                Rprintf("using BLOSUM80\n");
            } else {
                msaInput.substitutionMatrix = 0; //Gonnet
                Rprintf("using Gonnet\n");
            }
        } else {
            msaInput.substitutionMatrix = 0; //Gonnet
            Rprintf("using Gonnet\n");
        }

        ClustalOmegaOutput msaOutput;
        executeClustalOmega(argc, argv, &msaInput, &msaOutput);

        for (int i = 0; i < argc; i++) {
            free(argv[i]);
        }
        free(argv);

        /*for (int i = 0; i < n; i++) {
        	delete[](msaInput.inputSeqs);
        }*/
        if (!inputSeqIsFileFlag) {
        	delete[](*(msaInput.inputSeqs));
        }


        vector<string> result;
        for (int i = 0; i < msaOutput.msa_c; i++) {
            result.push_back(msaOutput.msa_v[i]);
            //Rprintf("free %i - %s\n", i, msaOutput.msa_v[i]);
            free(msaOutput.msa_v[i]);
        }
        free((msaOutput.msa_v));

        retList = Rcpp::List::create(Rcpp::Named("msa") =
                Rcpp::CharacterVector(result.begin(), result.end()));

    } catch(int i) {
        if (i == 0) {
            Rprintf("ClustalOmega finished successfully");
        } else {
            Rf_error("ClustalOmega finished with errors");
        }
    } catch( std::exception &ex ) {
    	Rf_error("ClustalOmega finished with errors");
        forward_exception_to_r(ex);
    } catch(...) {
        Rf_error("ClustalOmega finished by an unknown reason");
    }
    return retList;
}
