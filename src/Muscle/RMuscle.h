#ifndef _RMuscle_R_MUSCLE_H
#define _RMuscle_R_MUSCLE_H

#include "muscle.h"
#include "seqvect.h"

#include <Rcpp.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the RMuscle function so that
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
RcppExport SEXP RMuscle(SEXP rInputSeq,
                        SEXP rCluster,
                        SEXP rGapOpening,
                        SEXP rGapExtension,
                        SEXP rMaxiters,
                        SEXP rSubstitutionMatrix,
                        SEXP rType,
                        SEXP rVerbose,
                        SEXP rParams);

struct MuscleInput {
    SeqVect inputSeqs;
    std::vector<std::string> seqNames;
    std::vector<std::string> colNames;
    bool hasSubstitutionMatrix;
    float substitutionMatrix[32][32];
};

struct MuscleOutput { //can be used for further result objects
    std::vector<std::string> msa; //multiple sequence alignment
};

void DoMuscle(MuscleInput *msaInput, MuscleOutput *msaOutput);
void Run(MuscleInput *msaInput, MuscleOutput *msaOutput);

#endif
