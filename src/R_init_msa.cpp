#include "RClustalOmega.h"
#include "RClustalW.h"
#include "RMuscle.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
    /* RMuscle.cpp */
    {"RMuscle", (DL_FUNC) &RMuscle, 9},
    /* RClustalW.cpp */
    {"RClustalW", (DL_FUNC) &RClustalW, 9},
    /* RClustalOmega.cpp */
    {"RClustalOmega", (DL_FUNC) &RClustalOmega, 9},
    {NULL, NULL, 0}
};

extern "C" 
{
    void R_init_msa(DllInfo *info) {
		/* Register routines, allocate resources. */
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    }

    void R_unload_msa(DllInfo *info) {
		/* Release resources. */
    }
}
