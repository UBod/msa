/* -*- mode: c; tab-width: 4; c-basic-offset: 4;  indent-tabs-mode: nil -*- */

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
 *  RCS $Id: clustal-omega.c 284 2013-06-12 10:10:11Z fabian $
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>

#include "clustal-omega.h"
#include "hhalign/general.h"
#include "clustal/hhalign_wrapper.h"

/* The following comment block contains the frontpage/mainpage of the doxygen
 *  documentation. Please add some more info. FIXME add more
 */

/**
 *
 * @mainpage Clustal-Omega Documentation
 *
 * @section intro_sec Introduction
 *
 * For more information see http://www.clustal.org/
 *
 * @section api_section API
 *
 * @subsection example_prog_subsection An Example Program
 *
 * To use libclustalo you will have to include the clustal-omega.h header and
 * link against libclustalo. For linking against libclustalo you will have to
 * use a C++ compiler, no matter if your program was written in C or C++. See
 * below (\ref pkgconfig_subsubsec)) on how to figure out compiler flags with
 * pkg-config.
 *
 * Assuming Clustal Omega was installed in system-wide default directory (e.g.
 * /usr), first compile (don't link yet) your source (for example code see
 * section \ref example_src_subsubsec) and then link against libclustalo:
 *
 * @code
 * $ gcc -c -ansi -Wall clustalo-api-test.c
 * $ g++  -ansi -Wall -o clustalo-api-test clustalo-api-test.o  -lclustalo
 * @endcode
 *
 * Voila! Now you have your own alignment program based on Clustal Omega which
 * can be run with
 *
 * @code
 * $ ./clustalo-api-test <your-sequence-input>
 * @endcode
 *  
 * It's best to use the same compiler that you used for compiling libclustal.
 * If libclustal was compiled with OpenMP support, you will have to use OpenMP
 * flags for you program as well.
 *
 *
 * @subsubsection pkgconfig_subsubsec Using pkg-config / Figuring out compiler flags
 *
 * Clustal Omega comes with support for <a
 * href="http://pkg-config.freedesktop.org">pkg-config</a>, which means you
 * can run
 *
 * @code
 * $ pkg-config --cflags --libs clustalo
 * @endcode
 *
 * to figure out cflags and library flags needed to compile and link against
 * libclustalo. This is especially handy if Clustal Omega was installed to a
 * non-standard directory.
 *  
 * You might have to change PKG_CONFIG_PATH. For example, if you used the prefix $HOME/local/ for
 * installation then you will first need to set PKG_CONFIG_PATH:
 *
 * @code
 * $ export PKG_CONFIG_PATH=$HOME/local/lib/pkgconfig
 * $ pkg-config --cflags --libs clustalo
 * @endcode
 *  
 *  
 * To compile your source use as above but this time using proper flags:
 *  
 * @code
 * $ export PKG_CONFIG_PATH=$HOME/local/lib/pkgconfig
 * $ gcc -c -ansi -Wall $(pkg-config --cflags clustalo) clustalo-api-test.c  
 * $ g++  -ansi -Wall -o clustalo-api-test $(pkg-config --libs clustalo) clustalo-api-test.o 
 * @endcode
 *
 *
 * @subsubsection example_src_subsubsec Example Source Code
 *
 * @include "clustalo-api-test.c"
 *
 *
 */


/* the following are temporary flags while the code is still under construction;
   had problems internalising hhmake, so as temporary crutch 
   write alignment to file and get external hmmer/hhmake via system call 
   to read alignment and convert into HMM
   All this will go, once hhmake is properly internalised */
#define INDIRECT_HMM 0 /* temp flag: (1) write aln to file, use system(hmmer/hhmake), (0) internal hhmake */
#define USEHMMER 1 /* temp flag: use system(hmmer) to build HMM */
#define USEHHMAKE (!USEHMMER) /* temp flag: use system(hhmake) to build HMM */


/* shuffle order of input sequences */
#define SHUFFLE_INPUT_SEQ_ORDER 0

/* sort input sequences by length */
#define SORT_INPUT_SEQS 0


int iNumberOfThreads;

/* broken, unused and lonely */
static const int ITERATION_SCORE_IMPROVEMENT_THRESHOLD = 0.01;


/**
 * @brief Print Long version information to pre-allocated char. 
 *
 * @note short version
 * information is equivalent to PACKAGE_VERSION
 *
 * @param[out] pcStr
 * char pointer to write to preallocated to hold iSize chars.
 * @param[in] iSize
 * size of pcStr
 */
void 
PrintLongVersion(char *pcStr, int iSize)
{
    snprintf(pcStr, iSize, "version %s; code-name '%s'; build date %s",
             PACKAGE_VERSION, PACKAGE_CODENAME, __DATE__);
}
/* end of PrintLongVersion() */



/**
 * @brief free aln opts members
 *
 */
void
FreeAlnOpts(opts_t *prAlnOpts) {
    if (NULL != prAlnOpts->pcGuidetreeInfile) {
        CKFREE(prAlnOpts->pcGuidetreeInfile);
    }
    if (NULL != prAlnOpts->pcGuidetreeOutfile) {
        CKFREE(prAlnOpts->pcGuidetreeOutfile);
    }
    if (NULL  != prAlnOpts->pcDistmatOutfile) {
        CKFREE(prAlnOpts->pcDistmatOutfile);
    }
    if (NULL != prAlnOpts->pcDistmatInfile) {
        CKFREE(prAlnOpts->pcDistmatInfile);
    }
}
/* end of FreeAlnOpts() */



/**
 * @brief Sets members of given user opts struct to default values
 *
 * @param[out] prOpts
 * User opt struct to initialise
 *
 */
void
SetDefaultAlnOpts(opts_t *prOpts) {
    prOpts->bAutoOptions = FALSE;
    
    prOpts->pcDistmatInfile = NULL;
    prOpts->pcDistmatOutfile = NULL;

    prOpts->iClustersizes = 100; /* FS, r274 -> */
    prOpts->pcClustfile = NULL; /* FS, r274 -> */
    prOpts->iClusteringType = CLUSTERING_UPGMA;
    prOpts->iPairDistType = PAIRDIST_KTUPLE;
    prOpts->bUseMbed = TRUE; /* FS, r250 -> */
    prOpts->bUseMbedForIteration = TRUE; /* FS, r250 -> */
    prOpts->pcGuidetreeOutfile = NULL;
    prOpts->pcGuidetreeInfile = NULL;
    prOpts->bPercID = FALSE;

    prOpts->ppcHMMInput = NULL;
    prOpts->iHMMInputFiles = 0;
		
    prOpts->iNumIterations = 0;
    prOpts->bIterationsAuto = FALSE;
    prOpts->iMaxGuidetreeIterations = INT_MAX;
    prOpts->iMaxHMMIterations = INT_MAX;

    SetDefaultHhalignPara(& prOpts->rHhalignPara);
 }
/* end of SetDefaultAlnOpts() */



/**
 * @brief Check logic of parsed user options. Will exit (call Log(&rLog, LOG_FATAL, ))
 * on Fatal logic error
 *
 * @param[in] prOpts
 * Already parsed user options
 * 
 */    
void
AlnOptsLogicCheck(opts_t *prOpts)
{
    /* guide-tree & distmat
     *
     */
    if (prOpts->pcDistmatInfile && prOpts->pcGuidetreeInfile) {
        Log(&rLog, LOG_FATAL, "Read distances *and* guide-tree from file doesn't make sense.");
    }

    if (prOpts->pcDistmatOutfile && prOpts->pcGuidetreeInfile) {
        Log(&rLog, LOG_FATAL, "Won't be able to save distances to file, because I got a guide-tree as input.");
    }

    /* combination of options that don't make sense when not iterating
     */
    if (prOpts->iNumIterations==0 && prOpts->bIterationsAuto != TRUE)  {
        
        if (prOpts->pcGuidetreeInfile && prOpts->pcGuidetreeOutfile) {
            Log(&rLog, LOG_FATAL, "Got a guide-tree as input and output which doesn't make sense when not iterating.");            
        }
        /*
          if (prOpts->pcGuidetreeInfile && prOpts->bUseMbed > 0) {
          Log(&rLog, LOG_FATAL, "Got a guide-tree as input and was requested to cluster with mBed, which doesn't make sense when not iterating.");            
          }
        */
        /*
          AW: bUseMbedForIteration default since at least R252
          if (prOpts->bUseMbedForIteration > 0) {
            Log(&rLog, LOG_FATAL, "No iteration requested, but mbed for iteration was set. Paranoia exit.");
          }        
        */
    }
    
    if (prOpts->rHhalignPara.iMacRamMB < 512) {
        Log(&rLog, LOG_WARN, "Memory for MAC Algorithm quite low, Viterbi Algorithm may be triggered.");
        if (prOpts->rHhalignPara.iMacRamMB < 1) {
            Log(&rLog, LOG_WARN, "Viterbi Algorithm always turned on, increase MAC-RAM to turn on MAC.");
        }
    }

    return; 
}
/* end of AlnOptsLogicCheck() */


/**
 * @brief FIXME doc
 */
void
PrintAlnOpts(FILE *prFile, opts_t *prOpts)
{
    int iAux;


    /* keep in same order as struct */
    fprintf(prFile, "option: auto-options = %d\n", prOpts->bAutoOptions);
    fprintf(prFile, "option: distmat-infile = %s\n", 
            NULL != prOpts->pcDistmatInfile? prOpts->pcDistmatInfile: "(null)");
	fprintf(prFile, "option: distmat-outfile = %s\n",
            NULL != prOpts->pcDistmatOutfile? prOpts->pcDistmatOutfile: "(null)");
	fprintf(prFile, "option: clustering-type = %d\n", prOpts->iClusteringType);
	fprintf(prFile, "option: pair-dist-type = %d\n", prOpts->iPairDistType);
	fprintf(prFile, "option: use-mbed = %d\n", prOpts->bUseMbed);
	fprintf(prFile, "option: use-mbed-for-iteration = %d\n", prOpts->bUseMbedForIteration);
	fprintf(prFile, "option: guidetree-outfile = %s\n", 
            NULL != prOpts->pcGuidetreeOutfile? prOpts->pcGuidetreeOutfile: "(null)");
	fprintf(prFile, "option: guidetree-infile = %s\n",
            NULL != prOpts->pcGuidetreeInfile? prOpts->pcGuidetreeInfile: "(null)");
    for (iAux=0; iAux<prOpts->iHMMInputFiles; iAux++) {
        fprintf(prFile, "option: hmm-input no %d = %s\n", iAux, prOpts->ppcHMMInput[iAux]);
    }
	fprintf(prFile, "option: hmm-input-files = %d\n", prOpts->iHMMInputFiles);
	fprintf(prFile, "option: num-iterations = %d\n", prOpts->iNumIterations);
	fprintf(prFile, "option: iterations-auto = %d\n", prOpts->bIterationsAuto);
	fprintf(prFile, "option: max-hmm-iterations = %d\n", prOpts->iMaxHMMIterations);
	fprintf(prFile, "option: max-guidetree-iterations = %d\n", prOpts->iMaxGuidetreeIterations);
    fprintf(prFile, "option: iMacRamMB = %d\n", prOpts->rHhalignPara.iMacRamMB);
    fprintf(prFile, "option: percent-id = %d\n", prOpts->bPercID);
    fprintf(prFile, "option: use-kimura = %d\n", prOpts->bUseKimura);

}
/* end of PrintAlnOpts() */



/**
 * @brief Returns major version of HMMER. Whichever hmmbuild version
 * is found first in your PATH will be used
 *
 * @return -1 on error, major hmmer version otherwise
 *
 */
int
HmmerVersion()
{
    char zcHmmerTestCall[] = "hmmbuild -h";
    FILE *fp = NULL;
    int iMajorVersion = 0;
    char zcLine[16384];

    if (NULL == (fp = popen(zcHmmerTestCall, "r"))) {
        Log(&rLog, LOG_ERROR, "Couldn't exec %s", zcHmmerTestCall);
        return -1;
    }
    while (fgets(zcLine, sizeof(zcLine), fp)) {
        char *pcLocate;
        if ((pcLocate = strstr(zcLine, "HMMER "))) {
            iMajorVersion = atoi(&pcLocate[6]);
            break;
        }
    }
    pclose(fp);

    return iMajorVersion;
}
/* end of HmmerVersion() */



/**
 * @brief Create a HHM file from aligned sequences
 *
 * @warning Should be eliminated in the future 
 * as building routine should not create intermediate files
 *
 * @param[in] prMSeq
 * Aligned mseq_t
 * @param[in] pcHMMOut
 * HMM output file name
 *
 * @return Non-zero on error
 *
 */
int
AlnToHHMFile(mseq_t *prMSeq, char *pcHMMOut)
{
    char *tmp_aln = NULL;
    int retcode = OK;

    int msa_c = 0;
    char **msa_v = NULL;

    assert(NULL!=prMSeq);
    assert(NULL!=pcHMMOut);

    if (FALSE == prMSeq->aligned) {
        Log(&rLog, LOG_ERROR, "Sequences need to be aligned to create an HMM");
        return FAILURE;
    }

    /* Convert alignment to a2m, and call hhmake
     *
     * can't be static templates, or mktemp fails (at least on os x
     * (with a bus error))
     * 
     * gcc says we should use mkstemp to avoid race conditions, 
     * but that returns a file descriptor, which is of no use to 
     * us  
     */
    /* NOTE: the following won't work on windows: missing /tmp/ */
    tmp_aln = CkStrdup("/tmp/clustalo_tmpaln_XXXXXX");
    if (NULL == mktemp(tmp_aln)) {
        Log(&rLog, LOG_ERROR, "Could not create temporary alignment filename");
        retcode = FAILURE;
        goto cleanup_and_return;
    }
#define LINE_WRAP 60

    if (WriteAlignment(prMSeq, tmp_aln, MSAFILE_A2M, LINE_WRAP, FALSE, &msa_c, &msa_v)) {
        Log(&rLog, LOG_ERROR, "Could not save alignment to %s", tmp_aln);
        retcode = FAILURE;
        goto cleanup_and_return;
    }

    if (HHMake_Wrapper(tmp_aln, pcHMMOut)){
        Log(&rLog, LOG_ERROR, "Could not convert alignment %s into HHM", tmp_aln);
        retcode = FAILURE;
        goto cleanup_and_return;
    }


 cleanup_and_return:

    if (NULL != tmp_aln) {
        if (FileExists(tmp_aln)) {
            if (remove(tmp_aln)) {
                Log(&rLog, LOG_WARN, "Removing %s failed. Continuing anyway", tmp_aln);
            }
        }
        CKFREE(tmp_aln);
    }

    return retcode; 

} /* end of AlnToHHMFile() */



/**
 * @brief Create a HMM file from aligned sequences
 *
 * @warning Should be replaced in the future by some internal HMM
 * building routine that does not call external programs
 *
 * @param[in] prMSeq
 * Aligned mseq_t
 * @param[in] pcHMMOut
 * HMM output file name
 *
 * @return Non-zero on error
 *

 */
int
AlnToHMMFile(mseq_t *prMSeq, const char *pcHMMOut)
{
    char *tmp_aln = NULL;
    char *tmp_hmm = NULL; /* only needed for hmmer3 to hmmer2 conversion */
    char cmdbuf[16384];
    int iHmmerVersion = 0;
    int retcode = OK;
    
    int msa_c = 0;
    char **msa_v = NULL;

    assert(NULL!=prMSeq);
    assert(NULL!=pcHMMOut);

    if (FALSE == prMSeq->aligned) {
        Log(&rLog, LOG_ERROR, "Sequences need to be aligned to create an HMM");
        return FAILURE;
    }

    iHmmerVersion = HmmerVersion();
    if (2 != iHmmerVersion && 3 != iHmmerVersion) {
        Log(&rLog, LOG_ERROR, "Could not find suitable HMMER binaries");
        return FAILURE;
    }

    /* Convert alignment to stockholm, call hmmbuild and then
     * either hmmconvert (hmmer3) or hmmcalibrate (hmmer2)
     *
     * can't be static templates, or mktemp fails (at least on os x
     * (with a bus error))
     *
     * gcc says we should use mkstemp to avoid race conditions,
     * but that returns a file descriptor, which is of no use to
     * us
     */
    /* NOTE: the following won't work on windows: missing /tmp/ */
    tmp_aln = CkStrdup("/tmp/clustalo_tmpaln_XXXXXX");
    if (NULL == mktemp(tmp_aln)) {
        Log(&rLog, LOG_ERROR, "Could not create temporary alignment filename");
        retcode = FAILURE;
        goto cleanup_and_return;
    }
#define LINE_WRAP 60
    if (WriteAlignment(prMSeq, tmp_aln, MSAFILE_STOCKHOLM, LINE_WRAP, FALSE, &msa_c, &msa_v)) {
        Log(&rLog, LOG_ERROR, "Could not save alignment to %s", tmp_aln);
        retcode = FAILURE;
        goto cleanup_and_return;
    }

    if (2 == iHmmerVersion) {
        sprintf(cmdbuf, "hmmbuild %s %s >/dev/null && hmmcalibrate %s >/dev/null",
                pcHMMOut, tmp_aln, pcHMMOut);
        if (system(cmdbuf)) {
            Log(&rLog, LOG_ERROR, "Command '%s' failed", cmdbuf);
            retcode = FAILURE;
            goto cleanup_and_return;
        }
    } else if (3 == iHmmerVersion) {
        /* NOTE: the following won't work on windows: missing /tmp/ */
        tmp_hmm = CkStrdup("/tmp/clustalo_tmphmm2_XXXXXX");
        if (NULL == mktemp(tmp_hmm)) {
            Log(&rLog, LOG_ERROR, "Could not create temporary hmm filename");
            retcode = FAILURE;
            goto cleanup_and_return;
        }
        sprintf(cmdbuf, "hmmbuild %s %s >/dev/null && hmmconvert -2 %s > %s",
                tmp_hmm, tmp_aln, tmp_hmm, pcHMMOut);
        if (system(cmdbuf)) {
            Log(&rLog, LOG_ERROR, "Command '%s' failed", cmdbuf);
            retcode = FAILURE;
            goto cleanup_and_return;
        }
    } else {
        CKFREE(tmp_aln);
        Log(&rLog, LOG_FATAL, "Internal error: Unknown Hmmer version %d", iHmmerVersion);
    }


 cleanup_and_return:

    if (NULL != tmp_aln) {
        if (FileExists(tmp_aln)) {
            if (remove(tmp_aln)) {
                Log(&rLog, LOG_WARN, "Removing %s failed. Continuing anyway", tmp_aln);
            }
        }
        CKFREE(tmp_aln);
    }
    if (NULL != tmp_hmm) {
        if (FileExists(tmp_hmm)) {
            if (remove(tmp_hmm)) {
                Log(&rLog, LOG_WARN, "Removing %s failed. Continuing anyway", tmp_hmm);
            }
        }
        CKFREE(tmp_hmm);
    }

    return retcode;
}
/* end of AlnToHMMFile() */



/**
 * @brief Convert a multiple sequence structure into a HMM
 *
 * @param[out] prHMM
 * Pointer to preallocted HMM which will be set here
 * @param[in] prMSeq
 * Pointer to an alignment
 *
 * @return 0 on error, non-0 otherwise
 *
 * @see AlnToHMMFile()
 * 
 */    
int
AlnToHMM(hmm_light *prHMM, mseq_t *prMSeq)
{
    char *pcHMM; /* temp hmm file */

    Log(&rLog, LOG_INFO,
        "Using HMMER version %d to calculate a new HMM.",
         HmmerVersion());
    /* FIXME replace all this with internal HMM computation (HHmake) */

    /**
     * @warning the following probably won't work on windows: missing
     * /tmp/. Should be ok on Cygwin though
     */
    pcHMM = CkStrdup("/tmp/clustalo-hmm-iter_XXXXXX");
    if (NULL == mktemp(pcHMM)) {
        Log(&rLog, LOG_ERROR, "Could not create temporary hmm filename");
        CKFREE(pcHMM);
        return FAILURE;
    }
    
    /* Create a HMM representing the current alignment
     */
#if USEHMMER
    if (AlnToHMMFile(prMSeq, pcHMM)) {
        Log(&rLog, LOG_ERROR, "AlnToHMMFile() on %s failed.", pcHMM);
        CKFREE(pcHMM);
        return FAILURE;
    }
#elif USEHHMAKE
    if (AlnToHHMFile(prMSeq, pcHMM)) {
        Log(&rLog, LOG_ERROR, "AlnToHHMFile() on %s failed.", pcHMM);
        CKFREE(pcHMM);
        return FAILURE;
    }
    /* Log(&rLog, LOG_FATAL, "Method to create HHM (HMM using hhmake) not installed yet"); */
#else
    Log(&rLog, LOG_FATAL, "Unknown method to create temporary HMM");
#endif
    
    /* Read HMM information
     */
    if (OK != readHMMWrapper(prHMM, pcHMM)){
        Log(&rLog, LOG_ERROR, "Processing of HMM file %s failed", pcHMM);
        CKFREE(pcHMM);
        return FAILURE;
    }

    if (remove(pcHMM)) {
        Log(&rLog, LOG_WARN, "Removing %s failed. Continuing anyway", pcHMM);
    }
    CKFREE(pcHMM);
        
    return OK; 
}
/* end of AlnToHMM() */



/**
 * @brief FIXME
 *
 */
void
InitClustalOmega(int iNumThreadsRequested)
{

#ifdef WIN32
    if (iNumThreadsRequested>1) {
        Log(&rLog, LOG_FATAL, "Cannot change number of threads to %d. %s was build without OpenMP support.",
              iNumThreadsRequested, PACKAGE_NAME);
    }
    iNumberOfThreads = 1; /* need to set this, even if build without support */
#elif HAVE_OPENMP
    iNumberOfThreads = iNumThreadsRequested;
    omp_set_num_threads(iNumberOfThreads);
#else
	if (iNumThreadsRequested>1) {
        Log(&rLog, LOG_FATAL, "Cannot change number of threads to %d. %s was build without OpenMP support.",
              iNumThreadsRequested, PACKAGE_NAME);
    }
    iNumberOfThreads = 1; /* need to set this, even if build without support */
#endif

    Log(&rLog, LOG_INFO, "Using %d threads",
         iNumberOfThreads);

}
/* end of InitClustalOmega() */



/**
 * @brief Defines an alignment order, which adds sequences
 * sequentially, i.e. one at a time starting with seq 1 & 2
 *
 * @param[out] piOrderLR_p
 * order in which nodes/profiles are to be merged/aligned
 * @param[in] iNumSeq
 * Number of sequences
 *
 * @see TraverseTree()
 *
 */
void
SequentialAlignmentOrder(int **piOrderLR_p, int iNumSeq)
{
    unsigned int uNodes = iNumSeq*2-1;
    unsigned int uNodeCounter = 0;
    unsigned int uSeqCounter = 0;

    Log(&rLog, LOG_FATAL, "FIXME: Untested...");
    
    (*piOrderLR_p) = (int *)CKCALLOC(DIFF_NODE * uNodes, sizeof(int));
    /* loop over merge nodes, which have per definition even indices
     * and set up children which have odd indices
     */
    uSeqCounter = 0;
    for (uNodeCounter=iNumSeq; uNodeCounter<uNodes; uNodeCounter+=1) {
         unsigned int uLeftChildNodeIndex = uNodeCounter-1;
         unsigned int uRightChildNodeIndex = uNodeCounter-iNumSeq+1;
         unsigned int uParentNodeIndex = uNodeCounter+1;

         /* merge node setup */
        (*piOrderLR_p)[DIFF_NODE*uNodeCounter+LEFT_NODE] = uLeftChildNodeIndex;
        (*piOrderLR_p)[DIFF_NODE*uNodeCounter+RGHT_NODE] = uRightChildNodeIndex;
        (*piOrderLR_p)[DIFF_NODE*uNodeCounter+PRNT_NODE] = uParentNodeIndex;
        /* only setup left child if at first merge node, all other left childs
         * should be merge nodes that are already set up. also correct
         * left node number here.
         */
        if (uNodeCounter==iNumSeq) {
            (*piOrderLR_p)[DIFF_NODE*uNodeCounter+LEFT_NODE] = 0;

            (*piOrderLR_p)[0+LEFT_NODE] = 0;
            (*piOrderLR_p)[0+RGHT_NODE] = 0;
            (*piOrderLR_p)[0+PRNT_NODE] = uNodeCounter;
            uSeqCounter++;

            Log(&rLog, LOG_FORCED_DEBUG, "Set up first leaf with node counter %d: left=%d right=%d parent=%d",
                      0,
                      (*piOrderLR_p)[DIFF_NODE*uLeftChildNodeIndex+LEFT_NODE],
                      (*piOrderLR_p)[DIFF_NODE*uLeftChildNodeIndex+RGHT_NODE],
                      (*piOrderLR_p)[DIFF_NODE*uLeftChildNodeIndex+PRNT_NODE]);
        }
        Log(&rLog, LOG_FORCED_DEBUG, "Set up merge node with node counter %d: left=%d right=%d parent=%d",
                  uNodeCounter, (*piOrderLR_p)[DIFF_NODE*uNodeCounter+LEFT_NODE],
                  (*piOrderLR_p)[DIFF_NODE*uNodeCounter+RGHT_NODE],
                  (*piOrderLR_p)[DIFF_NODE*uNodeCounter+PRNT_NODE]);
        
        /* right child */
        (*piOrderLR_p)[DIFF_NODE*uRightChildNodeIndex+LEFT_NODE] = uSeqCounter;
        (*piOrderLR_p)[DIFF_NODE*uRightChildNodeIndex+RGHT_NODE] = uSeqCounter;
        (*piOrderLR_p)[DIFF_NODE*uRightChildNodeIndex+PRNT_NODE] = uNodeCounter;
        uSeqCounter++;

        Log(&rLog, LOG_FORCED_DEBUG, "Set up leaf with node counter %d: left=%d right=%d parent=%d",
                  uRightChildNodeIndex, (*piOrderLR_p)[DIFF_NODE*uRightChildNodeIndex+LEFT_NODE],
                  (*piOrderLR_p)[DIFF_NODE*uRightChildNodeIndex+RGHT_NODE],
                  (*piOrderLR_p)[DIFF_NODE*uRightChildNodeIndex+PRNT_NODE]);
    }
}
/* end of SequentialAlignmentOrder() */



/**
 * @brief Defines the alignment order by calculating a guide tree. In
 * a first-step pairwise distances will be calculated (or read from a
 * file). In a second step those distances will be clustered and a
 * guide-tree created. Steps 1 and 2 will be skipped if a guide-tree
 * file was given, in which case the guide-tree will be just read from
 * the file.
 *
 * @param[out] piOrderLR_p
 * order in which nodes/profiles are to be merged/aligned
 * @param[out] pdSeqWeights_p
 * Sequence weights
 * @param[out] pdSeqWeights_p
 * Sequence weights
 * @param[in] prMSeq
 * The sequences from which the alignment order is to be calculated
 * @param[in] iPairDistType
 * Method of pairwise distance comparison
 * @param[in] pcDistmatInfile
 * If not NULL distances will be read from this file instead of being
 * calculated
 * @param[in] pcDistmatOutfile
 * If not NULL computed pairwise distances will be written to this file
 * @param[in] iClusteringType
 * Clustering method to be used to cluster the pairwise distances
 * @param[in] pcGuidetreeInfile
 * If not NULL guidetree will be read from this file. Skips pairwise
 * distance and guidetree computation
 * @param[in] pcGuidetreeOutfile
 * If not NULL computed guidetree will be written to this file
 * @param[in] bUseMbed
 * If TRUE, fast mBed guidetree computation will be employed
 *
 * @return Non-zero on error
 *
 */
int
AlignmentOrder(int **piOrderLR_p, double **pdSeqWeights_p, mseq_t *prMSeq,
               int iPairDistType, char *pcDistmatInfile, char *pcDistmatOutfile,
               int iClusteringType, int iClustersizes, 
               char *pcGuidetreeInfile, char *pcGuidetreeOutfile, char *pcClusterFile, 
               bool bUseMbed, bool bPercID)
{
    /* pairwise distance matrix (tmat in 1.83) */
    symmatrix_t *distmat = NULL;
    /* guide tree */
    tree_t *prTree = NULL;
    int i = 0;

    
    /* Shortcut for only two sequences: Do not compute k-tuple
     * distances. Use the same logic as in TraverseTree() to setup
     * piOrderLR_p. Changes there will have to be reflected here as
     * well. */
    if (2==prMSeq->nseqs) {
        Log(&rLog, LOG_VERBOSE,
            "Have only two sequences: No need to compute pairwise score and compute a tree.");

        if (NULL != pcDistmatOutfile){
            Log(&rLog, LOG_WARN, "Have only two sequences: Will not calculate/print distance matrix.");
        }

        (*piOrderLR_p) = (int*) CKMALLOC(DIFF_NODE * 3 * sizeof(int));
        (*piOrderLR_p)[DIFF_NODE*0+LEFT_NODE] = 0;
        (*piOrderLR_p)[DIFF_NODE*0+RGHT_NODE] = 0;
        (*piOrderLR_p)[DIFF_NODE*0+PRNT_NODE] = 0;

        (*piOrderLR_p)[DIFF_NODE*1+LEFT_NODE] = 1;
        (*piOrderLR_p)[DIFF_NODE*1+RGHT_NODE] = 1;
        (*piOrderLR_p)[DIFF_NODE*1+PRNT_NODE] = 1;

        /* root */
        (*piOrderLR_p)[DIFF_NODE*2+LEFT_NODE] = 0;
        (*piOrderLR_p)[DIFF_NODE*2+RGHT_NODE] = 1;
        (*piOrderLR_p)[DIFF_NODE*2+PRNT_NODE] = 2;

        /* Same logic as CalcClustalWeights(). Changes there will
           have to be reflected here as well. */
#if USE_WEIGHTS
         (*pdWeights_p) = (double *) CKMALLOC(uNodeCount * sizeof(double));
         (*pdWeights_p)[0] = 0.5;
         (*pdWeights_p)[1] = 0.5;
#endif
        
        return OK;
    }

        
    /* compute distance & guide tree, alternatively read distances or
     * guide tree from file
     *
     */
    if (NULL != pcGuidetreeInfile) {
        Log(&rLog, LOG_INFO, "Reading guide-tree from %s", pcGuidetreeInfile);
        if (GuideTreeFromFile(&prTree, prMSeq, pcGuidetreeInfile)) {
            Log(&rLog, LOG_ERROR, "Reading of guide tree %s failed.", pcGuidetreeInfile);
            return FAILURE;
        }

    } else {

        if (bUseMbed) {
            if (Mbed(&prTree, prMSeq, iPairDistType, pcGuidetreeOutfile, iClustersizes, pcClusterFile)) {
                Log(&rLog, LOG_ERROR, "mbed execution failed.");
                return FAILURE;
            }
            Log(&rLog, LOG_INFO, "Guide-tree computation (mBed) done.");
            if (NULL != pcDistmatOutfile) {
                Log(&rLog, LOG_INFO,
                    "Ignoring request to write distance matrix (am in mBed mode)");
            }
        } else {

            if (PairDistances(&distmat, prMSeq, iPairDistType, bPercID, 
                              0, prMSeq->nseqs, 0, prMSeq->nseqs,
                              pcDistmatInfile, pcDistmatOutfile)) {
                Log(&rLog, LOG_ERROR, "Couldn't compute pair distances");
                return FAILURE;
            }

            /* clustering of distances to get guide tree
             */
            if (CLUSTERING_UPGMA == iClusteringType) {
                char **labels;
                labels = (char**) CKMALLOC(prMSeq->nseqs * sizeof(char*));
                for (i=0; i<prMSeq->nseqs; i++) {
                    labels[i] = prMSeq->sqinfo[i].name;
                }

                GuideTreeUpgma(&prTree, labels, distmat, pcGuidetreeOutfile);
                Log(&rLog, LOG_INFO, "Guide-tree computation done.");

                CKFREE(labels);
            } else {
                Log(&rLog, LOG_FATAL, "INTERNAL ERROR %s",
                      "clustering method should have been checked before");
            }
        } /* did not use mBed */
    } /* had to calculate tree (did not read from file) */


#if USE_WEIGHTS
    /* derive sequence weights from tree
     *
     */
    Log(&rLog, LOG_INFO, "Calculating sequence weights");
    CalcClustalWeights(pdSeqWeights_p, prTree);
    for (i = 0; i < GetLeafCount(prTree); i++) {
        Log(&rLog, LOG_VERBOSE,
            "Weight for seq no %d: %s = %f",
             i, prMSeq->sqinfo[i].name, (*pdSeqWeights_p)[i]);
    }
#else
    Log(&rLog, LOG_DEBUG, "Not using weights");
#endif

    
    /* define traversing order of tree
     *
     */
    TraverseTree(piOrderLR_p, prTree, prMSeq);
    if (rLog.iLogLevelEnabled <= LOG_DEBUG) {
        /* FIXME: debug only, FS */
        uint uNodeIndex;
        FILE *fp = LogGetFP(&rLog, LOG_INFO);
        Log(&rLog, LOG_DEBUG, "left/right order after tree traversal");
        for (uNodeIndex = 0; uNodeIndex < GetNodeCount(prTree); uNodeIndex++) {
            fprintf(fp, "%3d:\t%2d/%2d -> %d\n", i,
                   (*piOrderLR_p)[DIFF_NODE*uNodeIndex+LEFT_NODE],
                   (*piOrderLR_p)[DIFF_NODE*uNodeIndex+RGHT_NODE],
                   (*piOrderLR_p)[DIFF_NODE*uNodeIndex+PRNT_NODE]);
        }
    }

    FreeMuscleTree(prTree);
    FreeSymMatrix(&distmat);

#if 0
    Log(&rLog, LOG_FATAL, "DEBUG EXIT before leaving %s", __FUNCTION__);
#endif
    return OK;
}
/* end of AlignmentOrder() */



/**
 * @brief Set some options automatically based on number of sequences. Might
 * overwrite some user-set options.
 *
 * @param[out] prOpts
 * Pointer to alignment options structure
 * @param[in] iNumSeq
 * Number of sequences to align
 */
void
SetAutoOptions(opts_t *prOpts, int iNumSeq) {

    Log(&rLog, LOG_INFO,
        "Setting options automatically based on input sequence characteristics (might overwrite some of your options).");

    /* AW: new version of mbed is always good (uses subclusters) */
    if (FALSE == prOpts->bUseMbed) {
        Log(&rLog, LOG_INFO, "Auto settings: Enabling mBed.");
        prOpts->bUseMbed = TRUE;
    }
    
    if (iNumSeq >= 1000) {
        if (0 != prOpts->iNumIterations) {
            Log(&rLog, LOG_INFO, "Auto settings: Disabling iterations.");
            prOpts->iNumIterations = 0;
        }
        
    } else if (iNumSeq < 1000) {
        if (1 != prOpts->iNumIterations) {
            Log(&rLog, LOG_INFO, "Auto settings: Setting iteration to 1.");
            prOpts->iNumIterations = 1;
        }        
    }
}
/* end of */



/**
 * @brief The main alignment function which wraps everything else.
 *
 * @param[out] prMSeq
 * *the* multiple sequences structure
 * @param[in] prMSeqProfile
 * optional profile to align against
 * @param[in] prOpts
 * alignment options to use
 *
 * @return 0 on success, -1 on failure
 *
 */
int
Align(mseq_t *prMSeq, 
      mseq_t *prMSeqProfile,
      opts_t *prOpts) {
   
    /* HMM
     */
    /* structs with pseudocounts etc; one for each HMM infile, i.e.
     * index range: 0..iHMMInputFiles */
    hmm_light *prHMMs = NULL;

    /* MSA order in which nodes/profiles are to be merged/aligned
       (order of nodes in guide tree (left/right)*/
    int *piOrderLR = NULL;

    /* weights per sequence */
    double *pdSeqWeights = NULL;

    /* Iteration
     */
    int iIterationCounter = 0;
    double dAlnScore;
    /* last dAlnScore for iteration */
    double dLastAlnScore = -666.666;

    int i, j; /* aux */

    assert(NULL != prMSeq);
    if (NULL != prMSeqProfile) {
        assert(TRUE == prMSeqProfile->aligned);
    }


    /* automatic setting of options
     *
     */
    if (prOpts->bAutoOptions) {
        SetAutoOptions(prOpts, prMSeq->nseqs);
    }


#if SHUFFLE_INPUT_SEQ_ORDER
    /* 
     * shuffle input:  only useful for testing/debugging 
     */
    Log(&rLog, LOG_WARN, "Shuffling input sequences! (Will also change output order)");
    ShuffleMSeq(prMSeq);
#endif
        

#if SORT_INPUT_SEQS
    /* 
     * sort input:
     *
     * would ensure we *always* (unless we get into the mbed k-means stage)
     * get the same answer. usually you don't, because most pairwise alignment
     * scores are in theory not symmetric, therefore sequence ordering might
     * have an effect on the guide-tree. Sorting by length should get rid of
     * this (and takes no time even for 100k seqs). Benchmark results on
     * Balibase show almost no difference after sorting.
     */
    Log(&rLog, LOG_WARN, "Sorting input seq by length! This will also change the output order");
    SortMSeqByLength(prMSeq, 'd');

#endif


    /* Read backgrounds HMMs and store in prHMMs
     *
     */
    if (0 < prOpts->iHMMInputFiles) {
        int iHMMInfileIndex;
        
        /**
         * @warning old structure used to be initialised like this:
         * hmm_light rHMM = {0};
         */
        prHMMs = (hmm_light *) CKMALLOC(prOpts->iHMMInputFiles * sizeof(hmm_light));
        
        for (iHMMInfileIndex=0; iHMMInfileIndex<prOpts->iHMMInputFiles; iHMMInfileIndex++) {
            char *pcHMMInput = prOpts->ppcHMMInput[iHMMInfileIndex];
            if (OK != readHMMWrapper(&prHMMs[iHMMInfileIndex], pcHMMInput)){
                Log(&rLog, LOG_ERROR, "Processing of HMM file %s failed", pcHMMInput);
                return -1;
            }
            
#if 0
            Log(&rLog, LOG_FORCED_DEBUG, "HMM length is %d", prHMMs[iHMMInfileIndex].L);
            Log(&rLog, LOG_FORCED_DEBUG, "n-display  is %d", prHMMs[iHMMInfileIndex].n_display);
            for (i = 0; NULL != prHMMs[prOpts->iHMMInputFiles].seq[i]; i++){
                printf("seq[%d]: %s\n", i, prHMMs[iHMMInfileIndex].seq[i]);
            }
            Log(&rLog, LOG_FORCED_DEBUG, "Neff_HMM   is %f", prHMMs[iHMMInfileIndex].Neff_HMM);
#endif
            if (rLog.iLogLevelEnabled <= LOG_DEBUG){
                Log(&rLog, LOG_DEBUG, "print frequencies");
                for (i = 0; i < prHMMs[iHMMInfileIndex].L; i++){
#define PRINT_TAIL 5
                    if ( (PRINT_TAIL+1 == i) && (prHMMs[iHMMInfileIndex].L-PRINT_TAIL != i) ){
                        printf("....\n");
                    }
                    if ( (i > PRINT_TAIL) && (i < prHMMs[iHMMInfileIndex].L-PRINT_TAIL) ){
                        continue;
                    }
                    printf("%3d:", i);
                    for (j = 0; j < 20; j++){
                        printf("\t%1.3f", prHMMs[iHMMInfileIndex].f[i][j]);
                    }
                    printf("\n");
                }
            } /* debug print block */
            
            CKFREE(prOpts->ppcHMMInput[iHMMInfileIndex]);
        } /* for each background HMM file */
        CKFREE(prOpts->ppcHMMInput);
    } /* there were background HMM files */
    

    
    /* If the input ("non-profile") sequences are aligned, then turn
     * the alignment into a HMM and add to the list of background HMMs
     *
     */
    if (TRUE == prMSeq->aligned) {
        /* FIXME: gcc warns about missing initialiser here (-Wall -Wextra -pedantic) */
        hmm_light rHMMLocal = {0};
        
        Log(&rLog, LOG_INFO,
            "Input sequences are aligned. Will turn alignment into HMM and add it to the user provided background HMMs.");
        /* certain gap parameters ('~' MSF) cause problems, 
           sanitise them; FS, r258 -> r259 */
        SanitiseUnknown(prMSeq);
        if (OK !=
#if INDIRECT_HMM 
            AlnToHMM(&rHMMLocal, prMSeq)
#else
            AlnToHMM2(&rHMMLocal, prMSeq->seq, prMSeq->nseqs)
#endif
            ) {
            Log(&rLog, LOG_ERROR, "Couldn't convert aligned input sequences to HMM. Will try to continue");
        } else {
            prHMMs = (hmm_light *) CKREALLOC(prHMMs, ((prOpts->iHMMInputFiles+1) * sizeof(hmm_light)));
            memcpy(&(prHMMs[prOpts->iHMMInputFiles]), &rHMMLocal, sizeof(hmm_light));
            prOpts->iHMMInputFiles++;
        }
    }

    
    /* If we have a profile turn it into a HMM and add to
     * the list of background HMMs.
     *
     */
    if (NULL != prMSeqProfile) {
        /* FIXME: gcc warns about missing initialiser here (-Wall -Wextra -pedantic) */
        hmm_light rHMMLocal = {0};
        Log(&rLog, LOG_INFO,
            "Turning profile1 into HMM and will use it during progressive alignment.");
        if (OK !=
#if INDIRECT_HMM 
            AlnToHMM(&rHMMLocal, prMSeqProfile)
#else 
            AlnToHMM2(&rHMMLocal, prMSeqProfile->seq, prMSeqProfile->nseqs)
#endif
            ) {
            Log(&rLog, LOG_ERROR, "Couldn't convert profile1 to HMM. Will try to continue");
        } else {
            prHMMs = (hmm_light *) CKREALLOC(prHMMs, ((prOpts->iHMMInputFiles+1) * sizeof(hmm_light)));
            memcpy(&(prHMMs[prOpts->iHMMInputFiles]), &rHMMLocal, sizeof(hmm_light));
            prOpts->iHMMInputFiles++;
        }
    }


    /* Now do a first alignment of the input sequences (prMSeq) adding
     * all collected background HMMs
     *
     */
    /* Determine progressive alignment order
     */
    if (TRUE == prMSeq->aligned) {
        if ( (SEQTYPE_PROTEIN == prMSeq->seqtype) && (TRUE == prOpts->bUseKimura) ){
            Log(&rLog, LOG_INFO, "%s %s",
                "Input sequences are aligned.",
                "Will use Kimura distances of aligned sequences.");
            prOpts->iPairDistType = PAIRDIST_SQUIDID_KIMURA;
        }
        else {
            prOpts->iPairDistType = PAIRDIST_SQUIDID;
        }
    }

#if 0
    Log(&rLog, LOG_WARN, "Using a sequential alignment order.");
    SequentialAlignmentOrder(&piOrderLR, prMSeq->nseqs);
#else
    if (OK != AlignmentOrder(&piOrderLR, &pdSeqWeights, prMSeq,
                             prOpts->iPairDistType,
                             prOpts->pcDistmatInfile, prOpts->pcDistmatOutfile,
                             prOpts->iClusteringType, prOpts->iClustersizes, 
                             prOpts->pcGuidetreeInfile, prOpts->pcGuidetreeOutfile, prOpts->pcClustfile, 
                             prOpts->bUseMbed, prOpts->bPercID)) {
        Log(&rLog, LOG_ERROR, "AlignmentOrder() failed. Cannot continue");
        return -1;
    }
#endif

    /* if max-hmm-iter is set < 0 then do not perform alignment 
     * there is a problem/feature(?) that the un-aligned sequences are output 
     */
    if (prOpts->iMaxHMMIterations < 0){
        Log(&rLog, LOG_VERBOSE,
            "iMaxHMMIterations < 0 (%d), will not perform alignment", prOpts->iMaxHMMIterations);
        return 0;
    }


    /* Progressive alignment of input sequences. Order defined by
     * branching of guide tree (piOrderLR). Use optional
     * background HMM information (prHMMs[0..prOpts->iHMMInputFiles-1])
     *
     */
    dAlnScore = HHalignWrapper(prMSeq, piOrderLR, pdSeqWeights,
                               2*prMSeq->nseqs -1/* nodes */,
                               prHMMs, prOpts->iHMMInputFiles, -1, prOpts->rHhalignPara);
    dLastAlnScore = dAlnScore;
    Log(&rLog, LOG_VERBOSE,
        "Alignment score for first alignment = %f", dAlnScore);        




    /* ------------------------------------------------------------
     *
     * prMSeq is aligned now. Now start iterations if requested and save the
     * alignment at the very end.
     *
     * @note We discard the background HMM information at this point,
     * because it was already used. Could consider to make this choice
     * optional.  FIXME
     *
     * ------------------------------------------------------------ */

    
    /* iteration after first alignment was computed (if not profile-profile
     * alignment)
     *
     */
    for (iIterationCounter=0;
         (iIterationCounter < prOpts->iNumIterations || prOpts->bIterationsAuto);
         iIterationCounter++) {

        hmm_light rHMMLocal = {0};
        /* FIXME Keep copy of old alignment in case new one sucks? */


        if (iIterationCounter >= prOpts->iMaxHMMIterations
            &&
            iIterationCounter >= prOpts->iMaxGuidetreeIterations) {
            Log(&rLog, LOG_VERBOSE, "Reached maximum number of HMM and guide-tree iterations");
            break;
        }
        
        if (! prOpts->bIterationsAuto) {
            Log(&rLog, LOG_INFO, "Iteration step %d out of %d", 
                 iIterationCounter+1, prOpts->iNumIterations);
        } else {
            Log(&rLog, LOG_INFO, "Iteration step %d out of <auto>", 
                 iIterationCounter+1);
        }
#if 0
        if (rLog.iLogLevelEnabled <= LOG_VERBOSE) {
            char zcIntermediate[1000] = {0};
            char *pcFormat = "fasta";
            sprintf(zcIntermediate, "clustalo-aln-iter~%d~", iIterationCounter);
#define LINE_WRAP 60
            if (WriteAlignment(prMSeq, zcIntermediate, MSAFILE_A2M, LINE_WRAP, &msa_c, &msa_v)) {
                Log(&rLog, LOG_ERROR, "Could not save alignment to %s", zcIntermediate);
                return -1;
            }
        }
#endif


        /* new guide-tree
         *
         */
        if (iIterationCounter < prOpts->iMaxGuidetreeIterations) {
            /* determine progressive alignment order
             *
             * few things are different now when calling AlignmentOrder:
             * - we have to ignore prOpts->pcDistmatInfile and pcGuidetreeInfile
             *   as they were used before
             * - the corresponding outfiles are still valid though
             */
            /* Free stuff that has already been allocated by or further
             * downstream of AlignmentOrder()
             */
            if (NULL != piOrderLR)
                CKFREE(piOrderLR);
            if (NULL != pdSeqWeights)
                CKFREE(pdSeqWeights);
            Log(&rLog, LOG_INFO, "Computing new guide tree (iteration step %d)");
            if (AlignmentOrder(&piOrderLR, &pdSeqWeights, prMSeq,
                               ((SEQTYPE_PROTEIN == prMSeq->seqtype) && (TRUE == prOpts->bUseKimura)) ? PAIRDIST_SQUIDID_KIMURA : PAIRDIST_SQUIDID, 
                               NULL, prOpts->pcDistmatOutfile,
                               prOpts->iClusteringType, prOpts->iClustersizes,
                               NULL, prOpts->pcGuidetreeOutfile, prOpts->pcClustfile, 
                               prOpts->bUseMbedForIteration, prOpts->bPercID)) {
                Log(&rLog, LOG_ERROR, "AlignmentOrder() failed. Cannot continue");
                return -1;
            }
        } else {
            Log(&rLog, LOG_INFO, "Skipping guide-tree iteration at iteration step %d (reached maximum)", 
                iIterationCounter);
        }


        /* new local hmm iteration
         *
         */
        /* non-residue/gap characters will crash AlnToHMM2(), 
           therefore sanitise unknown characters, FS, r259 -> r260 */
        SanitiseUnknown(prMSeq);
        if (iIterationCounter < prOpts->iMaxHMMIterations) {
            Log(&rLog, LOG_INFO, "Computing HMM from alignment");
            
            if (OK !=
#if INDIRECT_HMM 
                AlnToHMM(&rHMMLocal, prMSeq)
#else
                AlnToHMM2(&rHMMLocal, prMSeq->seq, prMSeq->nseqs)
#endif
                ) {
                Log(&rLog, LOG_ERROR, "Couldn't convert alignment to HMM. Will stop iterating now...");
                break;
            }
        } else {
            Log(&rLog, LOG_INFO, "Skipping HMM iteration at iteration step %d (reached maximum)", 
                 iIterationCounter);
        }

        
        /* align the sequences (again)
         */
        dAlnScore = HHalignWrapper(prMSeq, piOrderLR, pdSeqWeights,
                                   2*prMSeq->nseqs -1/* nodes */, &rHMMLocal, 1, -1, 
                                   prOpts->rHhalignPara);
        Log(&rLog, LOG_VERBOSE,
            "Alignment score for alignmnent in hmm-iteration no %d = %f (last score = %f)",
             iIterationCounter+1, dAlnScore, dLastAlnScore);
        
        
        FreeHMMstruct(&rHMMLocal);

#if 0
        /* FIXME: need a better score for automatic iteration */
        if (prOpts->bIterationsAuto) {
            /* automatic iteration: break if score improvement was not
             * big enough
             */
            double dScoreImprovement = (dAlnScore-dLastAlnScore)/dLastAlnScore;
            if (dScoreImprovement < ITERATION_SCORE_IMPROVEMENT_THRESHOLD) {
                Log(&rLog, LOG_INFO,
                    "Stopping after %d guide-tree iterations. No further alignment score improvement achieved.",
                     iIterationCounter+1);
                /* use previous alignment */
                FreeMSeq(&prMSeq);
                Log(&rLog, LOG_FORCED_DEBUG, "FIXME: %s", "CopyMSeq breaks things in this context");
                CopyMSeq(&prMSeq, prMSeqCopy);
                /* FIXME: prOpts->pcDistmatOutfile and pcGuidetreeOutfile
                 * might have been updated, but then discarded here?
                 */
                break;
            } else {
                Log(&rLog, LOG_INFO,
                    "Got a %d%% better score in iteration step %d",
                     (int)dScoreImprovement*100, iIterationCounter+1);
                FreeMSeq(&prMSeqCopy);
            }
        }
        dLastAlnScore = dAlnScore;
#endif

    }
    /* end of iterations */



    /* Last step: if a profile was also provided then align now-aligned mseq
     * with this profile
     *
     * Don't use the backgrounds HMMs anymore and don't iterate.
     * (which was done before).
     *
     */
    if (NULL != prMSeqProfile) {
        if (AlignProfiles(prMSeq, prMSeqProfile, prOpts->rHhalignPara)) {
            Log(&rLog, LOG_ERROR, "An error occured during the profile/profile alignment");
            return -1;
        }
    }

    
    if (NULL != piOrderLR) {
        CKFREE(piOrderLR);
    }
    if (NULL != pdSeqWeights) {
        CKFREE(pdSeqWeights);
    }
    if (0 < prOpts->iHMMInputFiles) {
        for (i=0; i<prOpts->iHMMInputFiles; i++) {
            FreeHMMstruct(&prHMMs[i]);
        }
        CKFREE(prHMMs);
    }

    return 0;
}
/* end of Align() */




/**
 * @brief Align two profiles, ie two sets of prealigned sequences. Already
 * aligned columns won't be changed.
 *
 * @param[out] prMSeqProfile1
 * First profile/aligned set of sequences. Merged alignment will be found in
 * here.
 * @param[in] prMSeqProfile2
 * First profile/aligned set of sequences
 * @param[in] rHhalignPara
 * FIXME
 *
 * @return 0 on success, -1 on failure
 *
 */
int
AlignProfiles(mseq_t *prMSeqProfile1, 
              mseq_t *prMSeqProfile2, hhalign_para rHhalignPara) {
   
    double dAlnScore;

    /* number of seqs in first half of joined profile */
    int iProfProfSeparator = prMSeqProfile1->nseqs;

    assert(TRUE == prMSeqProfile1->aligned);
    assert(TRUE == prMSeqProfile2->aligned);

    Log(&rLog, LOG_INFO, "Performing profile/profile alignment");

    /* Combine the available mseqs into prMSeq
     * which will be aligned afterwards.
     */
    JoinMSeqs(&prMSeqProfile1, prMSeqProfile2);

        
    /* set alignment flag explicitly to FALSE */
    prMSeqProfile1->aligned = FALSE;
                    
    dAlnScore = HHalignWrapper(prMSeqProfile1,
                               NULL, /* no order */
                               NULL, /* no weights */
                               3, /* nodes: root+2profiles */
                               NULL, 0 /* no bg-hmms */,
                               iProfProfSeparator, rHhalignPara);
        
    Log(&rLog, LOG_VERBOSE, "Alignment score is = %f", dAlnScore);

    return 0;
}
/* end of AlignProfiles() */


