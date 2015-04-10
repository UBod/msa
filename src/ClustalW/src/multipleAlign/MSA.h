/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef MSA_H
#define MSA_H

#include <iomanip>
#include "ProfileBase.h"
#include "../alignment/Alignment.h"
//#include "../calcAlignSteps/Tree.h"
#include "../tree/AlignmentSteps.h"
#include "ProfileAlignAlgorithm.h"

namespace clustalw
{
//using tree::AlignmentSteps;

class MSA
{
    public:
        /* Functions */
        int multiSeqAlign(Alignment* alnPtr, DistMatrix* distMat, 
            vector<int>* seqWeight, AlignmentSteps* progSteps, int iStart);
        int seqsAlignToProfile(Alignment* alnPtr, DistMatrix* distMat, vector<int>* seqWeight, int iStart, 
                              string phylipName);
        int calcPairwiseForProfileAlign(Alignment* alnPtr, DistMatrix* distMat);
        int doProfileAlign(Alignment* alnPtr, DistMatrix* distMat, vector<int>* prof1Weight,
                           vector<int>* prof2Weight);

        /* Attributes */

    private:
        /* Functions */
        
        /* Attributes */
};

}
#endif
