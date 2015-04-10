/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef PROFILEWITHSUB_H
#define PROFILEWITHSUB_H

#include "../alignment/Alignment.h"
#include "ProfileBase.h"

namespace clustalw
{

class ProfileWithSub : public ProfileBase
{
    public:
        /* Functions */
        ProfileWithSub(int prfLen, int firstS, int lastS);
        void resetPrf1();
        void calcProfileWithSub(SeqArray* seqArray, vector<int>* gaps, 
                                int matrix[NUMRES][NUMRES], vector<int>* seqWeight); 

        /* Attributes */

    private:
        /* Functions */

        /* Attributes */
};

}

#endif
