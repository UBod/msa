/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef PROFILESTANDARD_H
#define PROFILESTANDARD_H

#include "../alignment/Alignment.h"
#include "ProfileBase.h"

namespace clustalw
{

class ProfileStandard : public ProfileBase
{
    public:
        /* Functions */
        ProfileStandard(int prfLen, int firstS, int lastS);
        void resetPrf2();
        void calcStandardProfile(SeqArray* alignment, vector<int>* seqWeight); 

        /* Attributes */

    private:
        /* Functions */

        /* Attributes */
};

}
#endif
