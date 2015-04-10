/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * NOTE: This profile is not used in the multiple alignment part. It is used for the
 * the clustal Qt part. It is used in the calculation of low scoring segments.
 */
#ifndef LOWSCORESEGPROFILE_H
#define LOWSCORESEGPROFILE_H

#include "../alignment/Alignment.h"

namespace clustalw
{

class LowScoreSegProfile
{
    public:
        /* Functions */
        LowScoreSegProfile(int prfLen, int firstS, int lastS);
        void calcLowScoreSegProfile(const SeqArray* seqArray, 
                                int matrix[NUMRES][NUMRES], vector<int>* seqWeight);
        const SeqArray* getProfilePtr(){return &profile;};                        
        /* Attributes */

    protected:
        /* Functions */

        /* Attributes */
        SeqArray profile; 
        int prfLength;
        int firstSeq, lastSeq;
    private:
        /* Functions */

        /* Attributes */

};

}
#endif
