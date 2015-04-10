/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "ProfileStandard.h"

namespace clustalw
{

/**
 * 
 * @param prfLen 
 * @param firstS 
 * @param lastS 
 * @return 
 */
ProfileStandard::ProfileStandard(int prfLen, int firstS, int lastS)
 : ProfileBase(prfLen, firstS, lastS)
{

}

/**
 * 
 */
void ProfileStandard::resetPrf2()
{
    profile.clear();
}


/**
 * 
 * @param seqArray 
 * @param seqWeight 
 */
void ProfileStandard::calcStandardProfile(SeqArray* seqArray, vector<int>* seqWeight)
{
    /** DONT FORGET TO CHECK THE SIZES ARE CORRECT */
    
    int sum1, sum2;
    int i, d;
    int r;
    int _maxAA = userParameters->getMaxAA();
    profile.resize(prfLength + 2, vector<int>(LENCOL + 2));
    int _gapPos1 = userParameters->getGapPos1();
    int _gapPos2 = userParameters->getGapPos2();
    
    for (r = 0; r < prfLength; r++)
    {
        /*
         * calculate sum2 = number of residues found in this column
         */
        sum2 = 0;
        for (i = firstSeq; i < lastSeq; i++)
        {
            sum2 += (*seqWeight)[i];
        }
        /*
         * only include matrix comparison scores for those residue types found in this
         * column
         */
        if (sum2 == 0)
        {
            for (d = 0; d <= _maxAA; d++)
            {
                profile[r + 1][d] = 0;
            }

            profile[r + 1][_gapPos1] = 0;
            profile[r + 1][_gapPos2] = 0;
        }
        else
        {
            for (d = 0; d <= _maxAA; d++)
            {
                sum1 = 0;
                for (i = firstSeq; i < lastSeq; i++)
                {
                    if (d == (*seqArray)[i][r])
                    {
                        sum1 += (*seqWeight)[i];
                    }
                }
                profile[r + 1][d] = (int)(10 *(float)sum1 / (float)sum2);
            }
            sum1 = 0;

            for (i = firstSeq; i < lastSeq; i++)
            {
                if (_gapPos1 == (*seqArray)[i][r])
                {
                    sum1 += (*seqWeight)[i];
                }
            }
            profile[r + 1][_gapPos1] = (int)(10 *(float)sum1 / (float)sum2);
            sum1 = 0;

            for (i = firstSeq; i < lastSeq; i++)
            {
                if (_gapPos2 == (*seqArray)[i][r])
                {
                    sum1 += (*seqWeight)[i];
                }
            }
            profile[r + 1][_gapPos2] = (int)(10 *(float)sum1 / (float)sum2);
        }
    }

}
                      

}
