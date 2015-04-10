/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "LowScoreSegProfile.h"

namespace clustalw
{

LowScoreSegProfile::LowScoreSegProfile(int prfLen, int firstS, int lastS)
    : prfLength(prfLen),
      firstSeq(firstS),
      lastSeq(lastS)
{
    profile.resize(prfLength + 2, vector<int>(LENCOL + 2));
}

void LowScoreSegProfile::calcLowScoreSegProfile(const SeqArray* seqArray, 
                                int matrix[NUMRES][NUMRES], vector<int>* seqWeight)
{
    vector<vector<int> > weighting; 
    int d, i, res; 
    int r, pos;
    int f;
    int _gapPos1 = userParameters->getGapPos1();
    int _gapPos2 = userParameters->getGapPos2();
    int _maxAA = userParameters->getMaxAA();

    weighting.resize(NUMRES + 2, vector<int>(prfLength + 2));
    
    for (r = 0; r < prfLength; r++)
    {
        for (d = 0; d <= _maxAA; d++)
        {
            weighting[d][r] = 0;
            for (i = firstSeq; i < lastSeq; i++)
            {
                if (r + 1 < (int)(*seqArray)[i + 1].size() - 1)
                {
                    if (d == (*seqArray)[i + 1][r + 1])
                    { 
                        weighting[d][r] += (*seqWeight)[i];
                    }
                }
            }
        }
        
        weighting[_gapPos1][r] = 0;
        
        for (i = firstSeq; i < lastSeq; i++)
        {
            if (r + 1 < (int)(*seqArray)[i + 1].size() - 1)
            {
                if (_gapPos1 == (*seqArray)[i + 1][r + 1])
                { 
                    weighting[_gapPos1][r] += (*seqWeight)[i];
                }
            }
        }
        
        weighting[_gapPos2][r] = 0;
        
        for (i = firstSeq; i < lastSeq; i++)
        {
            if (r + 1 < (int)(*seqArray)[i + 1].size() - 1)
            {
                if (_gapPos2 == (*seqArray)[i + 1][r + 1])
                { 
                    weighting[_gapPos2][r] += (*seqWeight)[i];
                }
            }
        }
    }

    for (pos = 0; pos < prfLength; pos++)
    {
        for (res = 0; res <= _maxAA; res++)
        {
            f = 0;
            
            for (d = 0; d <= _maxAA; d++)
            {
                f += (weighting[d][pos] * matrix[d][res]);
            }
            
            f += (weighting[_gapPos1][pos] * matrix[_gapPos1][res]);
            f += (weighting[_gapPos2][pos] * matrix[_gapPos2][res]);
            profile[pos + 1][res] = f;
        }
        f = 0;
        
        for (d = 0; d <= _maxAA; d++)
        {
            f += (weighting[d][pos] * matrix[d][_gapPos1]);
        }   
        
        f += (weighting[_gapPos1][pos] * matrix[_gapPos1][_gapPos1]);
        f += (weighting[_gapPos2][pos] * matrix[_gapPos2][_gapPos1]);
        profile[pos + 1][_gapPos1] = f;
        f = 0;
           
        for (d = 0; d <= _maxAA; d++)
        {
            f += (weighting[d][pos] * matrix[d][_gapPos2]);
        }   
        f += (weighting[_gapPos1][pos] * matrix[_gapPos1][_gapPos2]);
        f += (weighting[_gapPos2][pos] * matrix[_gapPos2][_gapPos2]);
        profile[pos + 1][_gapPos2] = f;
    }
}                                

}
