/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "ProfileWithSub.h"

namespace clustalw
{

/**
 * 
 * @param prfLen 
 * @param firstS 
 * @param lastS 
 * @return 
 */
ProfileWithSub::ProfileWithSub(int prfLen, int firstS, int lastS)
 : ProfileBase(prfLen, firstS, lastS)
{
}

/**
 * 
 */
void ProfileWithSub::resetPrf1()
{
    profile.clear();
}


/**
 * 
 * @param seqArray 
 * @param gaps 
 * @param matrix[][] 
 * @param seqWeight 
 */
void ProfileWithSub::calcProfileWithSub(SeqArray* seqArray, vector<int>* gaps, 
                              int matrix[NUMRES][NUMRES], vector<int>* seqWeight)
{
    vector<vector<int> > weighting;
    int sum2, aa, seq, res;
    int _numSeq;
    int col, pos;
    int f;
    float scale;
    int _maxAA = userParameters->getMaxAA();
    int _gapPos1 = userParameters->getGapPos1();
    int _gapPos2 = userParameters->getGapPos2();
    
    weighting.resize(NUMRES + 2, vector<int>(prfLength + 2));
    
    _numSeq = lastSeq - firstSeq;

    sum2 = 0;
    for (seq = firstSeq; seq < lastSeq; seq++)
    {
        sum2 += (*seqWeight)[seq];
    }

    for (col = 0; col < prfLength; col++)
    {
        for (aa = 0; aa <= _maxAA; aa++)
        {
            weighting[aa][col] = 0;

            for (seq = firstSeq; seq < lastSeq; seq++)
                if (aa == (*seqArray)[seq][col])
                {
                    weighting[aa][col] += (*seqWeight)[seq];
                }
        }
        weighting[_gapPos1][col] = 0;

        for (seq = firstSeq; seq < lastSeq; seq++)
        {
            if (_gapPos1 == (*seqArray)[seq][col])
            {
                weighting[_gapPos1][col] += (*seqWeight)[seq];
            }
        }

        weighting[_gapPos2][col] = 0;

        for (seq = firstSeq; seq < lastSeq; seq++)
        {
            if (_gapPos2 == (*seqArray)[seq][col])
            {
                weighting[_gapPos2][col] += (*seqWeight)[seq];
            }
        }
    }

    for (pos = 0; pos < prfLength; pos++)
    {
        if ((*gaps)[pos] == _numSeq) // If all gaps
        {
            for (res = 0; res <= _maxAA; res++)
            {
                profile[pos + 1][res] = matrix[res][_gapPos1];
            }
            profile[pos + 1][_gapPos1] = matrix[_gapPos1][_gapPos1];
            profile[pos + 1][_gapPos2] = matrix[_gapPos2][_gapPos1];
        }
        else
        {
            scale = (float)(_numSeq - (*gaps)[pos]) / (float)_numSeq;
            for (res = 0; res <= _maxAA; res++)
            {
                f = 0;

                for (aa = 0; aa <= _maxAA; aa++)
                {
                    f += (weighting[aa][pos] * matrix[aa][res]);
                }

                f += (weighting[_gapPos1][pos] * matrix[_gapPos1][res]);
                f += (weighting[_gapPos2][pos] * matrix[_gapPos2][res]);
                profile[pos + 1][res] = (int)(((float)f / (float)sum2) * scale);
            }
            f = 0;

            for (aa = 0; aa <= _maxAA; aa++)
            {
                f += (weighting[aa][pos] * matrix[aa][_gapPos1]);
            }

            f += (weighting[_gapPos1][pos] * matrix[_gapPos1][_gapPos1]);
            f += (weighting[_gapPos2][pos] * matrix[_gapPos2][_gapPos1]);
            profile[pos + 1][_gapPos1] = (int)(((float)f / (float)sum2) * scale);
            f = 0;

            for (aa = 0; aa <= _maxAA; aa++)
            {
                f += (weighting[aa][pos] * matrix[aa][_gapPos2]);
            }

            f += (weighting[_gapPos1][pos] * matrix[_gapPos1][_gapPos2]);
            f += (weighting[_gapPos2][pos] * matrix[_gapPos2][_gapPos2]);
            profile[pos + 1][_gapPos2] = (int)(((float)f / (float)sum2) * scale);
        }
    }
}
                      

}
