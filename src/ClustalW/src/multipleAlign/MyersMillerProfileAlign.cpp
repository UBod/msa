/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 * Changes: Mark 20-6-07, I added a call to calculateMaxlengths
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "MyersMillerProfileAlign.h"
#include "Iteration.h"
#include <math.h>

namespace clustalw
{

/**
 * 
 * @return 
 */
MyersMillerProfileAlign::MyersMillerProfileAlign()
 : _gapPos1(userParameters->getGapPos1()),
   _gapPos2(userParameters->getGapPos2())
{

}

/**
 * 
 * @param alnPtr 
 * @param distMat 
 * @param group 
 * @param aligned 
 * @return 
 */
int MyersMillerProfileAlign::profileAlign(Alignment* alnPtr, DistMatrix* distMat,
                                          vector<int>* group, int* aligned)
{
    bool negative = false;
    int i = 0, j = 0, count = 0;
    int numSeq = 0;
    int seqNum = 0;
    int len = 0, len1 = 0, len2 = 0, is = 0, minLen = 0;
    int se1 = 0, se2 = 0, sb1 = 0, sb2 = 0;
    int maxRes = 0;
    int c = 0;
    int score = 0;
    double logmin = 0.0, logdiff = 0.0;
    double pcid = 0.0;
    int matrix[NUMRES][NUMRES];
    int numSeqsProf1 = 0, numSeqsProf2 = 0;
    // Initialise the matrix. To get rid of a valgrind error.
    for(int i = 0; i < NUMRES; i++)
    {
        for(int j = 0; j < NUMRES; j++)
        {
            matrix[i][j] = 0;
        }
    }
    
    numSeq = alnPtr->getNumSeqs();
    seqArray.resize(numSeq);
    alnWeight.resize(numSeq);
    
    for (i = 0; i < numSeq; i++)
    {
        if (aligned[i + 1] == 0)
        {
            (*group)[i + 1] = 0;
        }
    }

    nseqs1 = nseqs2 = 0;
    for (i = 0; i < numSeq; i++)
    {
        if ((*group)[i + 1] == 1)
        {
            nseqs1++;
        }
        else if ((*group)[i + 1] == 2)
        {
            nseqs2++;
        }
    }

    if ((nseqs1 == 0) || (nseqs2 == 0))
    {
        return (0);
    }
    numSeqsProf1 = nseqs1;
    numSeqsProf2 = nseqs2;
    
    if (nseqs2 > nseqs1)
    {
        switchProfiles = true;
        for (i = 0; i < numSeq; i++)
        {
            if ((*group)[i + 1] == 1)
            {
                (*group)[i + 1] = 2;
            }
            else if ((*group)[i + 1] == 2)
            {
                (*group)[i + 1] = 1;
            }
        }
    }
    else
    {
        switchProfiles = false;
    }
    
    // calculate the mean of the sequence pc identities between the two groups
    
    count = 0;
    pcid = 0.0;
    negative = userParameters->getUseNegMatrix();
    
    for (i = 0; i < numSeq; i++)
    {
        if ((*group)[i + 1] == 1)
        {
            for (j = 0; j < numSeq; j++)
            {
                if ((*group)[j + 1] == 2)
                {
                    count++;
                    pcid += (*distMat)(i + 1, j + 1);
                }
            }
        }
    }

    pcid = pcid / (float)count;

    // Make the first profile.
    
    prfLength1 = 0;
    for (i = 0; i < numSeq; i++)
    {
        if ((*group)[i + 1] == 1)
        {    
            if (alnPtr->getSeqLength(i + 1) > prfLength1)
            {
                prfLength1 = alnPtr->getSeqLength(i + 1);
            }
        }
    }
    nseqs1 = 0;

    for (i = 0; i < numSeq; i++)
    {
        if ((*group)[i + 1] == 1)
        {
            len = alnPtr->getSeqLength(i + 1);
            
            // initialise with the other vector!
            seqArray[nseqs1] = vector<int>(alnPtr->getSequence(i + 1)->begin() + 1,
                                            alnPtr->getSequence(i + 1)->end());
            seqArray[nseqs1].resize(prfLength1 + extraEndElemNum);

            for (j = len; j < prfLength1; j++)
            {
                seqArray[nseqs1][j] = userParameters->getGapPos1();
            }
            alnWeight[nseqs1] = alnPtr->getSeqWeight(i);
            nseqs1++;
        }
    }


    // Make the second profile.
    
    prfLength2 = 0;
    for (i = 0; i < numSeq; i++)
    {
        if ((*group)[i + 1] == 2)
        {
            if (alnPtr->getSeqLength(i + 1) > prfLength2)
            {
                prfLength2 = alnPtr->getSeqLength(i + 1);
            }
        }
    }
    nseqs2 = 0;

    for (i = 0; i < numSeq; i++)
    {
        if ((*group)[i + 1] == 2)
        {
            len = alnPtr->getSeqLength(i + 1);

            seqArray[nseqs1 + nseqs2] = vector<int>(alnPtr->getSequence(i + 1)->begin() + 1,
                                                     alnPtr->getSequence(i + 1)->end());
            seqArray[nseqs1 + nseqs2].resize(prfLength2 + extraEndElemNum);
                        
            for (j = len; j < prfLength2; j++)
            {
                seqArray[nseqs1 + nseqs2][j] = userParameters->getGapPos1();
            }
            
            seqArray[nseqs1 + nseqs2][j] = ENDALN;
            alnWeight[nseqs1 + nseqs2] = alnPtr->getSeqWeight(i);
            //cout << "weight " << nseqs1 + nseqs2 << alnWeight[nseqs1 + nseqs2] << "\n";
            nseqs2++;
        }
    }
    
    // Change the Max alignment length in the Alignment Object! 
    alnPtr->setMaxAlnLength(prfLength1 + prfLength2 + 2);

    // calculate real length of profiles - removing gaps!
    
    len1 = 0;
    for (i = 0; i < nseqs1; i++)
    {
        is = 0;
        for (j = 0; j < utilityObject->MIN(static_cast<int>(seqArray[i].size()), prfLength1); j++)
        {
            c = seqArray[i][j];
            if ((c != _gapPos1) && (c != _gapPos2))
            {
                is++;
            }
        }
        len1 += is;
    }
    len1 = static_cast<int>(len1 / (float)nseqs1);

    len2 = 0;
    for (i = nseqs1; i < nseqs2 + nseqs1; i++)
    {
        is = 0;
        for (j = 0; j < utilityObject->MIN(static_cast<int>(seqArray[i].size()), prfLength2);
             j++)
        {
            c = seqArray[i][j];
            if ((c != _gapPos1) && (c != _gapPos2))
            {
                is++;
            }
        }
        len2 += is;
    }
    len2 = static_cast<int>(len2 / (float)nseqs2);

    PrfScaleValues scaleVals;
    scaleVals.scale = 1.0;
    scaleVals.intScale = 100.0;
    
    int matAvgScore = 0;
    minLen = utilityObject->MIN(len1, len2);
    maxRes = 0;

    // Get the substitution matrix that will be stored in 'matrix'
    maxRes = subMatrix->getProfileAlignMatrix(matrix, pcid, minLen, scaleVals, matAvgScore);

    if (maxRes == 0 || maxRes == -1)
    {
        return -1;
    }    
    if (userParameters->getDNAFlag())
    {
        gapcoef1 = gapcoef2 = static_cast<int>(100.0 * userParameters->getGapOpen() * scaleVals.scale);
        lencoef1 = lencoef2 = static_cast<int>(100.0 * userParameters->getGapExtend() * scaleVals.scale);
    }
    else
    {
        if (len1 == 0 || len2 == 0)
        {
            logmin = 1.0;
            logdiff = 1.0;
        }
        else
        {
            logmin = 1.0 / log10((double)minLen);
            if (len2 < len1)
            {
                logdiff = 1.0 + 0.5 * log10((double)((float)len2 / (float)len1));
            }
            else if (len1 < len2)
            {
                logdiff = 1.0 + 0.5 * log10((double)((float)len1 / (float)len2));
            }
            else
            {
                logdiff = 1.0;
            }
            if (logdiff < 0.9)
            {
                logdiff = 0.9;
            }
        }


        if (negative)
        {
            gapcoef1 = gapcoef2 = static_cast<int>(100.0 *(float)(userParameters->getGapOpen()));
            lencoef1 = lencoef2 = static_cast<int>(100.0 * userParameters->getGapExtend());
        }
        else
        {
            if (matAvgScore <= 0)
            {
                gapcoef1 = gapcoef2 = static_cast<int>(100.0 *(float)(userParameters->getGapOpen() + logmin));
            }
            else
            {
                gapcoef1 = gapcoef2 = static_cast<int>(scaleVals.scale * matAvgScore * 
                            (float)(userParameters->getGapOpen() / (logdiff * logmin)));
            }
            lencoef1 = lencoef2 = static_cast<int>(100.0 * userParameters->getGapExtend());
        }
    }
    // We need one profile with substitution matrix information and one without!
    // But this will change when we have the LE scoring function.
    profileWithSub = new ProfileWithSub(prfLength1, 0, nseqs1);
    profileStandard = new ProfileStandard(prfLength2, nseqs1, nseqs1 + nseqs2);

    // Calculate the profile array with Substitution matrix info 
    // calculate the Gap Coefficients.
    
    gaps.resize(alnPtr->getMaxAlnLength() + 1);
    
    bool profile1Pen; 
    if (switchProfiles == false)
    {
        profile1Pen = userParameters->getStructPenalties1() && userParameters->getUseSS1();
        profileWithSub->calcGapCoeff(&seqArray, &gaps, profile1Pen,
                                     alnPtr->getGapPenaltyMask1(), gapcoef1, lencoef1);
    }
    else
    {
        profile1Pen = userParameters->getStructPenalties2() && userParameters->getUseSS2();
        profileWithSub->calcGapCoeff(&seqArray, &gaps, profile1Pen,
                                     alnPtr->getGapPenaltyMask2(), gapcoef1, lencoef1);
    }
    // calculate the profile matrix.
    profileWithSub->calcProfileWithSub(&seqArray, &gaps, matrix, &alnWeight);
    profile1 = profileWithSub->getProfilePtr();
    
    if (userParameters->getDebug() > 4)
    {
        string aminoAcidCodes = userParameters->getAminoAcidCodes();
        for (j = 0; j <= userParameters->getMaxAA(); j++)
        {
            cout << aminoAcidCodes[j] << "    ";
        }
        cout << "\n";
        for (i = 0; i < prfLength1; i++)
        {
            for (j = 0; j <= userParameters->getMaxAA(); j++)
            {
                cout << (*profile1)[i + 1][j] << " ";
            }
            cout << (*profile1)[i + 1][_gapPos1]<< " ";
            cout << (*profile1)[i + 1][_gapPos2]<< " ";
            cout << (*profile1)[i + 1][GAPCOL] << " " << (*profile1)[i + 1][LENCOL]
                 << "\n";
        }
    }    
    // Calculate the standard profile array 
    // calculate the Gap Coefficients.
    bool profile2Pen;
    if (switchProfiles == false)
    {
        profile2Pen = userParameters->getStructPenalties2() && userParameters->getUseSS2();
        profileStandard->calcGapCoeff(&seqArray, &gaps, profile2Pen,
                                      alnPtr->getGapPenaltyMask2(), gapcoef2, lencoef2);
    }
    else
    {
        profile2Pen = userParameters->getStructPenalties1() && userParameters->getUseSS1();
        profileStandard->calcGapCoeff(&seqArray, &gaps, profile2Pen,
                                      alnPtr->getGapPenaltyMask1(), gapcoef2, lencoef2);
    }

    // calculate the profile matrix.
    
    profileStandard->calcStandardProfile(&seqArray, &alnWeight);
    profile2 = profileStandard->getProfilePtr();
    
    if (userParameters->getDebug() > 4)
    {
        string aminoAcidCodes = userParameters->getAminoAcidCodes();
        for (j = 0; j <= userParameters->getMaxAA(); j++)
        {
            cout << aminoAcidCodes[j] << "    ";
        }
        cout << "\n";
        for (i = 0; i < prfLength2; i++)
        {
            for (j = 0; j <= userParameters->getMaxAA(); j++)
            {
                cout << (*profile2)[i + 1][j] << " ";
            }
            cout << (*profile2)[i + 1][_gapPos1]<< " ";
            cout << (*profile2)[i + 1][_gapPos2]<< " ";
            cout << (*profile2)[i + 1][GAPCOL] << " " << (*profile2)[i + 1][LENCOL]
                 << "\n";
        }
    }
    alnWeight.clear();

    int _maxAlnLength = alnPtr->getMaxAlnLength();
    
    alnPath1.resize(_maxAlnLength + 1);
    alnPath2.resize(_maxAlnLength + 1);

    //align the profiles
    
    // use Myers and Miller to align two sequences 

    lastPrint = 0;
    printPtr = 1;

    sb1 = sb2 = 0;
    se1 = prfLength1;
    se2 = prfLength2;

    HH.resize(_maxAlnLength + 1);
    DD.resize(_maxAlnLength + 1);
    RR.resize(_maxAlnLength + 1);
    SS.resize(_maxAlnLength + 1);
    gS.resize(_maxAlnLength + 1);
    displ.resize(_maxAlnLength + 1);

    score = progDiff(sb1, sb2, se1 - sb1, se2 - sb2, (*profile1)[0][GAPCOL],
                    (*profile1)[prfLength1][GAPCOL]);

    // May not need if we recreate the Profile every time!
    HH.clear();
    DD.clear();
    RR.clear();
    SS.clear();
    gS.clear();

    alignmentLength = progTracepath();

    displ.clear();

    addGGaps(alnPtr, &seqArray);

    profileStandard->resetProfile();
    profileWithSub->resetProfile();
    
    prfLength1 = alignmentLength;

    alnPath1.clear();
    alnPath2.clear();
    
    if(userParameters->getDoRemoveFirstIteration() == TREE)
    {
        Iteration iterateObj;
        iterateObj.iterationOnTreeNode(numSeqsProf1, numSeqsProf2, prfLength1, prfLength2,
                                       &seqArray);           
    }
        
    //  Now we resize the SeqArray that holds the sequences in the alignment class
    //    and also update it with the new aligned sequences 
         
    SeqArray* newSequences = alnPtr->getSeqArrayForRealloc();
    seqNum = 0;
    for (j = 0; j < numSeq; j++)
    {
        if ((*group)[j + 1] == 1)
        {
            (*newSequences)[j + 1].clear();
            (*newSequences)[j + 1].resize(prfLength1 + 1);
            for (i = 0; i < prfLength1; i++)
            {
                (*newSequences)[j + 1][i + 1] = seqArray[seqNum][i];
            }
            seqNum++;
        }
    }
    for (j = 0; j < numSeq; j++)
    {
        if ((*group)[j + 1] == 2)
        {
            (*newSequences)[j + 1].clear();
            (*newSequences)[j + 1].resize(prfLength1 + 1);
            for (i = 0; i < prfLength1; i++)
            {
                (*newSequences)[j + 1][i + 1] = seqArray[seqNum][i];
            }
            seqNum++;
        }
    }
    
    alnPtr->calculateMaxLengths(); // Mark change 20-6-07
    
    gaps.clear();
    seqArray.clear();
    
    delete profileWithSub;
    delete profileStandard;

    int retScore = (score / 100);

    return retScore;
    
}

/** ****************************************************************************************
 *                          Private functions                                              *
 *******************************************************************************************/

/**
 * 
 * @param A 
 * @param B 
 * @param M 
 * @param N 
 * @param go1 
 * @param go2 
 * @return 
 */
int MyersMillerProfileAlign::progDiff(int A, int B, int M, int N, int go1, int go2)
{
    int midi, midj, type;
    int midh;

    static int t, tl, g, h;

    {
        static int i, j;
        static int hh, f, e, s;

        /* Boundary cases: M <= 1 or N == 0 */
        if (userParameters->getDebug() > 2)
        {
            cout << "A " << A << " B " << B << " M " << M << " N " << N 
                 << " midi " << M / 2 << " go1 " << go1 << " go2 " << go2<< "\n";
        }

        /* if sequence B is empty....                                            */

        if (N <= 0)
        {

            /* if sequence A is not empty....                                        */

            if (M > 0)
            {

                /* delete residues A[1] to A[M]                                          */

                progDel(M); 
            }
            return ( - gapPenalty1(A, B, M));
        }

        /* if sequence A is empty....                                            */

        if (M <= 1)
        {
            if (M <= 0)
            {

                /* insert residues B[1] to B[N]                                          */

                progAdd(N);
                return ( - gapPenalty2(A, B, N));
            }

            /* if sequence A has just one residue....                                */

            if (go1 == 0)
            {
                midh =  - gapPenalty1(A + 1, B + 1, N);
            }
            else
            {
                midh =  - gapPenalty2(A + 1, B, 1) - gapPenalty1(A + 1, B + 1,
                    N);
            }
            midj = 0;
            for (j = 1; j <= N; j++)
            {
                hh =  - gapPenalty1(A, B + 1, j - 1) + prfScore(A + 1, B + j)
                    - gapPenalty1(A + 1, B + j + 1, N - j);
                if (hh > midh)
                {
                    midh = hh;
                    midj = j;
                }
            }

            if (midj == 0)
            {
                progAdd(N);
                progDel(1);
            }
            else
            {
                if (midj > 1)
                {
                    progAdd(midj - 1);
                }
                progAlign();
                if (midj < N)
                {
                    progAdd(N - midj);
                }
            }
            return midh;
        }


        /* Divide sequence A in half: midi */

        midi = M / 2;

        /* In a forward phase, calculate all HH[j] and HH[j] */

        HH[0] = 0;
        t =  - openPenalty1(A, B + 1);
        tl =  - extPenalty1(A, B + 1);
        for (j = 1; j <= N; j++)
        {
            HH[j] = t = t + tl;
            DD[j] = t - openPenalty2(A + 1, B + j);
        }

        if (go1 == 0)
        {
            t = 0;
        }
        else
        {
            t =  - openPenalty2(A + 1, B);
        }
        tl =  - extPenalty2(A + 1, B);
        for (i = 1; i <= midi; i++)
        {
            s = HH[0];
            HH[0] = hh = t = t + tl;
            f = t - openPenalty1(A + i, B + 1);

            for (j = 1; j <= N; j++)
            {
                g = openPenalty1(A + i, B + j);
                h = extPenalty1(A + i, B + j);
                if ((hh = hh - g - h) > (f = f - h))
                {
                    f = hh;
                }
                g = openPenalty2(A + i, B + j);
                h = extPenalty2(A + i, B + j);
                if ((hh = HH[j] - g - h) > (e = DD[j] - h))
                {
                    e = hh;
                }
                hh = s + prfScore(A + i, B + j);
                if (f > hh)
                {
                    hh = f;
                }
                if (e > hh)
                {
                    hh = e;
                }

                s = HH[j];
                HH[j] = hh;
                DD[j] = e;

            }
        }

        DD[0] = HH[0];

        /* In a reverse phase, calculate all RR[j] and SS[j] */

        RR[N] = 0;
        tl = 0;
        for (j = N - 1; j >= 0; j--)
        {
            g =  - openPenalty1(A + M, B + j + 1);
            tl -= extPenalty1(A + M, B + j + 1);
            RR[j] = g + tl;
            SS[j] = RR[j] - openPenalty2(A + M, B + j);
            gS[j] = openPenalty2(A + M, B + j);
        }

        tl = 0;
        for (i = M - 1; i >= midi; i--)
        {
            s = RR[N];
            if (go2 == 0)
            {
                g = 0;
            }
            else
            {
                g =  - openPenalty2(A + i + 1, B + N);
            }
            tl -= extPenalty2(A + i + 1, B + N);
            RR[N] = hh = g + tl;
            t = openPenalty1(A + i, B + N);
            f = RR[N] - t;

            for (j = N - 1; j >= 0; j--)
            {
                g = openPenalty1(A + i, B + j + 1);
                h = extPenalty1(A + i, B + j + 1);
                if ((hh = hh - g - h) > (f = f - h - g + t))
                {
                    f = hh;
                }
                t = g;
                g = openPenalty2(A + i + 1, B + j);
                h = extPenalty2(A + i + 1, B + j);
                hh = RR[j] - g - h;
                if (i == (M - 1))
                {
                    e = SS[j] - h;
                }
                else
                {
                    e = SS[j] - h - g + openPenalty2(A + i + 2, B + j);
                    gS[j] = g;
                }
                if (hh > e)
                {
                    e = hh;
                }
                hh = s + prfScore(A + i + 1, B + j + 1);
                if (f > hh)
                {
                    hh = f;
                }
                if (e > hh)
                {
                    hh = e;
                }

                s = RR[j];
                RR[j] = hh;
                SS[j] = e;

            }
        }
        SS[N] = RR[N];
        gS[N] = openPenalty2(A + midi + 1, B + N);

        /* find midj, such that HH[j]+RR[j] or DD[j]+SS[j]+gap is the maximum */

        midh = HH[0] + RR[0];
        midj = 0;
        type = 1;
        for (j = 0; j <= N; j++)
        {
            hh = HH[j] + RR[j];
            if (hh >= midh)
            if (hh > midh || (HH[j] != DD[j] && RR[j] == SS[j]))
            {
                midh = hh;
                midj = j;
            }
        }

        for (j = N; j >= 0; j--)
        {
            hh = DD[j] + SS[j] + gS[j];
            if (hh > midh)
            {
                midh = hh;
                midj = j;
                type = 2;
            }
        }
    }

    /* Conquer recursively around midpoint                                   */


    if (type == 1)
    {
         /* Type 1 gaps  */
        if (userParameters->getDebug() > 2)
        {
            cout << "Type 1,1: midj " << midj << "\n";
        }
        progDiff(A, B, midi, midj, go1, 1);
        if (userParameters->getDebug() > 2)
        {
            cout << "Type 1,2: midj " << midj << "\n";
        }
        progDiff(A + midi, B + midj, M - midi, N - midj, 1, go2);
    }
    else
    {
        if (userParameters->getDebug() > 2)
        {
            cout << "Type 2,1: midj " << midj << "\n";
        }
        progDiff(A, B, midi - 1, midj, go1, 0);
        progDel(2);
        if (userParameters->getDebug() > 2)
        {
            cout << "Type 2,2: midj " << midj << "\n";
        }
        progDiff(A + midi + 1, B + midj, M - midi - 1, N - midj, 0, go2);
    }

    return midh; /* Return the score of the best alignment */
} 
 
/**
 * 
 * @param alnPtr 
 * @param seqArray 
 */
void MyersMillerProfileAlign::addGGaps(Alignment* alnPtr, SeqArray* seqArray)
{
    int j;
    int i, ix;
    int len;
    vector<int> ta;

    ta.resize(alignmentLength + 1);

    for (j = 0; j < nseqs1; j++)
    {
        ix = 0;
        for (i = 0; i < alignmentLength; i++)
        {
            if (alnPath1[i] == 2)
            {
                if (ix < ((int)(*seqArray)[j].size() - extraEndElemNum))
                {
                    ta[i] = (*seqArray)[j][ix];
                }
                else
                {
                    ta[i] = ENDALN;
                }
                ix++;
            }
            else if (alnPath1[i] == 1)
            {
                /*
                insertion in first alignment...
                 */
                ta[i] = _gapPos1;
            }
            else
            {
                cerr << "Error in aln_path\n";
            }
        }
        ta[i] = ENDALN;

        len = alignmentLength;

        (*seqArray)[j].resize(len + 2);
        
        for (i = 0; i < len; i++)
        {
            (*seqArray)[j][i] = ta[i];
        }
        (*seqArray)[j][len] = ENDALN;

    }

    for (j = nseqs1; j < nseqs1 + nseqs2; j++)
    {
        ix = 0;
        for (i = 0; i < alignmentLength; i++)
        {
            if (alnPath2[i] == 2)
            {
                if (ix < ((int)(*seqArray)[j].size() - extraEndElemNum))
                {
                    ta[i] = (*seqArray)[j][ix];
                }
                else
                {
                    ta[i] = ENDALN;
                }
                ix++;
            }
            else if (alnPath2[i] == 1)
            {
                /*
                insertion in second alignment...
                 */
                ta[i] = _gapPos1;
            }
            else
            {
                cerr << "Error in alnPath\n";
            }
        }
        ta[i] = ENDALN;

        len = alignmentLength;

        (*seqArray)[j].resize(len + 2);
        
        for (i = 0; i < len; i++)
        {
            (*seqArray)[j][i] = ta[i];
        }
        (*seqArray)[j][len] = ENDALN;
    }


    if (userParameters->getStructPenalties1() != NONE)
    {
        addGGapsMask(alnPtr->getGapPenaltyMask1(), alignmentLength,
                     &alnPath1, &alnPath2);
    }
    if (userParameters->getStructPenalties1() == SECST)
    {
        addGGapsMask(alnPtr->getSecStructMask1(), alignmentLength,
                     &alnPath1, &alnPath2);
    }

    if (userParameters->getStructPenalties2() != NONE)
    {
        addGGapsMask(alnPtr->getGapPenaltyMask2(), alignmentLength,
                     &alnPath2, &alnPath1);
    }
    if (userParameters->getStructPenalties2() == SECST)
    {
        addGGapsMask(alnPtr->getSecStructMask2(), alignmentLength,
                     &alnPath2, &alnPath1);
    }
}


/**
 * 
 * @param mask 
 * @param len 
 * @param path1 
 * @param path2 
 */
void MyersMillerProfileAlign::addGGapsMask(vector<char>* mask, int len, vector<int>* path1,
                                           vector<int>* path2)
{
    int i, ix;
    char *ta;

    ta = new char[len + 1];

    ix = 0;
    if (switchProfiles == false)
    {
        for (i = 0; i < len; i++)
        {
            if ((*path1)[i] == 2)
            {
                ta[i] = (*mask)[ix];
                ix++;
            }
            else if ((*path1)[i] == 1)
            {
                ta[i] = _gapPos1;
            }
        }
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            if ((*path2)[i] == 2)
            {
                ta[i] = (*mask)[ix];
                ix++;
            }
            else if ((*path2)[i] == 1)
            {
                ta[i] = _gapPos1;
            }
        }
    }
    mask->resize(len + 2);
    
    for (i = 0; i < len; i++)
    {
        (*mask)[i] = ta[i];
    }


    delete [] ta;
    ta  = NULL;
}


/**
 * 
 * @param n 
 * @param m 
 * @return 
 */
inline int MyersMillerProfileAlign::prfScore(int n, int m)
{
    int ix;
    int score;

    score = 0;
    int _maxAA = userParameters->getMaxAA(); // NOTE Change here!
    for (ix = 0; ix <= _maxAA; ix++)
    {
        score += ((*profile1)[n][ix] * (*profile2)[m][ix]);
    }
    score += ((*profile1)[n][_gapPos1] * (*profile2)[m][_gapPos1]);
    score += ((*profile1)[n][_gapPos2] * (*profile2)[m][_gapPos2]);
    return (score / 10);
}


/**
 * 
 * @return 
 */
int MyersMillerProfileAlign::progTracepath()
{
    int i, j, k, pos, toDo;
    int alignLen;
    pos = 0;

    toDo = printPtr - 1;

    for (i = 1; i <= toDo; ++i)
    {
        if (userParameters->getDebug() > 1)
        {
            cout << displ[i] << " ";
        }
        if (displ[i] == 0)
        {
            alnPath1[pos] = 2;
            alnPath2[pos] = 2;
            ++pos;
        }
        else
        {
            if ((k = displ[i]) > 0)
            {
                for (j = 0; j <= k - 1; ++j)
                {
                    alnPath2[pos + j] = 2;
                    alnPath1[pos + j] = 1;
                }
                pos += k;
            }
            else
            {
                k = (displ[i] < 0) ? displ[i] * - 1: displ[i];
                for (j = 0; j <= k - 1; ++j)
                {
                    alnPath1[pos + j] = 2;
                    alnPath2[pos + j] = 1;
                }
                pos += k;
            }
        }
    }
    if (userParameters->getDebug() > 1)
    {
        cout << "\n";
    }

    alignLen = pos;
    return alignLen;
}


/**
 * 
 * @param k 
 */
void MyersMillerProfileAlign::progDel(int k)
{
    if (lastPrint < 0)
    {
        lastPrint = displ[printPtr - 1] -= k;
    }
    else
    {
        lastPrint = displ[printPtr++] =  - (k);
    }
}


/**
 * 
 * @param k 
 */
void MyersMillerProfileAlign::progAdd(int k)
{
    if (lastPrint < 0)
    {
        displ[printPtr - 1] = k;
        displ[printPtr++] = lastPrint;
    }
    else
    {
        lastPrint = displ[printPtr++] = k;
    }
}


/**
 * 
 */
void MyersMillerProfileAlign::progAlign()
{
    displ[printPtr++] = lastPrint = 0;
}


/**
 * 
 * @param i 
 * @param j 
 * @return 
 */
inline int MyersMillerProfileAlign::openPenalty1(int i, int j) // NOTE Change here!
{
    int g;
    if (!userParameters->getEndGapPenalties() && (i == 0 || i == prfLength1))
    {
        return (0);
    }

    g = (*profile2)[j][GAPCOL] + (*profile1)[i][GAPCOL];
    return (g);
}


/**
 * 
 * @param i 
 * @param j 
 * @return 
 */
inline int MyersMillerProfileAlign::extPenalty1(int i, int j) // NOTE Change here!
{
    int h;

    if (!userParameters->getEndGapPenalties() && (i == 0 || i == prfLength1))
    {
        return (0);
    }

    h = (*profile2)[j][LENCOL];
    return (h);
}


/**
 * 
 * @param i 
 * @param j 
 * @param k 
 * @return 
 */
int MyersMillerProfileAlign::gapPenalty1(int i, int j, int k)
{
    int ix;
    int gp;
    int g, h = 0;

    if (k <= 0)
    {
        return (0);
    }
    if (!userParameters->getEndGapPenalties() && (i == 0 || i == prfLength1))
    {
        return (0);
    }

    g = (*profile2)[j][GAPCOL] + (*profile1)[i][GAPCOL];
    for (ix = 0; ix < k && ix + j < prfLength2; ix++)
    {
        h += (*profile2)[ix + j][LENCOL];
    }

    gp = g + h;
    return (gp);
}


/**
 * 
 * @param i 
 * @param j 
 * @return 
 */
inline int MyersMillerProfileAlign::openPenalty2(int i, int j) // NOTE Change here!
{
    int g;

    if (!userParameters->getEndGapPenalties() && (j == 0 || j == prfLength2))
    {
        return (0);
    }

    g = (*profile1)[i][GAPCOL] + (*profile2)[j][GAPCOL];
    return (g);
}


/**
 * 
 * @param i 
 * @param j 
 * @return 
 */
inline int MyersMillerProfileAlign::extPenalty2(int i, int j) // NOTE Change here!
{
    int h;

    if (!userParameters->getEndGapPenalties() && (j == 0 || j == prfLength2))
    {
        return (0);
    }

    h = (*profile1)[i][LENCOL];
    return (h);
}


/**
 * 
 * @param i 
 * @param j 
 * @param k 
 * @return 
 */
int MyersMillerProfileAlign::gapPenalty2(int i, int j, int k)
{
    int ix;
    int gp;
    int g, h = 0;

    if (k <= 0)
    {
        return (0);
    }
    if (!userParameters->getEndGapPenalties() && (j == 0 || j == prfLength2))
    {
        return (0);
    }

    g = (*profile1)[i][GAPCOL] + (*profile2)[j][GAPCOL];
    for (ix = 0; ix < k && ix + i < prfLength1; ix++)
    {
        h += (*profile1)[ix + i][LENCOL];
    }

    gp = g + h;
    return (gp);
}
                     
}
