/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "FullPairwiseAlign.h"
#include <math.h>

namespace clustalw
{

FullPairwiseAlign::FullPairwiseAlign()
: _maxAlnLength(0),
  intScale(0),
  mmScore(0),
  printPtr(0),
  lastPrint(0),
  _gapOpen(0),
  _gapExtend(0),
  seq1(0),
  seq2(0),
  maxScore(0),
  sb1(0),
  sb2(0),
  se1(0),
  se2(0)
{

}

void FullPairwiseAlign::pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, int iEnd, int jStart, int jEnd)
{
    int si, sj, i;
    int n, m, len1, len2;
    int maxRes;
    int _matAvgScore;
    int res;
    double _score;
    float gapOpenScale, gapExtendScale;
    
    try
    {
        
        if(distMat->getSize() != alignPtr->getNumSeqs() + 1)
        {
            cerr << "The distance matrix is not the right size!\n"
                 << "Need to terminate program.\n";
            throw 1;
        }
        if((iStart < 0) || (iEnd < iStart) || (jStart < 0) || (jEnd < jStart))
        {
            cerr << "The range for pairwise Alignment is incorrect.\n"
                 << "Need to terminate program.\n";
            throw 1;
        }
        
        _maxAlnLength = alignPtr->getMaxAlnLength();
    
        int _numSeqs = alignPtr->getNumSeqs();
        if(_numSeqs == 0)
        {
            return;
        }
    
        int num = (2 * _maxAlnLength) + 1;
        bool _DNAFlag = userParameters->getDNAFlag();
        float _pwGapOpen, _pwGapExtend;
        _pwGapOpen = userParameters->getPWGapOpen();
        _pwGapExtend = userParameters->getPWGapExtend();
        
        displ.resize(num);
        HH.resize(_maxAlnLength);
        DD.resize(_maxAlnLength);
        RR.resize(_maxAlnLength);
        SS.resize(_maxAlnLength);
        // Note these 2 lines replace the stuff above because it is all done in the SubMatrix
        PairScaleValues scaleValues;
        maxRes = subMatrix->getPairwiseMatrix(matrix, scaleValues, _matAvgScore);
        if (maxRes == 0)
        {
            cerr << "Could not get the substitution matrix\n";
            return;
        }
        
        intScale = scaleValues.intScale;
        gapOpenScale = scaleValues.gapOpenScale;
        gapExtendScale = scaleValues.gapExtendScale;
    
        int _gapPos1, _gapPos2;
        _gapPos1 = userParameters->getGapPos1();
        _gapPos2 = userParameters->getGapPos2();
        const SeqArray* _ptrToSeqArray = alignPtr->getSeqArray(); //This is faster! 
    
        for (si = utilityObject->MAX(0, iStart); si < _numSeqs && si < iEnd; si++)
        {
            n = alignPtr->getSeqLength(si + 1);
            len1 = 0;
            for (i = 1; i <= n; i++)
            {
                res = (*_ptrToSeqArray)[si + 1][i];
                if ((res != _gapPos1) && (res != _gapPos2))
                {
                    len1++;
                }
            }

            for (sj = utilityObject->MAX(si+1, jStart+1); sj < _numSeqs && sj < jEnd; sj++)
            {
                m = alignPtr->getSeqLength(sj + 1);
                if (n == 0 || m == 0)
                {
                    distMat->SetAt(si + 1, sj + 1, 1.0);
                    distMat->SetAt(sj + 1, si + 1, 1.0);
                    continue;
                }
                len2 = 0;
                for (i = 1; i <= m; i++)
                {
                    res = (*_ptrToSeqArray)[sj + 1][i];
                    if ((res != _gapPos1) && (res != _gapPos2))
                    {
                        len2++;
                    }
                }

                if (_DNAFlag)
                {
                    _gapOpen = static_cast<int>(2 * _pwGapOpen * intScale *
                                    gapOpenScale);
                    _gapExtend = static_cast<int>(_pwGapExtend * intScale * gapExtendScale);
                }
                else
                {
                    if (_matAvgScore <= 0)
                    {
                        _gapOpen = 2 * static_cast<int>((_pwGapOpen +
                               log(static_cast<double>(utilityObject->MIN(n, m)))) * intScale);
                    }
                    else
                    {
                        _gapOpen = static_cast<int>(2 * _matAvgScore * (_pwGapOpen +
                        log(static_cast<double>(utilityObject->MIN(n, m)))) * gapOpenScale);
                    }
                    _gapExtend = static_cast<int>(_pwGapExtend * intScale);
                }
                // align the sequences
            
                seq1 = si + 1;
                seq2 = sj + 1;

                _ptrToSeq1 = alignPtr->getSequence(seq1);
                _ptrToSeq2 = alignPtr->getSequence(seq2);
            
                forwardPass(_ptrToSeq1, _ptrToSeq2, n, m);
                reversePass(_ptrToSeq1, _ptrToSeq2);

                lastPrint = 0;
                printPtr = 1;

                // use Myers and Miller to align two sequences 

                maxScore = diff(sb1 - 1, sb2 - 1, se1 - sb1 + 1, se2 - sb2 + 1,
                    (int)0, (int)0);

                // calculate percentage residue identity

                mmScore = tracePath(sb1, sb2);

                if (len1 == 0 || len2 == 0)
                {
                    mmScore = 0;
                }
                else
                {
                    mmScore /= (float)utilityObject->MIN(len1, len2);
                }

                _score = ((float)100.0 - mmScore) / (float)100.0;
                distMat->SetAt(si + 1, sj + 1, _score);
                distMat->SetAt(sj + 1, si + 1, _score);
                
                if(userParameters->getDisplayInfo())
                {
                    utilityObject->info("Sequences (%d:%d) Aligned. Score:  %d",
                                        si+1, sj+1, (int)mmScore);     
                }
            }
        }

        displ.clear();
        HH.clear();
        DD.clear();
        RR.clear();
        SS.clear();
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the FullPairwiseAlign class.\n"
             << e.what() << "\n";
        throw 1;
    }
}

void FullPairwiseAlign::add(int v)
{
    if (lastPrint < 0)
    {
        displ[printPtr - 1] = v;
        displ[printPtr++] = lastPrint;
    }
    else
    {
        lastPrint = displ[printPtr++] = v;
    }
}

inline int FullPairwiseAlign::calcScore(int iat, int jat, int v1, int v2)
{
    return matrix[(*_ptrToSeq1)[v1 + iat]][(*_ptrToSeq2)[v2 + jat]];
}

float FullPairwiseAlign::tracePath(int tsb1, int tsb2)
{
    int res1, res2;
    int i1, i2;
    int i, k, pos, toDo;
    int count;
    float score;

    toDo = printPtr - 1;
    i1 = tsb1;
    i2 = tsb2;

    pos = 0;
    count = 0;
    for (i = 1; i <= toDo; ++i)
    {
        if (displ[i] == 0)
        {
            res1 = (*_ptrToSeq1)[i1];
            res2 = (*_ptrToSeq2)[i2];

            if ((res1 != userParameters->getGapPos1()) && 
                (res2 != userParameters->getGapPos2()) && (res1 == res2))
            {
                count++;
            }
            ++i1;
            ++i2;
            ++pos;
        }
        else
        {
            if ((k = displ[i]) > 0)
            {
                i2 += k;
                pos += k;
            }
            else
            {
                i1 -= k;
                pos -= k;
            }
        }
    }
    
    score = 100.0 *(float)count;
    return (score);
}

void FullPairwiseAlign::forwardPass(const vector<int>* seq1, const vector<int>* seq2, int n, int m)
{
    int i, j;
    int f, hh, p, t;

    maxScore = 0;
    se1 = se2 = 0;
    for (i = 0; i <= m; i++)
    {
        HH[i] = 0;
        DD[i] =  -_gapOpen;
    }

    for (i = 1; i <= n; i++)
    {
        hh = p = 0;
        f =  -_gapOpen;

        for (j = 1; j <= m; j++)
        {

            f -= _gapExtend;
            t = hh - _gapOpen - _gapExtend;
            if (f < t)
            {
                f = t;
            }

            DD[j] -= _gapExtend;
            t = HH[j] - _gapOpen - _gapExtend;
            if (DD[j] < t)
            {
                DD[j] = t;
            }

            hh = p + matrix[(*seq1)[i]][(*seq2)[j]];
            if (hh < f)
            {
                hh = f;
            }
            if (hh < DD[j])
            {
                hh = DD[j];
            }
            if (hh < 0)
            {
                hh = 0;
            }

            p = HH[j];
            HH[j] = hh;

            if (hh > maxScore)
            {
                maxScore = hh;
                se1 = i;
                se2 = j;
            }
        }
    }

}

void FullPairwiseAlign::reversePass(const vector<int>* seq1, const vector<int>* seq2)
{
    int i, j;
    int f, hh, p, t;
    int cost;

    cost = 0;
    sb1 = sb2 = 1;
    for (i = se2; i > 0; i--)
    {
        HH[i] =  - 1;
        DD[i] =  - 1;
    }

    for (i = se1; i > 0; i--)
    {
        hh = f =  - 1;
        if (i == se1)
        {
            p = 0;
        }
        else
        {
            p =  - 1;
        }

        for (j = se2; j > 0; j--)
        {

            f -= _gapExtend;
            t = hh - _gapOpen - _gapExtend;
            if (f < t)
            {
                f = t;
            }

            DD[j] -= _gapExtend;
            t = HH[j] - _gapOpen - _gapExtend;
            if (DD[j] < t)
            {
                DD[j] = t;
            }

            hh = p + matrix[(*seq1)[i]][(*seq2)[j]];
            if (hh < f)
            {
                hh = f;
            }
            if (hh < DD[j])
            {
                hh = DD[j];
            }

            p = HH[j];
            HH[j] = hh;

            if (hh > cost)
            {
                cost = hh;
                sb1 = i;
                sb2 = j;
                if (cost >= maxScore)
                {
                    break;
                }
            }
        }
        if (cost >= maxScore)
        {
            break;
        }
    }

}

int FullPairwiseAlign::diff(int A, int B, int M, int N, int tb, int te)
{
    int type;
    int midi, midj, i, j;
    int midh;
    static int f, hh, e, s, t;

    if (N <= 0)
    {
        if (M > 0)
        {
            del(M);
        }

        return ( -(int)tbgap(M, tb));
    }

    if (M <= 1)
    {
        if (M <= 0)
        {
            add(N);
            return ( -(int)tbgap(N, tb));
        }

        midh =  - (tb + _gapExtend) - tegap(N, te);
        hh =  - (te + _gapExtend) - tbgap(N, tb);
        if (hh > midh)
        {
            midh = hh;
        }
        midj = 0;
        for (j = 1; j <= N; j++)
        {
            hh = calcScore(1, j, A, B) - tegap(N - j, te) - tbgap(j - 1, tb);
            if (hh > midh)
            {
                midh = hh;
                midj = j;
            }
        }

        if (midj == 0)
        {
            del(1);
            add(N);
        }
        else
        {
            if (midj > 1)
            {
                add(midj - 1);
            }
            displ[printPtr++] = lastPrint = 0;
            if (midj < N)
            {
                add(N - midj);
            }
        }
        return midh;
    }

    // Divide: Find optimum midpoint (midi,midj) of cost midh

    midi = M / 2;
    HH[0] = 0;
    t =  - tb;
    for (j = 1; j <= N; j++)
    {
        HH[j] = t = t - _gapExtend;
        DD[j] = t - _gapOpen;
    }

    t =  - tb;
    for (i = 1; i <= midi; i++)
    {
        s = HH[0];
        HH[0] = hh = t = t - _gapExtend;
        f = t - _gapOpen;
        for (j = 1; j <= N; j++)
        {
            if ((hh = hh - _gapOpen - _gapExtend) > (f = f - _gapExtend))
            {
                f = hh;
            }
            if ((hh = HH[j] - _gapOpen - _gapExtend) > (e = DD[j] - _gapExtend))
            {
                e = hh;
            }
            hh = s + calcScore(i, j, A, B);
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

    RR[N] = 0;
    t =  - te;
    for (j = N - 1; j >= 0; j--)
    {
        RR[j] = t = t - _gapExtend;
        SS[j] = t - _gapOpen;
    }

    t =  - te;
    for (i = M - 1; i >= midi; i--)
    {
        s = RR[N];
        RR[N] = hh = t = t - _gapExtend;
        f = t - _gapOpen;

        for (j = N - 1; j >= 0; j--)
        {

            if ((hh = hh - _gapOpen - _gapExtend) > (f = f - _gapExtend))
            {
                f = hh;
            }
            if ((hh = RR[j] - _gapOpen - _gapExtend) > (e = SS[j] - _gapExtend))
            {
                e = hh;
            }
            hh = s + calcScore(i + 1, j + 1, A, B);
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
        hh = DD[j] + SS[j] + _gapOpen;
        if (hh > midh)
        {
            midh = hh;
            midj = j;
            type = 2;
        }
    }

    // Conquer recursively around midpoint 


    if (type == 1)
    {
        // Type 1 gaps
        diff(A, B, midi, midj, tb, _gapOpen);
        diff(A + midi, B + midj, M - midi, N - midj, _gapOpen, te);
    }
    else
    {
        diff(A, B, midi - 1, midj, tb, 0);
        del(2);
        diff(A + midi + 1, B + midj, M - midi - 1, N - midj, 0, te);
    }

    return midh; // Return the score of the best alignment
}

void FullPairwiseAlign::del(int k)
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

int FullPairwiseAlign::gap(int k)
{
    if(k <= 0)
    {
        return 0;
    }
    else
    {
        return _gapOpen + _gapExtend * k;
    }
}

int FullPairwiseAlign::tbgap(int k, int tb)
{
    if(k <= 0)
    {
        return 0;
    }
    else
    {
        return tb + _gapExtend * k;
    }
}

int FullPairwiseAlign::tegap(int k, int te)
{
    if(k <= 0)
    {
        return 0;
    }
    else
    {
        return te + _gapExtend * k;
    }
}

}
