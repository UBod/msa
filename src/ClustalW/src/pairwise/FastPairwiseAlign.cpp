/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <math.h>
#include "FastPairwiseAlign.h"

namespace clustalw
{

FastPairwiseAlign::FastPairwiseAlign()
{
    _maxAlnLength = 0;
}

void FastPairwiseAlign::pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, 
                                      int iEnd, int jStart, int jEnd)
{
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
            cout << "The range for pairwise Alignment is incorrect.\n"
                 << "Need to terminate program.\n";
            throw 1;
        }
    
        int i, j, dsr;
        double calcScore;
        bool _DNAFlag = userParameters->getDNAFlag();
        _maxAlnLength = alignPtr->getMaxAlnLength();
        int num = (2 * _maxAlnLength) + 1;
        accum.ResizeRect(5, num);
    
        displ.resize(num);
        slopes.resize(num);
        diagIndex.resize(num);

        zza.resize(_maxAlnLength + 1);
        zzb.resize(_maxAlnLength + 1);
        zzc.resize(_maxAlnLength + 1);
        zzd.resize(_maxAlnLength + 1);
    
        if (_DNAFlag)
        {
            userParameters->setDNAParams();
        }
        else
        {
            userParameters->setProtParams();
        }

        cout << "\n\n";
    
        for (i = iStart + 1; i <= iEnd; ++i)
        {
            const vector<int>* _seqIPtr = alignPtr->getSequence(i);
            int _seqILength = alignPtr->getSeqLength(i);
            if (_DNAFlag)
            {
                makeNPtrs(zza, zzc, _seqIPtr, _seqILength);
            }
            else
            {
                makePPtrs(zza, zzc, _seqIPtr, _seqILength);
            }
            double _score;
	    for (j = jStart + 2 > i+1 ? jStart + 2 : i+1; j <= jEnd; ++j)
            {
                const vector<int>* _seqJPtr = alignPtr->getSequence(j);
                int _seqJLength = alignPtr->getSeqLength(j);
                if (_DNAFlag)
                {
                    makeNPtrs(zzb, zzd, _seqJPtr, _seqJLength);
                }
                else
                {
                    makePPtrs(zzb, zzd, _seqJPtr, _seqJLength);
                }
                pairAlign(_seqIPtr, _seqILength, _seqJLength);
                if (!maxSoFar)
                {
                    calcScore = 0.0;
                }
                else
                {
                    calcScore = (double)accum[0][maxSoFar];
                    if (userParameters->getPercent())
                    {
                        dsr = (_seqILength < _seqJLength) ? _seqILength
                            : _seqJLength;
                        calcScore = (calcScore / (double)dsr) *100.0;
                    }
                }
                _score = (100.0 - calcScore) / 100.0;
                distMat->SetAt(i, j, _score);
                //distMat->SetAt(j, i, _score); /* distMat symmetric, FS, 2009-04-06 */

            
                if(userParameters->getDisplayInfo())
                {
                    if (calcScore > 0.1)
                    {
                        utilityObject->info("Sequences (%d:%d) Aligned. Score: %lg",
                                            i, j, calcScore);     
                    }
                    else
                    {
                        utilityObject->info("Sequences (%d:%d) Not Aligned", i, j);
                    }
                }
            }
        }
        accum.clearArray();    
        displ.clear();
        slopes.clear();
        diagIndex.clear();

        zza.clear();
        zzb.clear();
        zzc.clear();
        zzd.clear();
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the FastPairwiseAlign class.\n"
             << e.what() << "\n";
        throw 1;
    }    
}


/*
 * Note: There is a problem with the treatment of DNA/RNA. 
 * During file reading all residues are encoded as AminoAcids, 
 * even before it has been established if they are AA or not. 
 * This is bad and will have to be changed (later). 
 * 'A' is assigned code 0, C is assigned 2, G = 6, T=18, U=19. 
 * However, the fast alignment routines require that 
 * A=0, C=1, G=2, T=U=3. In the best case the results of the 
 * fast alignment (of DNA) will simply be meaningless, the 
 * worst case is a core dump. 
 * As a quick fix I implemented the following mask, that 
 * (for DNA/RNA) translates 0->0, 2->1, 6->2, 18->3, 19->3. 
 * This is awfull (it is not OO compliant) but it was quick, 
 * uses (much) less memory than a second DNA array, and is 
 * (much) faster than calling a translation function. 
 * Ideally this will be removed, but this requires changes 
 * to (i) the sequence encoding during file-reading AND 
 * (ii) all the DNA substitution matrices. (Fabian, 2009-02-25)
 *
 *                A  B  C  D  E  F  G  H  I  K  L  M  0  N  P  Q  R  S 
 *                T  U  W  X  Y  Z (goodmeasuregoodmeasuregood) */
int ziAA2DNA[] = {0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		  3, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
const int ciMaxResID = 3;




void FastPairwiseAlign::pairAlign(const vector<int>* seq, int l1, int l2)
{
    int pot[8], i, j, l, m, limit, pos, tl1, vn1, vn2, flen, osptr, fs;
    int tv1, tv2, encrypt, subt1, subt2, rmndr;
    bool flag;
    int residue; 
    bool _DNAFlag = userParameters->getDNAFlag();
    int _ktup = userParameters->getKtup();
    int _maxAA = userParameters->getMaxAA();
    int _windowGap = userParameters->getWindowGap();
    int _window = userParameters->getWindow();
    int _signif = userParameters->getSignif();
        
    if (_DNAFlag)
    {
        for (i = 1; i <= _ktup; ++i)
        {
            pot[i] = (int)pow((double)4, (double)(i - 1));
        }
        limit = (int)pow((double)4, (double)_ktup);
    }
    else
    {
        for (i = 1; i <= _ktup; i++)
        {
            pot[i] = (int)pow((double)(_maxAA + 1), (double)(i - 1));
        }
        limit = (int)pow((double)(_maxAA + 1), (double)_ktup);
    }
    tl1 = (l1 + l2) - 1;

    for (i = 1; i <= tl1; ++i)
    {
        slopes[i] = displ[i] = 0;
        diagIndex[i] = i;
    }


    // increment diagonal score for each k_tuple match 

    for (i = 1; i <= limit; ++i)
    {
        vn1 = zzc[i];
        while (true)
        {
            if (!vn1)
            {
                break;
            }
            vn2 = zzd[i];
            while (vn2 != 0)
            {
                osptr = vn1 - vn2 + l2;
                ++displ[osptr];
                vn2 = zzb[vn2];
            }
            vn1 = zza[vn1];
        }
    }

    // choose the top SIGNIF diagonals 

    desQuickSort(displ, diagIndex, tl1);

    j = tl1 - _signif + 1;
    if (j < 1)
    {
        j = 1;
    }

    // flag all diagonals within WINDOW of a top diagonal

    for (i = tl1; i >= j; i--)
    if (displ[i] > 0)
    {
        pos = diagIndex[i];
        l = (1 > pos - _window) ? 1 : pos - _window;
        m = (tl1 < pos + _window) ? tl1 : pos + _window;
        for (; l <= m; l++)
        {
            slopes[l] = 1;
        }
    }

    for (i = 1; i <= tl1; i++)
    {
        displ[i] = 0;
    }


    currFrag = maxSoFar = 0;

    for (i = 1; i <= (l1 - _ktup + 1); ++i)
    {
        encrypt = flag = 0;
	if (_DNAFlag){
	  for (j = 1; j <= _ktup; ++j)
	    {
	      residue = ziAA2DNA[(*seq)[i + j - 1]];
	      if ((residue < 0) || (residue > ciMaxResID))
		{
		  flag = true;
		  break;
		}
	      encrypt += ((residue) * pot[j]);
	    }
	}
	else {
	  for (j = 1; j <= _ktup; ++j)
	    {
	      residue = (*seq)[i + j - 1];
	      if ((residue < 0) || (residue > _maxAA))
		{
		  flag = true;
		  break;
		}
	      encrypt += ((residue) * pot[j]);
	    }
	}

        if (flag)
        {
            continue;
        }
        ++encrypt;

        vn2 = zzd[encrypt];

        flag = false;
        while (true)
        {
            if (!vn2)
            {
                flag = true;
                break;
            }
            osptr = i - vn2 + l2;
            if (slopes[osptr] != 1)
            {
                vn2 = zzb[vn2];
                continue;
            }
            flen = 0;
            fs = _ktup;
            next = maxSoFar;

            
            //  A-loop
             

            while (true)
            {
                if (!next)
                {
                    ++currFrag;
                    if (currFrag >= 2 * _maxAlnLength)
                    {
                        utilityObject->info("(Partial alignment)");
                        vatend = 1;
                        return ;
                    }
                    displ[osptr] = currFrag;
                    putFrag(fs, i, vn2, flen);
                }
                else
                {
                    tv1 = accum[1][next];
                    tv2 = accum[2][next];
                    if (fragRelPos(i, vn2, tv1, tv2))
                    {
                        if (i - vn2 == accum[1][next] - accum[2][next])
                        {
                            if (i > accum[1][next] + (_ktup - 1))
                            {
                                fs = accum[0][next] + _ktup;
                            }
                            else
                            {
                                rmndr = i - accum[1][next];
                                fs = accum[0][next] + rmndr;
                            }
                            flen = next;
                            next = 0;
                            continue;
                        }
                        else
                        {
                            if (displ[osptr] == 0)
                            {
                                subt1 = _ktup;
                            }
                            else
                            {
                                if (i > accum[1][displ[osptr]] + (_ktup - 1))
                                {
                                    subt1 = accum[0][displ[osptr]] + _ktup;
                                }
                                else
                                {
                                    rmndr = i - accum[1][displ[osptr]];
                                    subt1 = accum[0][displ[osptr]] + rmndr;
                                }
                            }
                            subt2 = accum[0][next] - _windowGap + _ktup;
                            if (subt2 > subt1)
                            {
                                flen = next;
                                fs = subt2;
                            }
                            else
                            {
                                flen = displ[osptr];
                                fs = subt1;
                            }
                            next = 0;
                            continue;
                        }
                    }
                    else
                    {
                        next = accum[4][next];
                        continue;
                    }
                }
                break;
            }
            // End of Aloop
            vn2 = zzb[vn2];
        }
    }
    vatend = 0;
}

void FastPairwiseAlign::makePPtrs(vector<int>& tptr, vector<int>& pl, const vector<int>* seq, int length)
{
    int a[10];
    int i, j, limit, code;
    bool flag;
    int residue;
    int _ktup = userParameters->getKtup();
    int _maxAA = userParameters->getMaxAA();
    
    for (i = 1; i <= _ktup; i++)
    {
        a[i] = (int)pow((double)(_maxAA + 1), (double)(i - 1));
    }

    limit = (int)pow((double)(_maxAA + 1), (double)_ktup);
    if(limit >= (int)pl.size())
    {
        pl.resize(limit + 1);
    }
    if(length >= (int)tptr.size())
    {
        tptr.resize(length + 1);
    }
    
    for (i = 1; i <= limit; ++i)
    {
        pl[i] = 0;
    }
    for (i = 1; i <= length; ++i) // NOTE changed this
    {
        tptr[i] = 0;
    }

    for (i = 1; i <= (length - _ktup + 1); ++i)
    {
        code = 0;
        flag = false;
        for (j = 1; j <= _ktup; ++j)
        {
            residue = (*seq)[i + j - 1];
            if ((residue < 0) || (residue > _maxAA))
            {
                flag = true;
                break;
            }
            code += ((residue) * a[j]);
        }
        if (flag)
        {
            continue;
        }
        ++code;
        if (pl[code] != 0)
        {
            tptr[i] = pl[code];
        }
        pl[code] = i; //Ktuple code is at position i in sequence
    }
}

void FastPairwiseAlign::makeNPtrs(vector<int>& tptr,vector<int>& pl, const vector<int>* seq, int length)
{
    int pot[] =
    {
        0, 1, 4, 16, 64, 256, 1024, 4096
    };
    int i, j, limit, code;
    bool flag;
    int residue;
    int _ktup = userParameters->getKtup();
    
    limit = (int)pow((double)4, (double)_ktup);

    if(limit >= (int)pl.size())
    {
        pl.resize(limit + 1);
    }
    if(length >= (int)tptr.size())
    {
        tptr.resize(length + 1);
    }
    
    for (i = 1; i <= limit; ++i)
    {
        pl[i] = 0;
    }
    for (i = 1; i <= length; ++i)
    {
        tptr[i] = 0;
    }

    for (i = 1; i <= length - _ktup + 1; ++i)
    {
        code = 0;
        flag = false;
	for (j = 1; j <= _ktup; ++j)
	  {
	    residue = ziAA2DNA[(*seq)[i + j - 1]];
	    if ((residue < 0) || (residue > ciMaxResID))
	      {
		flag = true;
		break;
	      }
	    code += ((residue) * pot[j]);
	  }
        if (flag)
        {
            continue;
        }
        ++code;
        if (pl[code] != 0)
        {
            tptr[i] = pl[code];
        }
        pl[code] = i;
    }

}

void FastPairwiseAlign::putFrag(int fs, int v1, int v2, int flen)
{
    int end;
    accum[0][currFrag] = fs;
    accum[1][currFrag] = v1;
    accum[2][currFrag] = v2;
    accum[3][currFrag] = flen;

    if (!maxSoFar)
    {
        maxSoFar = 1;
        accum[4][currFrag] = 0;
        return ;
    }

    if (fs >= accum[0][maxSoFar])
    {
        accum[4][currFrag] = maxSoFar;
        maxSoFar = currFrag;
        return ;
    }
    else
    {
        next = maxSoFar;
        while (true)
        {
            end = next;
            next = accum[4][next];
            if (fs >= accum[0][next])
            {
                break;
            }
        }
        accum[4][currFrag] = next;
        accum[4][end] = currFrag;
    }
}

inline int FastPairwiseAlign::fragRelPos(int a1, int b1, int a2, int b2)
{
    int ret;
    int _ktup = userParameters->getKtup();
    
    ret = false;
    if (a1 - b1 == a2 - b2)
    {
        if (a2 < a1)
        {
            ret = true;
        }
    }
    else
    {
        if (a2 + _ktup - 1 < a1 && b2 + _ktup - 1 < b1)
        {
            ret = true;
        }
    }
    return ret;
}

void FastPairwiseAlign::desQuickSort(vector<int>& array1, vector<int>& array2, int arraySize)
{
    // Quicksort routine, adapted from chapter 4, page 115 of software tools 
    // by Kernighan and Plauger, (1986) 
    // Sort the elements of array1 and sort the 
    // elements of array2 accordingly
    int temp1, temp2;
    int p, pivlin;
    int i, j;
    int lst[50], ust[50]; // the maximum no. of elements must be
                            // < log(base2) of 50 

    lst[1] = 1;
    ust[1] = arraySize - 1;
    p = 1;

    while (p > 0)
    {
        if (lst[p] >= ust[p])
        {
            p--;
        }
        else
        {
            i = lst[p] - 1;
            j = ust[p];
            pivlin = array1[j];
            while (i < j)
            {
                for (i = i + 1; array1[i] < pivlin; i++)
                    ;
                for (j = j - 1; j > i; j--)
                    if (array1[j] <= pivlin)
                    {
                        break;
                    }
                if (i < j)
                {
                    temp1 = array1[i];
                    array1[i] = array1[j];
                    array1[j] = temp1;

                    temp2 = array2[i];
                    array2[i] = array2[j];
                    array2[j] = temp2;
                }
            }

            j = ust[p];

            temp1 = array1[i];
            array1[i] = array1[j];
            array1[j] = temp1;

            temp2 = array2[i];
            array2[i] = array2[j];
            array2[j] = temp2;

            if (i - lst[p] < ust[p] - i)
            {
                lst[p + 1] = lst[p];
                ust[p + 1] = i - 1;
                lst[p] = i + 1;
            }
            else
            {
                lst[p + 1] = i + 1;
                ust[p + 1] = ust[p];
                ust[p] = i - 1;
            }
            p = p + 1;
        }
    }
    return ;

}

}
