/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "ProfileBase.h"

namespace clustalw
{

/**
 * 
 * @param prfLen 
 * @param firstS 
 * @param lastS 
 */
ProfileBase::ProfileBase(int prfLen, int firstS, int lastS)
 : vwindow(5),
   vll(50),
   reducedGap(1.0),
   prfLength(prfLen),
   firstSeq(firstS),
   lastSeq(lastS)
{
    vlut.resize(numLetters, vector<int>(numLetters));
    for(int i = 0; i < numLetters; i++)
    {
        vlut[i][i] = 1;
    }
    
    pascarellaRes = "ACDEFGHKILMNPQRSTVYW";
    int pasprob[] = { 87, 87, 104, 69, 80, 139, 100, 104, 68, 79,
                     71, 137, 126, 93, 128, 124, 111, 75, 100, 77};
    pascarellaProb = vector<int>(pasprob, pasprob + 20);

    profile.resize(prfLength + 2, vector<int>(LENCOL + 2));
}


/**
 * 
 * @param seqArray 
 * @param gaps 
 * @param useStructPenalties 
 * @param gapPenaltyMask 
 * @param gapCoef 
 * @param lenCoef 
 */
void ProfileBase::calcGapCoeff(SeqArray* seqArray, vector<int>* gaps,  
         bool useStructPenalties, vector<char>* gapPenaltyMask, int gapCoef, int lenCoef)
{
    int c;
    int i, j;
    int is, ie;
    int _numSeq, val, pcid;
    vector<int> gapPos;
    vector<int> vWeight, resWeight, hydWeight;
    float scale;
    int _maxAA = userParameters->getMaxAA();
    
    _numSeq = lastSeq - firstSeq;
    if(_numSeq == 2)
    {
        pcid = static_cast<int>(percentId(&(*seqArray)[firstSeq], &(*seqArray)[firstSeq + 1]));
    }
    else 
        pcid = 0;

    for (j = 0; j < prfLength; j++)
    {
        (*gaps)[j] = 0;
    }

    // Check for a gap penalty mask

    if (useStructPenalties != false)
    {
        nVarPen = nHydPen = nPrefPen = true;
        gdist = 0;
    }
    else if (userParameters->getNoVarPenalties() == false && pcid > 60)
    {            
        nHydPen = nPrefPen = true;
        nVarPen = false;
    }
    else
    {
        nVarPen = true;
        nHydPen = userParameters->getNoHydPenalties();
        nPrefPen = userParameters->getNoPrefPenalties();
        gdist = userParameters->getGapDist();
    }                  
     
    for (i = firstSeq; i < lastSeq; i++)
    {
        // Include end gaps as gaps ?
        is = 0;
        ie = prfLength;
        if (userParameters->getUseEndGaps() == false && 
            userParameters->getEndGapPenalties() == false)
        {
            for (j = 0; j < prfLength; j++)
            {
                c = (*seqArray)[i][j];
                if ((c < 0) || (c > _maxAA))
                    is++;
                else
                    break;
            }
            for (j = prfLength - 1; j >= 0; j--)
            {
                c = (*seqArray)[i][j];
                if ((c < 0) || (c > _maxAA))
                    ie--;
                else
                    break;
            }
        }

        for (j = is; j < ie; j++)
        {
            if (((*seqArray)[i][j] < 0) || ((*seqArray)[i][j] > _maxAA))
            {
                (*gaps)[j]++;
            }
        }
    }

    int _DNAFlag = userParameters->getDNAFlag();
    
    if ((!_DNAFlag) && (nVarPen == false))
    {
        vWeight.resize(prfLength + 2);
        calcVPenalties(seqArray, &vWeight);
    }


    if ((!_DNAFlag) && (nPrefPen == false))
    {
        resWeight.resize(prfLength + 2);
        calcResidueSpecificPen(seqArray, &resWeight);
    }

    if ((!_DNAFlag) && (nHydPen == false))
    {
        hydWeight.resize(prfLength + 2);
        calcHydrophilicPen(seqArray, &hydWeight);
    }

    gapPos.resize(prfLength + 2);

    // mark the residues close to an existing gap (set gaps[i] = -ve)
  
    if (_DNAFlag || (gdist <= 0))
    {
        for (i = 0; i < prfLength; i++)
        { 
            gapPos[i] = (*gaps)[i];
        }
    }
    else
    {
        i = 0;
        while (i < prfLength)
        {
            if ((*gaps)[i] <= 0)
            {
                gapPos[i] = (*gaps)[i];
                i++;
            }
            else 
            {
                for (j = -gdist + 1; j < 0; j++)
                {
                    if ((i + j >= 0) && (i + j < prfLength) &&
                       (((*gaps)[i + j] == 0) || ((*gaps)[i + j] < j)))
                    {     
                        gapPos[i + j] = j;
                    }
                }
                while ((*gaps)[i] > 0)
                {
                    if (i >= prfLength)
                    { 
                        break;
                    }
                    gapPos[i] = (*gaps)[i];
                    i++;
                }
                for (j = 0; j < gdist; j++)
                {
                    if ((*gaps)[i + j] > 0)
                    { 
                        break;
                    }
                    if ((i + j >= 0) && (i + j < prfLength) && 
                       (((*gaps)[i + j] == 0) || ((*gaps)[i + j] < -j)))
                    {     
                         gapPos[i + j] = -j-1;
                    }
                }
                i += j;
            }
        }
    }
    
    for (j = 0; j < prfLength; j++)
    {          
        if (gapPos[j] <= 0)
        {

         // apply residue-specific and hydrophilic gap penalties.
         
            if (!_DNAFlag) 
            {
                profile[j + 1][GAPCOL] = localPenalty(gapCoef, j, &resWeight, &hydWeight,
                                                       &vWeight);
                profile[j+1][LENCOL] = lenCoef;
            }
            else 
            {
                profile[j + 1][GAPCOL] = gapCoef;
                profile[j + 1][LENCOL] = lenCoef;
            }

         // increase gap penalty near to existing gaps.
         
            if (gapPos[j] < 0)
            {
                profile[j + 1][GAPCOL] = static_cast<int>((profile[j + 1][GAPCOL] * (2.0 + 2.0 * (gdist + gapPos[j]) / gdist)));
            }
        }
        else
        {
            scale = ((float)(_numSeq - (*gaps)[j]) / (float)_numSeq) * reducedGap;
            profile[j + 1][GAPCOL] = static_cast<int>(scale * gapCoef);
            profile[j + 1][LENCOL] = static_cast<int>(0.5 * lenCoef);
        }

        // apply the gap penalty mask
        
        if (useStructPenalties != NONE)
        {
            val = (*gapPenaltyMask)[j] - '0';
            if (val > 0 && val < 10)
            {
                profile[j + 1][GAPCOL] *= val;
                profile[j + 1][LENCOL] *= val;
            }
        }

        // make sure no penalty is zero - even for all-gap positions
        
        if (profile[j + 1][GAPCOL] <= 0)
        { 
            profile[j + 1][GAPCOL] = 1;
        }
        if (profile[j + 1][LENCOL] <= 0)
        { 
            profile[j + 1][LENCOL] = 1;
        }
    }

    // set the penalties at the beginning and end of the profile
    if(userParameters->getEndGapPenalties() == true)
    {
        profile[0][GAPCOL] = gapCoef;
        profile[0][LENCOL] = lenCoef;
    }
    else
    {
        profile[0][GAPCOL] = 0;
        profile[0][LENCOL] = 0;
        profile[prfLength][GAPCOL] = 0;
        profile[prfLength][LENCOL] = 0;
    }
    if (userParameters->getDebug() > 0)
    {
        cout << "Opening penalties:\n";
        for(i = 0; i <= prfLength; i++)
        { 
            cout <<" " << i << ":" << profile[i][GAPCOL]<< " ";
        }
        cout << "\n";
    }
    if (userParameters->getDebug() > 0)
    {
        cout << "Extension penalties:\n";
        for(i = 0; i <= prfLength; i++) 
        {
            cout << i << ":" << profile[i][LENCOL] << " ";
        }
        cout << "\n";
    }    
}

/** **************************************************************************************
 *                               Protected functions                                     *
 *****************************************************************************************/
 

/**
 * 
 * @param aln 
 * @param weight 
 */
void ProfileBase::calcVPenalties(SeqArray* aln, vector<int>* weight)
{
    int ix1, ix2;
    int i, j, t;
    int _maxAA = userParameters->getMaxAA();
    int aminoCodeix1, aminoCodeix2;
    
    for (i = 0; i < prfLength; i++)
    {
        (*weight)[i] = 0;
        t = 0;
        for(j = i - vwindow;j < i + vwindow; j++)
        {
            if(j >= 0 && j < prfLength)
            {
                ix1 = (*aln)[firstSeq][j];
                ix2 = (*aln)[firstSeq + 1][j];
                if ((ix1 < 0) || (ix1 > _maxAA) || (ix2 < 0) || (ix2 > _maxAA))
                { 
                    continue;
                }
                aminoCodeix1 = userParameters->getAminoAcidCode(ix1);
                aminoCodeix2 = userParameters->getAminoAcidCode(ix2);
                (*weight)[i] += vlut[aminoCodeix1 - 'A'][aminoCodeix2 - 'A'];
                t++;
            } 
        }
        /* now we have a weight -t < w < t */
        (*weight)[i] +=t;
        if(t > 0)
            (*weight)[i] = ((*weight)[i] * 100)/(2 * t);
        else
            (*weight)[i] = 100;
        /* now we have a weight vll < w < 100 */
        if ((*weight)[i] < vll) 
            (*weight)[i] = vll;
    }
}


/**
 * 
 * @param aln 
 * @param weight 
 */
void ProfileBase::calcResidueSpecificPen(SeqArray* aln, vector<int>* weight)
{
    char ix;
    int j, k, _numSeq;
    int i;
    int _maxAA = userParameters->getMaxAA();
    int _pascarellaNumRes = pascarellaRes.size();
    
    _numSeq = lastSeq - firstSeq;
    for (i = 0; i < prfLength; i++)
    {
        (*weight)[i] = 0;
        for (k = firstSeq; k < lastSeq; k++)
        {
            for (j = 0; j < _pascarellaNumRes; j++)
            {
                ix = (*aln)[k][i];
                if ((ix < 0) || (ix > _maxAA)) 
                    continue;
                if (userParameters->getAminoAcidCode(ix) == pascarellaRes[j])
                {
                    (*weight)[i] += (180 - pascarellaProb[j]);
                    break;
                }
            }
        }
        (*weight)[i] /= _numSeq;
    }
}

/**
 * 
 * @param aln 
 * @param weight 
 */
void ProfileBase::calcHydrophilicPen(SeqArray* aln, vector<int>* weight)
{
    int res;
    int numHydResidues, j, k;
    int i, e, s;
    vector<int> hyd;
    float scale;
    int _maxAA = userParameters->getMaxAA();
    
    hyd.resize(prfLength + 2);
    string _hydResidues(userParameters->getHydResidues());
    numHydResidues = _hydResidues.size();
        
    for (i = 0; i < prfLength; i++)
    {
        (*weight)[i] = 0;
    }

    for (k = firstSeq; k < lastSeq; k++)
    {
        for (i = 0; i < prfLength; i++)
        {
            hyd[i] = 0;
            for (j = 0; j < numHydResidues; j++)
            {
                res = (*aln)[k][i];
                if ((res < 0) || (res > _maxAA)) 
                    continue;
                if (userParameters->getAminoAcidCode(res) == _hydResidues[j])
                {
                    hyd[i] = 1;
                    break;
                }
            }
        }
        i = 0;
        while (i < prfLength)
        {
            if (hyd[i] == 0) 
                i++;
            else
            {
                s = i;
                while ((hyd[i] != 0) && (i < prfLength))
                { 
                    i++;
                }
                e = i;
                if (e - s > 3)
                {
                    for (j = s; j < e; j++)
                    { 
                        (*weight)[j] += 100;
                    }
                }
            }
        }
    }

    scale = lastSeq - firstSeq;
    for (i = 0; i < prfLength; i++)
    {
        (*weight)[i] = static_cast<int>(((*weight)[i] / scale)); // Mark change 17-5-07
    }
}


/**
 * 
 * @param penalty 
 * @param n 
 * @param resWeight 
 * @param hydWeight 
 * @param vWeight 
 * @return 
 */
int ProfileBase::localPenalty(int penalty, int n, vector<int>* resWeight, 
                              vector<int>* hydWeight, vector<int>* vWeight)
{
    bool h = false;
    float gw;

    if (userParameters->getDNAFlag())
    { 
        return(1);
    }

    gw = 1.0;
    if (nVarPen == false)
    {
        gw *= (*vWeight)[n] / 100.0;
    }

    if (nHydPen == false)
    {
        if ((*hydWeight)[n] > 0)
        {
            gw *= 0.5;
            h = true;
        }
    }
    if ((nPrefPen == false) && (h == false))
    {
         gw *= ((*resWeight)[n] / 100.0);
    }

    gw *= penalty;
    return (int)gw;

}

/** *********************************************************************
 * Note dont think this will work. Dont have -3. Need to use lengths!   *
 ************************************************************************/

/**
 * 
 * @param s1 
 * @param s2 
 * @return 
 */
float ProfileBase::percentId(vector<int>* s1, vector<int>* s2)
{
    int i;
    int count, total;
    float score;

    count = total = 0;
    for (i = 0; i < prfLength; i++) 
    {
        if (((*s1)[i] >= 0) && ((*s1)[i] < userParameters->getMaxAA())) 
        {
            total++;
            if ((*s1)[i] == (*s2)[i])
            { 
                count++;
            }
        }
        if ((*s1)[i]==(-3) || (*s2)[i]==(-3))
        { 
            break; // I dont have -3 at the end!
        }

    }

    if(total == 0)
    { 
        score = 0;
    }
    else
    {
        score = 100.0 * (float)count / (float)total;
    }
    return (score);
}
                         

}
