/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <stdio.h>
#include <string>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include "SubMatrix.h"
#include "matrices.h"
#include "../general/InvalidCombination.cpp"
#include "../general/debuglogObject.h"

namespace clustalw
{
using namespace std;

/**
 * 
 * @param log 
 */
SubMatrix::SubMatrix()
: sizenAAMatrix(276),
  sizeDNAMatrix(153),
  matrixAvgScore(0),
  QTDNAHistMatNum(DNAIUB),
  QTAAHistMatNum(AAHISTGONNETPAM250),
  QTsegmentDNAMatNum(DNAIUB),
  QTsegmentAAMatNum(QTAASEGGONNETPAM250)
{
    userSeries = false;
         
    setUpCrossReferences();

    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("Creating the SubMatrix object\n");
        }
    #endif
    
    try
    {              
        /* Set up the vectors with the matrices defined in matrices.h 
         * The matrices are intially defined as a short array, as there is no way
         * to intiailize a vector with a {....} list.  
         */
        blosum30mtVec = new Matrix(blosum30mt, blosum30mt + sizenAAMatrix);
        blosum40mtVec = new Matrix(blosum40mt, blosum40mt + sizenAAMatrix);
        blosum45mtVec = new Matrix(blosum45mt, blosum45mt + sizenAAMatrix);
        blosum62mt2Vec = new Matrix(blosum62mt2, blosum62mt2 + sizenAAMatrix);
        blosum80mtVec = new Matrix(blosum80mt, blosum80mt + sizenAAMatrix);
        pam20mtVec = new Matrix(pam20mt, pam20mt + sizenAAMatrix);
        pam60mtVec = new Matrix(pam60mt, pam60mt + sizenAAMatrix);
        pam120mtVec = new Matrix(pam120mt, pam120mt + sizenAAMatrix);
        pam350mtVec = new Matrix(pam350mt, pam350mt + sizenAAMatrix);
        idmatVec = new Matrix(idmat, idmat + sizenAAMatrix);
        gon40mtVec = new Matrix(gon40mt, gon40mt + sizenAAMatrix);
        gon80mtVec = new Matrix(gon80mt, gon80mt + sizenAAMatrix);
        gon120mtVec = new Matrix(gon120mt, gon120mt + sizenAAMatrix);
        gon160mtVec = new Matrix(gon160mt, gon160mt + sizenAAMatrix);
        gon250mtVec = new Matrix(gon250mt, gon250mt + sizenAAMatrix);
        gon350mtVec = new Matrix(gon350mt, gon350mt + sizenAAMatrix);
        clustalvdnamtVec = new Matrix(clustalvdnamt, clustalvdnamt + sizeDNAMatrix);
        swgapdnamtVec = new Matrix(swgapdnamt, swgapdnamt + sizeDNAMatrix);
        
        /*
         * Set up the vectors for user defined types.
         * Probably dont need to do this, as they may not be used. It would be better
         * to initialise their size if thye are used.
         */
        userMat.resize(NUMRES * NUMRES);
        pwUserMat.resize(NUMRES * NUMRES);
        userDNAMat.resize(NUMRES * NUMRES);
        pwUserDNAMat.resize(NUMRES * NUMRES);
        QTscoreUserMatrix.resize(NUMRES * NUMRES);
        QTscoreUserDNAMatrix.resize(NUMRES * NUMRES);
        QTsegmentDNAMatrix.resize(NUMRES * NUMRES);
        QTsegmentAAMatrix.resize(NUMRES * NUMRES);
        
        userMatSeries.resize(MAXMAT); // Maximum number of matrices        
        vector<Matrix>::iterator firstM = userMatSeries.begin();
        vector<Matrix>::iterator lastM = userMatSeries.end();

        while(firstM != lastM)
        {
            firstM->resize(NUMRES * NUMRES);
            firstM++;
        }
        
        // Maybe I should put this in with the other xref intialisations!
        AAXrefseries.resize(MAXMAT);
        vector<Xref>::iterator firstX = AAXrefseries.begin();
        vector<Xref>::iterator lastX = AAXrefseries.end();
        
        while(firstX != lastX)
        {
            firstX->resize(NUMRES + 1);
            firstX++;
        }
        
        // Set the defaults.
        matrixNum = 3;
        matrixName = new string("gonnet");
        DNAMatrixNum = 1;
        DNAMatrixName = new string("iub");
        pwMatrixNum = 3;
        pwMatrixName = new string("gonnet");
        pwDNAMatrixNum = 1;
        pwDNAMatrixName = new string("iub");
    }
    catch(const exception &ex)
    {
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue\n";
        throw 1;
    }  
}

void SubMatrix::setValuesToDefault()
{
    matrixAvgScore = 0;
    QTDNAHistMatNum = DNAIUB;
    QTAAHistMatNum = AAHISTGONNETPAM250;
    QTsegmentDNAMatNum = DNAIUB;
    QTsegmentAAMatNum = QTAASEGGONNETPAM250;
    userSeries = false;
    matrixNum = 3;
    DNAMatrixNum = 1;
    pwMatrixNum = 3;
    pwDNAMatrixNum = 1;      
}

/**
 * The destructor frees up any dynamically allocated memory.
 */
SubMatrix::~SubMatrix()
{
    delete blosum30mtVec;
    delete blosum40mtVec;
    delete blosum45mtVec;
    delete blosum62mt2Vec;
    delete blosum80mtVec;
    delete pam20mtVec;
    delete pam60mtVec;
    delete pam120mtVec;
    delete pam350mtVec;
    delete idmatVec;
    delete gon40mtVec;
    delete gon80mtVec;
    delete gon120mtVec;
    delete gon160mtVec;
    delete gon250mtVec;
    delete gon350mtVec;
    delete clustalvdnamtVec;
    delete swgapdnamtVec;
    delete matrixName;
    delete DNAMatrixName;
    delete pwMatrixName;
    delete pwDNAMatrixName;
}

/**
 * This function sets up the initial cross references. 
 */
void SubMatrix::setUpCrossReferences()
{
    char c1, c2;
    short i, j, maxRes;  
    defaultAAXref.resize(NUMRES + 1);
    defaultDNAXref.resize(NUMRES + 1);

    string aminoAcidOrder = "ABCDEFGHIKLMNPQRSTVWXYZ";
    string nucleicAcidOrder = "ABCDGHKMNRSTUVWXY"; 
    /* 
     * I also need to resize the user defined xrefs.
     */
    DNAXref.resize(NUMRES + 1);
    AAXref.resize(NUMRES + 1);
    pwAAXref.resize(NUMRES + 1);
    pwDNAXref.resize(NUMRES + 1);
    QTscoreXref.resize(NUMRES + 1);
    QTscoreDNAXref.resize(NUMRES + 1);
    QTsegmentDNAXref.resize(NUMRES + 1);
    QTsegmentAAXref.resize(NUMRES + 1);
    /*
     * set up cross-reference for default matrices hard-coded in matrices.h
     */
    for (i = 0; i < NUMRES; i++)
    {
        defaultAAXref[i] = -1;
    }
    for (i = 0; i < NUMRES; i++)
    {
        defaultDNAXref[i] = -1;
    }

    maxRes = 0;

    for (i = 0; (c1 = aminoAcidOrder[i]); i++)
    {
        for (j = 0; (c2 = userParameters->getAminoAcidCode(j)); j++)
        {
            if (c1 == c2)
            {
                defaultAAXref[i] = j;

                maxRes++;
                break;
            }
        }
        if ((defaultAAXref[i] ==  - 1) && (aminoAcidOrder[i] != '*'))
        {
            utilityObject->error("residue %c in matrices.h is not recognised",
                aminoAcidOrder[i]);
        }
    } 

    maxRes = 0;
    for (i = 0; (c1 = nucleicAcidOrder[i]); i++)
    {
        for (j = 0; (c2 = userParameters->getAminoAcidCode(j)); j++)
        {
            if (c1 == c2)
            {
                defaultDNAXref[i] = j;
                maxRes++;
                break;
            }
        }
        if ((defaultDNAXref[i] ==  - 1) && (nucleicAcidOrder[i] != '*'))
        {
            utilityObject->error("nucleic acid %c in matrices.h is not recognised",
                nucleicAcidOrder[i]);
        }
    }
    
}

/**
 * The function getPairwiseMatrix is called by the user to get the correct sub matrix for the
 * pairwise alignment stage. It calls getMatrix to actually calculate the matrix.
 * This function provides an interface for the user. It allows the user to get the matrix 
 * they wish to use for the pairwise alignment part. 
 * @param matrix[][] 
 * @param scale 
 * @param matAvg 
 * @return 
 */
int SubMatrix::getPairwiseMatrix(int matrix[NUMRES][NUMRES], PairScaleValues& scale,
                                 int& matAvg)
{
    int _maxRes; // Local copy.
    /* Pointers to Matrix and xref to use in calculation */
    Matrix* _matPtr;
    Xref* _matXref;
    
    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("In the function getPairwiseMatrix: \n");
        }    
    #endif
    
    string matrixPointer;
    string xrefPointer;
        
    #ifdef OS_MAC
        scale.intScale = 10;
    #else
        scale.intScale = 100;
    #endif
    scale.gapOpenScale = scale.gapExtendScale = 1.0;
    
    if (userParameters->getDNAFlag())
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                string msg = "    (DNA AND Pairwise) " + *pwDNAMatrixName + "\n";
                logObject->logMsg(msg);
            }
        #endif
        if (pwDNAMatrixName->compare("iub") == 0)
        {
            matrixPointer = "swgapdnamtVec"; xrefPointer = "defaultDNAXref";
            _matPtr = swgapdnamtVec;
            _matXref = &defaultDNAXref;
        }
        else if (pwDNAMatrixName->compare("clustalw") == 0)
        {
            matrixPointer = "clustalvdnamtVec"; xrefPointer = "defaultDNAXref";
            _matPtr = clustalvdnamtVec;
            _matXref = &defaultDNAXref;
            scale.gapOpenScale = 0.6667;
            scale.gapExtendScale = 0.751;
        }
        else
        {
            matrixPointer = "pwUserDNAMat"; xrefPointer = "pwDNAXref";
            _matPtr = &pwUserDNAMat;
            _matXref = &pwDNAXref;
        }
        _maxRes = getMatrix(_matPtr, _matXref, matrix, true, scale.intScale);
        
        if (_maxRes == 0)
        {
            return ((int) -1);
        }
        float _transitionWeight = userParameters->getTransitionWeight();
        // Mark change 17-5-07
        matrix[0][4] = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[4][0] = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[2][11] = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[11][2] = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[2][12] = static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[12][2] = static_cast<int>(_transitionWeight * matrix[0][0]);
        
    }
    else
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                string msg = "    (Protein AND Pairwise) " + *pwMatrixName + "\n";
                logObject->logMsg(msg);
            }
        #endif
        
        if (pwMatrixName->compare("blosum") == 0)
        {
            matrixPointer = "blosum30mtVec"; xrefPointer = "defaultAAXref";
            _matPtr = blosum30mtVec;
            _matXref = &defaultAAXref;
        }
        else if (pwMatrixName->compare("pam") == 0)
        {
            matrixPointer = "pam350mtVec"; xrefPointer = "defaultAAXref";
            _matPtr = pam350mtVec;
            _matXref = &defaultAAXref;
        }
        else if (pwMatrixName->compare("gonnet") == 0)
        {
            matrixPointer = "gon250mtVec"; xrefPointer = "defaultAAXref";
            _matPtr = gon250mtVec;
            scale.intScale /= 10;
            _matXref = &defaultAAXref;
        }
        else if (pwMatrixName->compare("id") == 0)
        {
            matrixPointer = "idmatVec"; xrefPointer = "defaultAAXref";
            _matPtr = idmatVec;
            _matXref = &defaultAAXref;
        }
        else
        {
            matrixPointer = "pwUserMat"; xrefPointer = "pwAAXref";
            _matPtr = &pwUserMat;
            _matXref = &pwAAXref;
        }

        _maxRes = getMatrix(_matPtr, _matXref, matrix, true, scale.intScale);
        
        if (_maxRes == 0)
        {
            return ((int) -1);
        }
    }
    
    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            ostringstream outs;
            outs << "    Called getMatrix with "
                 << matrixPointer << " and " << xrefPointer << ".\n"
                 << "    intScale = " << scale.intScale << ", gapOpenScale = "
                 << scale.gapOpenScale << ", gapExtendScale = " << scale.gapExtendScale 
                 << "\n\n";
            logObject->logMsg(outs.str());
        }
    #endif
    matAvg = matrixAvgScore;  
    return _maxRes; 
}

/**
 * The function getProfileAlignMatrix provides an interface for the user to get
 * the matrix to be used in the profile alignment. This depends on the matrix series
 * that was chosen, and also on the percent identity.
 * @param matrix[][] 
 * @param pcid 
 * @param minLen 
 * @param scaleParam 
 * @param matAvg 
 * @return 
 */
int SubMatrix::getProfileAlignMatrix(int matrix[NUMRES][NUMRES], double pcid, int minLen, 
                                  PrfScaleValues& scaleParam, int& matAvg)
{
    bool found = false;
    bool errorGiven = false;
    bool _negMatrix = userParameters->getUseNegMatrix();
    int i = 0, j = 0;
    int _maxRes = 0;
    scaleParam.intScale = 100;
    string matrixPointer;
    string xrefPointer;
    
    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("In the function getProfileAlignMatrix: \n");
        }    
    #endif    
    if (userParameters->getDNAFlag())
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                ostringstream outs;
                outs << "    (DNA AND Multiple align) "<< DNAMatrixName->c_str() << "\n";
                logObject->logMsg(outs.str());
            }    
        #endif        
        scaleParam.scale = 1.0;
        if (DNAMatrixName->compare("iub") == 0)
        {
            _matPtr = swgapdnamtVec;
            _matXref = &defaultDNAXref;
            matrixPointer = "swgapdnamtVec"; xrefPointer = "defaultDNAXref";
        }
        else if (DNAMatrixName->compare("clustalw") == 0)
        {
            _matPtr = clustalvdnamtVec;
            _matXref = &defaultDNAXref;
            scaleParam.scale = 0.66;
            matrixPointer = "clustalvdnamtVec"; xrefPointer = "defaultDNAXref";
        }
        else
        {
            _matPtr = &userDNAMat;
            _matXref = &DNAXref;
            matrixPointer = "userDNAMat"; xrefPointer = "DNAXref";
        }
        
        _maxRes = getMatrix(_matPtr, _matXref, matrix, _negMatrix,
                            static_cast<int>(scaleParam.intScale)); // Mark change 17-5-07
        if (_maxRes == 0)
        {
            return ((int) - 1);
        }
        
        float _transitionWeight = userParameters->getTransitionWeight();
        // fix suggested by Chanan Rubin at Compugen 
        matrix[(*_matXref)[0]][(*_matXref)[4]] = 
                            static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[(*_matXref)[4]][(*_matXref)[0]] = 
                            static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[(*_matXref)[2]][(*_matXref)[11]] = 
                            static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[(*_matXref)[11]][(*_matXref)[2]] = 
                            static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[(*_matXref)[2]][(*_matXref)[12]] = 
                            static_cast<int>(_transitionWeight * matrix[0][0]);
        matrix[(*_matXref)[12]][(*_matXref)[2]] = 
                            static_cast<int>(_transitionWeight * matrix[0][0]);

    }
    else // Amino acid alignment!!!!
    {   
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                ostringstream outs;
                outs << "    (Protein AND Multiple align) "<< matrixName->c_str() << "\n";
                logObject->logMsg(outs.str());
            }    
        #endif   
             
        scaleParam.scale = 0.75;
        if (matrixName->compare("blosum") == 0)
        {
            scaleParam.scale = 0.75;
            if (_negMatrix || userParameters->getDistanceTree() == false)
            {
                _matPtr = blosum40mtVec;
                matrixPointer = "blosum40mtVec";
            }
            else if (pcid > 80.0)
            {
                _matPtr = blosum80mtVec;
                matrixPointer = "blosum80mtVec";
            }
            else if (pcid > 60.0)
            {
                _matPtr = blosum62mt2Vec;
                matrixPointer = "blosum62mt2Vec";
            }
            else if (pcid > 40.0)
            {
                _matPtr = blosum45mtVec;
                matrixPointer = "blosum45mtVec";
            }
            else if (pcid > 30.0)
            {
                scaleParam.scale = 0.5;
                _matPtr = blosum45mtVec;
                matrixPointer = "blosum45mtVec";
            }
            else if (pcid > 20.0)
            {
                scaleParam.scale = 0.6;
                _matPtr = blosum45mtVec;
                matrixPointer = "blosum45mtVec";
            }
            else
            {
                scaleParam.scale = 0.6;
                _matPtr = blosum30mtVec;
                matrixPointer = "blosum30mtVec";
            }
            _matXref = &defaultAAXref;
            xrefPointer = "defaultAAXref";
        }
        else if (matrixName->compare("pam") == 0)
        {
            scaleParam.scale = 0.75;
            if (_negMatrix || userParameters->getDistanceTree() == false)
            {
                _matPtr = pam120mtVec;
                matrixPointer = "pam120mtVec";
            }
            else if (pcid > 80.0)
            {
                _matPtr = pam20mtVec;
                matrixPointer = "pam20mtVec";
            }
            else if (pcid > 60.0)
            {
                _matPtr = pam60mtVec;
                matrixPointer = "pam60mtVec";
            }
            else if (pcid > 40.0)
            {
                _matPtr = pam120mtVec;
                matrixPointer = "pam120mtVec";
            }
            else
            {
                _matPtr = pam350mtVec;
                matrixPointer = "pam350mtVec";
            }
            _matXref = &defaultAAXref;
            xrefPointer = "defaultAAXref";
        }
        else if (matrixName->compare("gonnet") == 0)
        {
            scaleParam.scale /= 2.0;
            if (_negMatrix || userParameters->getDistanceTree() == false)
            {
                _matPtr = gon250mtVec;
                matrixPointer = "gon250mtVec";
            }
            else if (pcid > 35.0)
            {
                _matPtr = gon80mtVec;
                scaleParam.scale /= 2.0;
                matrixPointer = "gon80mtVec";
            }
            else if (pcid > 25.0)
            {
                if (minLen < 100)
                {
                    _matPtr = gon250mtVec;
                    matrixPointer = "gon250mtVec";
                }
                else
                {
                    _matPtr = gon120mtVec;
                    matrixPointer = "gon120mtVec";
                }
            }
            else
            {
                if (minLen < 100)
                {
                    _matPtr = gon350mtVec;
                    matrixPointer = "gon350mtVec";
                }
                else
                {
                    _matPtr = gon160mtVec;
                    matrixPointer = "gon160mtVec";
                }
            }
            _matXref = &defaultAAXref;
            xrefPointer = "defaultAAXref";
            scaleParam.intScale /= 10;
        }
        else if (matrixName->compare("id") == 0)
        {
            _matPtr = idmatVec;
            _matXref = &defaultAAXref;
            xrefPointer = "defaultAAXref";
            matrixPointer = "idmatVec";
        }
        else if (userSeries)
        {
            _matPtr = NULL;
            found = false;
            for (i = 0; i < matSeries.nmat; i++)
            {
                if (pcid >= matSeries.mat[i].llimit && pcid <=
                    matSeries.mat[i].ulimit)
                {
                    j = i;
                    found = true;
                    break;
                }
            }
            if (found == false)
            {
                if (!errorGiven)
                {
                    utilityObject->warning(
                        "\nSeries matrix not found for sequence percent identity = %d.\n""(Using first matrix in series as a default.)\n""This alignment may not be optimal!\n""SUGGESTION: Check your matrix series input file and try again.", (int)pcid);
                }
                errorGiven = true;
                j = 0;
            }

            _matPtr = matSeries.mat[j].matptr;
            _matXref = matSeries.mat[j].AAXref;
            // this gives a scale of 0.5 for pcid=llimit and 1.0 for pcid=ulimit 
            scaleParam.scale = 0.5 + (pcid - matSeries.mat[j].llimit) / (
                (matSeries.mat[j].ulimit - matSeries.mat[j].llimit) *2.0);
            xrefPointer = "matSeries.mat[j].AAXref";
            matrixPointer = "matSeries.mat[j].matptr";
        }
        else
        {
            _matPtr = &userMat;
            _matXref = &AAXref;
            xrefPointer = "AAXref";
            matrixPointer = "userMat";
        }
        
        _maxRes = getMatrix(_matPtr, _matXref, matrix, _negMatrix,
                            static_cast<int>(scaleParam.intScale));
        if (_maxRes == 0)
        {
            cerr << "Error: matrix " << matrixName << " not found\n";
            return ( - 1);
        }       
    }
    
    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            ostringstream outs;
            outs << "    Called getMatrix with "
                 << matrixPointer << " and " << xrefPointer << ".\n"
                 << "    intScale = " << scaleParam.intScale << ", scale = "
                 << scaleParam.scale << ", pcid = " << pcid << "\n\n";
            logObject->logMsg(outs.str());
        }    
    #endif    
    
    matAvg = matrixAvgScore;
    return _maxRes;
}

/**
 * The function getMatrix is what the other parts of the code call to get a useable 
 * substitution matrix. This is stored in matrix[NUMRES][NUMRES].
 * @param matptr 
 * @param xref 
 * @param matrix[][] 
 * @param negFlag 
 * @param scale 
 * @return 
 */
int SubMatrix::getMatrix(Matrix* matptr, Xref* xref, int matrix[NUMRES][NUMRES], 
                         bool negFlag, int scale, bool minimise)
{
    int ggScore = 0;
    int grScore = 0;
    int i, j, k, ix = 0;
    int ti, tj;
    int maxRes;
    int av1, av2, av3, min, max;  

    for (i = 0; i < NUMRES; i++)
    {
        for (j = 0; j < NUMRES; j++)
        {
            matrix[i][j] = 0;
        }
    }
    
    ix = 0;
    maxRes = 0;
    for (i = 0; i <= userParameters->getMaxAA(); i++)
    {
        ti = (*xref)[i];
        for (j = 0; j <= i; j++)
        {
            tj = (*xref)[j];
            if ((ti !=  - 1) && (tj !=  - 1))
            {
                k = (*matptr)[ix];
                if (ti == tj)
                {
                    matrix[ti][ti] = k * scale;
                    maxRes++;
                }
                else
                {
                    matrix[ti][tj] = k * scale;
                    matrix[tj][ti] = k * scale;
                }
                ix++;
            }
        }
    }

    --maxRes;

    av1 = av2 = av3 = 0;
    for (i = 0; i <= userParameters->getMaxAA(); i++)
    {
        for (j = 0; j <= i; j++)
        {
            av1 += matrix[i][j];
            if (i == j)
            {
                av2 += matrix[i][j];
            }
            else
            {
                av3 += matrix[i][j];
            }
        }
    }

    av1 /= (maxRes * maxRes) / 2;
    av2 /= maxRes;
    av3 = static_cast<int>(av3 / (((float)(maxRes * maxRes - maxRes)) / 2));
    matrixAvgScore =  - av3;

    min = max = matrix[0][0];
    for (i = 0; i <= userParameters->getMaxAA(); i++)
    for (j = 1; j <= i; j++)
    {
        if (matrix[i][j] < min)
        {
            min = matrix[i][j];
        }
        if (matrix[i][j] > max)
        {
            max = matrix[i][j];
        }
    }
    //cout << "MAX = " << max << "\n";
    if(!minimise)
    {
        /*
            if requested, make a positive matrix - add -(lowest score) to every entry
        */
        if (negFlag == false)
        {
            if (min < 0)
            {
                for (i = 0; i <= userParameters->getMaxAA(); i++)
                {
                    ti = (*xref)[i];
                    if (ti !=  - 1)
                    {
                        for (j = 0; j <= userParameters->getMaxAA(); j++)
                        {
                            tj = (*xref)[j];

                            if (tj !=  - 1)
                            {
                                matrix[ti][tj] -= min;
                            }
                        }
                    }
                }
            }
        }
    
        // local copies of the gap positions
        int _gapPos1 = userParameters->getGapPos1();
        int _gapPos2 = userParameters->getGapPos2();

        for (i = 0; i < userParameters->getGapPos1(); i++)
        {
            matrix[i][_gapPos1] = grScore;
            matrix[_gapPos1][i] = grScore;
            matrix[i][_gapPos2] = grScore;
            matrix[_gapPos2][i] = grScore;
        }
        matrix[_gapPos1][_gapPos1] = ggScore;
        matrix[_gapPos2][_gapPos2] = ggScore;
        matrix[_gapPos2][_gapPos1] = ggScore;
        matrix[_gapPos1][_gapPos2] = ggScore;
    }
    else
    {
        // DO THE SAGA MATRIX
        for (i = 0; i <= userParameters->getMaxAA(); i++)
        {
            for (j = 0; j <= userParameters->getMaxAA(); j++)
            {
                matrix[i][j] = max - matrix[i][j];
            }
        }    
    }
    maxRes += 2;

    return (maxRes);
}

/**
 * The function getUserMatFromFile is used to read in a user defined matrix.
 * @param str
 * @param alignResidueType
 * @param alignType
 * @return
 */
bool SubMatrix::getUserMatFromFile(char *str, int alignResidueType, int alignType)
{
    int maxRes;

    FILE *infile;
    // Need to check if the values are a valid combination!
    checkResidueAndAlignType(alignResidueType, alignType);

    if(userParameters->getMenuFlag())
    {
        utilityObject->getStr(string("Enter name of the matrix file"), line2);
    }
    else
    {
        line2 = string(str);
    }

    if(line2.size() == 0)
    {
        return false;
    }

    if((infile = fopen(line2.c_str(), "r")) == NULL)
    {
        utilityObject->error("Cannot find matrix file [%s]", line2.c_str());
        return false;
    }

    strcpy(str, line2.c_str());
    // Find out which part of the code we are reading the matrix in for.
    mat = getUserMatAddress(alignResidueType, alignType);
    xref = getUserXrefAddress(alignResidueType, alignType);

    if ((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        // Try read in a matrix series.
        maxRes = readMatrixSeries(str, userMat, AAXref);
    }
    else // Read in a single matrix!!!
    {
        maxRes = readUserMatrix(str, *mat, *xref);
    }

    if (maxRes <= 0) return false;

    return true;
}

bool SubMatrix::getUserMatFromR(Rcpp::NumericMatrix substitutionMatrix, int alignResidueType, int alignType) {
	int maxRes;
  //  // Need to check if the values are a valid combination!
    checkResidueAndAlignType(alignResidueType, alignType);

    // Find out which part of the code we are reading the matrix in for.
    mat = getUserMatAddress(alignResidueType, alignType);
    xref = getUserXrefAddress(alignResidueType, alignType);
    
    if ((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        // Try read in a matrix series.
        maxRes = readMatrixSeriesFromR(substitutionMatrix, userMat, AAXref);
    }
    else // Read in a single matrix!!!
    {
        maxRes = readUserMatrixFromR(substitutionMatrix, *mat, *xref);
    }
    
    if (maxRes <= 0) {
    	return false;
    }

    return true;
}

/**
 * The function compareMatrices is used to compare 2 matrices that have been read in.
 * It will compare them in a given region. It will not compare all of them, as some of it
 * will be random memory.
 * @param mat1[][] 
 * @param mat2[][] 
 */
void SubMatrix::compareMatrices(int mat1[NUMRES][NUMRES], int mat2[NUMRES][NUMRES])
{
    int same = 1;
    for(int row = 0; row < NUMRES; row++)
    {
        for(int col = 0; col < NUMRES; col++)
        {
            if(mat1[row][col] != mat2[row][col])
            {
                same = 0;
                cout << "The row is " << row << ". The column is " << col << endl;
                break; // It is not the same. End the loop.
            }
        }
    }
    
    if(same == 0)
    {
        cout << "It was not the same\n";
    }
    else
    {
        cout << "It is the same\n";
    }
}

/**
 * This function is simply to display the results of the getMatrix function.
 * This is so that I can compare it to the original clustal version of it.
 * @param mat[][] 
 */
void SubMatrix::printGetMatrixResults(int mat[NUMRES][NUMRES]) {
    ofstream outfile("getmatrix.out");
    
    if(!outfile)
        cerr<<"oops failed to open !!!\n";
    
    for(int row = 0; row < NUMRES; row++)
    {
        for(int col = 0; col < NUMRES; col++)
        {
            if((mat[row][col] > 9) || (mat[row][col] < 0))
            {
                outfile <<" "<< mat[row][col]<<",";
            }
            else
            {
                outfile <<"  "<< mat[row][col]<<",";
            }
        }
        outfile<<"\n";
    }
}

/**
 * This function is from the old interface. It simply calls the readMatrixSeries
 * function. It seems to only be called from the commandline part.
 * @param str 
 * @return 
 */
bool SubMatrix::getUserMatSeriesFromFile(char *str)
{
    int maxRes;

    FILE *infile;

    if(userParameters->getMenuFlag()) // Check if we are using menu!
    {
        utilityObject->getStr(string("Enter name of the matrix file"), line2);
    }
    else
    {
        //strcpy(lin2,str);
        line2 = string(str);
    }

    if(line2.size() == 0) return false;

    if((infile = fopen(line2.c_str(), "r"))==NULL)
    {
        utilityObject->error("Cannot find matrix file [%s]",line2.c_str());
        return false;
    }

    strcpy(str, line2.c_str());

    maxRes = readMatrixSeries(str, userMat, AAXref);
    if (maxRes <= 0) return false;

    return true;
}

/*
 * The function readMatrixSeries is used to read in a series of matrices from a file.
 * It calls readUserMatrix to read in the individual matrices. The matrices are stored
 * in userMatSeries. 
 */
int SubMatrix::readMatrixSeries(const char *fileName, Matrix& userMat, Xref& xref)
{
    FILE *fd = NULL;
    char mat_fileName[FILENAMELEN];
    char inline1[1024];
    int maxRes = 0;
    int nmat;
    int n, llimit, ulimit;

    if (fileName[0] == '\0')
    {
        utilityObject->error("comparison matrix not specified");
        return ((int)0);
    }
    if ((fd = fopen(fileName, "r")) == NULL)
    {
        utilityObject->error("cannot open %s", fileName);
        return ((int)0);
    }

    /* check the first line to see if it's a series or a single matrix */
    while (fgets(inline1, 1024, fd) != NULL)
    {
        if (commentline(inline1))
        {
            continue;
        }
        if (utilityObject->lineType(inline1, "CLUSTAL_SERIES"))
        {
            userSeries = true;
        }
        else
        {
            userSeries = false;
        }
        break;
    }

    /* it's a single matrix */
    if (userSeries == false)
    {
        fclose(fd);
        maxRes = readUserMatrix(fileName, userMat, xref);
        return (maxRes);
    }

    /* it's a series of matrices, find the next MATRIX line */
    nmat = 0;
    matSeries.nmat = 0;
    while (fgets(inline1, 1024, fd) != NULL)
    {
        if (commentline(inline1))
        {
            continue;
        }
        if (utilityObject->lineType(inline1, "MATRIX")) // Have found a matrix
        {
            if (sscanf(inline1 + 6, "%d %d %s", &llimit, &ulimit, mat_fileName)
                != 3)
            {
                utilityObject->error("Bad format in file %s\n", fileName);
                fclose(fd);
                return ((int)0);
            }
            if (llimit < 0 || llimit > 100 || ulimit < 0 || ulimit > 100)
            {
                utilityObject->error("Bad format in file %s\n", fileName);
                fclose(fd);
                return ((int)0);
            }
            if (ulimit <= llimit)
            {
                utilityObject->error("in file %s: lower limit is greater than upper (%d-%d)\n",
                    fileName, llimit, ulimit);
                fclose(fd);
                return ((int)0);
            }
            
            n = readUserMatrix(mat_fileName, userMatSeries[nmat],
                AAXrefseries[nmat]);
            
            //cout << "Read in matrix number " << nmat << "\n"; // NOTE Testing!!!!! 
            char nameOfFile[] = "matrix"; 
            printInFormat(userMatSeries[nmat], nameOfFile);
            if (n <= 0)
            {
                utilityObject->error("Bad format in matrix file %s\n", mat_fileName);
                fclose(fd);
                return ((int)0);
            }
            matSeries.mat[nmat].llimit = llimit;
            matSeries.mat[nmat].ulimit = ulimit;
            matSeries.mat[nmat].matptr = &userMatSeries[nmat];
            matSeries.mat[nmat].AAXref = &AAXrefseries[nmat];
            nmat++;
            
            if(nmat >= MAXMAT)
            {
               // We have read in all the matrices that we can read into this vector
               // Write a message to the screen, and break out of loop.
               cerr << "The matrix series file has more entries than allowed in \n"
                    << "a user defined series. The most that are allowed is "
                    << MAXMAT << ".\n"
                    << "The first " << MAXMAT << " have been read in and will be used.\n";
               break; // Get out of the loop! 
            }
        }
    }
    fclose(fd);
    matSeries.nmat = nmat;

    maxRes = n;
    return (maxRes);
}

/*
 * This function is used to read a single user matrix from a file.
 * It can be called repeatedly if there are multiple matrices in the file.
 */
int SubMatrix::readUserMatrix(const char *fileName, Matrix& userMat, Xref& xref)
{
    double f;
    FILE *fd;
    int numargs, farg;
    int i, j, k = 0;
    char codes[NUMRES];
    char inline1[1024];
    char *args[NUMRES + 4];
    char c1, c2;
    int ix1, ix = 0;
    int maxRes = 0;
    float scale;

    if (fileName[0] == '\0')
    {
        utilityObject->error("comparison matrix not specified");
        return ((int)0);
    }

    if ((fd = fopen(fileName, "r")) == NULL)
    {
        utilityObject->error("cannot open %s", fileName);
        return ((int)0);
    }
    maxRes = 0;
    while (fgets(inline1, 1024, fd) != NULL)
    {
        if (commentline(inline1))
        {
            continue;
        }
        if (utilityObject->lineType(inline1, "CLUSTAL_SERIES"))
        {
            utilityObject->error("in %s - single matrix expected.", fileName);
            fclose(fd);
            return ((int)0);
        }
        /*
        read residue characters.
         */
        k = 0;
        for (j = 0; j < (int)strlen(inline1); j++)
        {
            if (isalpha((int)inline1[j]))
            {
                codes[k++] = inline1[j];
            }
            if (k > NUMRES)
            {
                utilityObject->error("too many entries in matrix %s", fileName);
                fclose(fd);
                return ((int)0);
            }
        }
        codes[k] = '\0';
        break;
    }

    if (k == 0)
    {
        utilityObject->error("wrong format in matrix %s", fileName);
        fclose(fd);
        return ((int)0);
    }

    /*
    cross-reference the residues
     */
    for (i = 0; i < NUMRES; i++)
    {
        xref[i] =  - 1;
    }

    maxRes = 0;
    for (i = 0; (c1 = codes[i]); i++)
    {
        for (j = 0; (c2 = userParameters->getAminoAcidCode(j)); j++)
        if (c1 == c2)
        {
            xref[i] = j;
            maxRes++;
            break;
        }
        if ((xref[i] ==  - 1) && (codes[i] != '*'))
        {
            utilityObject->warning("residue %c in matrix %s not recognised", codes[i],
                fileName);
        }
    }


    /*
    get the weights
     */

    ix = ix1 = 0;
    while (fgets(inline1, 1024, fd) != NULL)
    {
        if (inline1[0] == '\n')
        {
            continue;
        }
        if (inline1[0] == '#' || inline1[0] == '!')
        {
            break;
        }
        numargs = getArgs(inline1, args, (int)(k + 1));
        if (numargs < maxRes)
        {
            utilityObject->error("wrong format in matrix %s", fileName);
            fclose(fd);
            return ((int)0);
        }
        if (isalpha(args[0][0]))
        {
            farg = 1;
        }
        else
        {
            farg = 0;
        }

        /* decide whether the matrix values are float or decimal */
        scale = 1.0;
        for (i = 0; i < (int)strlen(args[farg]); i++)
        if (args[farg][i] == '.')
        {
            /* we've found a float value */
            scale = 10.0;
            break;
        }

        for (i = 0; i <= ix; i++)
        {
            if (xref[i] !=  - 1)
            {
                f = atof(args[i + farg]);
                userMat[ix1++] = (short)(f *scale);
            }
        }
        ix++;
    }
    if (ix != k + 1)
    {
        utilityObject->error("wrong format in matrix %s", fileName);
        fclose(fd);
        return ((int)0);
    }

    userMat.resize(ix1 + 1);

    maxRes += 2;
    fclose(fd);

    return (maxRes);
}

bool SubMatrix::getUserMatSeriesFromR(Rcpp::NumericMatrix substitutionMatrix) {
	int maxRes;

    maxRes = readMatrixSeriesFromR(substitutionMatrix, userMat, AAXref);
    if (maxRes <= 0) {
    	return false;
    }

    return true;
}

int SubMatrix::readMatrixSeriesFromR(Rcpp::NumericMatrix substitutionMatrix, Matrix& userMat, Xref& xref) {
    char inline1[1024];
    int maxRes = 0;
    int nmat;
    int n, llimit, ulimit;

    maxRes = readUserMatrixFromR(substitutionMatrix, userMat, xref);
    return (maxRes);
}

int SubMatrix::readUserMatrixFromR(Rcpp::NumericMatrix substitutionMatrix, Matrix& userMat, Xref& xref) {
	vector<string> matrixRows;
	stringstream line;

	Rcpp::Function rownames("rownames");
	vector<string> names = Rcpp::as<vector<string> >(rownames(substitutionMatrix));
	for (int i = 0, n = names.size(); i < n; i++) {
		if (i > 0) {
			line << " ";
		}
		line << names[i];
	}
	matrixRows.push_back(line.str());
	line.str("");

	int nrows = substitutionMatrix.nrow();
	int ncols = substitutionMatrix.ncol();
	for (int i = 0; i < ncols; i++) {
		//line << names[i];
		for (int j = 0; j < nrows; j++) {
			if (j > 0) {
				line << " ";
			}
			line << substitutionMatrix(j,i);
		}
		matrixRows.push_back(line.str());
		line.str("");
	}

    double f;
    // FILE *fd;
    int numargs, farg;
    int i, j, k = 0;
    char codes[NUMRES];
    string curline;
    char *args[NUMRES + 4];
    char c1, c2;
    int ix1, ix = 0;
    int maxRes = 0;
    float scale;

	/*
	read residue characters.
	 */
    for (int i = 0, n = matrixRows.size(); i < n; i++) {
    	curline = matrixRows[i];
		k = 0;
		for (j = 0; j < curline.length(); j++) {
			if (isalpha((int)curline.at(j))) {
				codes[k++] = curline.at(j);
			}
			if (k > NUMRES) {
				utilityObject->error("too many entries in matrix");
				return ((int)0);
			}
		}
		codes[k] = '\0';
		break;
    }

	if (k == 0) {
        utilityObject->error("wrong format in matrix");
        return ((int)0);
    }

    /*
    cross-reference the residues
     */
    for (i = 0; i < NUMRES; i++) {
        xref[i] =  - 1;
    }

    maxRes = 0;
    for (i = 0; (c1 = codes[i]); i++) {
        for (j = 0; (c2 = userParameters->getAminoAcidCode(j)); j++)
        if (c1 == c2) {
            xref[i] = j;
            maxRes++;
            break;
        }
        if ((xref[i] ==  - 1) && (codes[i] != '*')) {
            utilityObject->warning("residue %c in matrix not recognised", codes[i]);
        }
    }

    /*
    get the weights
     */

    ix = ix1 = 0;
    for (int i = 0, n = matrixRows.size(); i < n; i++) {
    	curline = matrixRows[i];
        if (curline[0] == '\n') {
            continue;
        }
        if (curline[0] == '#' || curline[0] == '!') {
            break;
        }
        numargs = getArgs(curline, args, (int)(k + 1));
        if (numargs < maxRes) {
            utilityObject->error("wrong format in matrix");
            return ((int)0);
        }
        if (isalpha(args[0][0])) {
            farg = 1;
        }
        else {
            farg = 0;
        }

        /* decide whether the matrix values are float or decimal */
        scale = 1.0;
        for (i = 0; i < (int)strlen(args[farg]); i++)
        if (args[farg][i] == '.') {
            /* we've found a float value */
            scale = 10.0;
            break;
        }

        for (i = 0; i <= ix; i++) {
            if (xref[i] !=  - 1) {
                f = atof(args[i + farg]);
                userMat[ix1++] = (short)(f *scale);
            }
        }
        ix++;
    }

    if (ix != k + 1) {
        utilityObject->error("wrong format in matrix");
        return ((int)0);
    }

    userMat.resize(ix1 + 1);
    
    maxRes += 2;
    //    fclose(fd);

    return (maxRes);
}

/**
 * 
 * @param inline1 
 * @param args[] 
 * @param max 
 * @return 
 */
int SubMatrix::getArgs(string line, char *args[],int max)
{
    char *inptr;
    int i;

    inptr = strdup(line.c_str());
    for (i = 0; i <= max; i++)
    {
        if ((args[i] = strtok(inptr, " \t\n")) == NULL)
        {
            break;
        }
        inptr = NULL;
    }

    return (i);
}


/**
 * 
 * @return 
 */
int SubMatrix::getMatrixNum()
{
    return matrixNum;
}

/**
 * 
 * @return 
 */
int SubMatrix::getDNAMatrixNum()
{
    return DNAMatrixNum;
}

/**
 * 
 * @return 
 */
int SubMatrix::getPWMatrixNum()
{
    return pwMatrixNum;
}

/**
 * 
 * @return 
 */
int SubMatrix::getPWDNAMatrixNum()
{
    return pwDNAMatrixNum;
}

/**
 * The function setCurrentNameAndNum is used to select a matrix series.
 * This will then be used for the alignment. The matrices will change, but the
 * series remains the same. We can set the series for pairwise/full for both 
 * protein and DNA. NOTE: The default can be set to user defined matrices.
 * @param _matrixName 
 * @param _matrixNum 
 * @param alignResidueType 
 * @param alignType 
 */
void SubMatrix::setCurrentNameAndNum(string _matrixName, int _matrixNum, 
                                          int alignResidueType,int alignType)
{       
    // Check if the values are valid.
    checkResidueAndAlignType(alignResidueType, alignType);
    
    string residue;
    string align;
    if((alignResidueType == Protein) && (alignType == Pairwise))
    {
        residue = "Protein"; align = "Pairwise";
        pwMatrixNum = _matrixNum;
        pwMatrixName = new string(_matrixName);
    }
    else if((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        residue = "Protein"; align = "MultipleAlign";
        matrixNum = _matrixNum;
        matrixName = new string(_matrixName); 
    }
    else if((alignResidueType == DNA) && (alignType == Pairwise))
    {
        residue = "DNA"; align = "Pairwise";
        pwDNAMatrixNum = _matrixNum;
        pwDNAMatrixName = new string(_matrixName);     
    }
    else if((alignResidueType == DNA) && (alignType == MultipleAlign))
    {
        residue = "DNA"; align = "MultipleAlign";    
        DNAMatrixNum = _matrixNum;
        DNAMatrixName = new string(_matrixName);
            
    }
    
    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            ostringstream outs;
            outs << "The matrix/matrix series has been changed for "
                 << "(" << residue << " AND " << align << ")."
                 << " New value: " << _matrixName << "\n\n";           
            logObject->logMsg(outs.str());
        }
    #endif    
}

/**
 * 
 * @param alignResidueType 
 * @param alignType 
 * @return 
 */
int SubMatrix::getMatrixNumForMenu(int alignResidueType, int alignType)
{
    checkResidueAndAlignType(alignResidueType, alignType);
    
    if((alignResidueType == Protein) && (alignType == Pairwise))
    {
        return pwMatrixNum;
    }
    else if((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        return matrixNum;  
    }
    else if((alignResidueType == DNA) && (alignType == Pairwise))
    {
        return pwDNAMatrixNum;   
    }
    else if((alignResidueType == DNA) && (alignType == MultipleAlign))
    {
        return DNAMatrixNum;   
    }
    else
        return -100; // NOTE NONE of these. I need to put in better error checking
}

/**
 * 
 * @param line 
 * @return 
 */
bool SubMatrix::commentline(char* line)
{
    int i;

    if (line[0] == '#')
    {
        return true;
    }
    for (i = 0; line[i] != '\n' && line[i] != EOS; i++)
    {
        if (!isspace(line[i]))
        {
            return false;
        }
    }
    return true;
}

/**
 * This function prints out the vector to the file specified by name. 
 * The vector is printed out in a triangular format, the same as the way
 * the arrays are displayed in matrices.h. This function is for testing purposes.
 * @param temp 
 * @param name 
 */
void SubMatrix::printInFormat(vector<short>& temp, const char* name)
{
    char nameOfFile[30];
    strcpy(nameOfFile, name);
    strcat(nameOfFile, ".out");

    ofstream outfile(nameOfFile);
    
    if(!outfile)
        cerr<<"oops failed to open !!!\n";
    
    outfile<<"short "<<name<<"[]{\n";
            
    int numOnCurrentLine = 1;
    int soFar = 0;
    for(int i = 0; i < (int)temp.size(); i++)
    {
        if(soFar == numOnCurrentLine)
        {
            outfile<<"\n";
            soFar = 0;
            numOnCurrentLine++;            
        }
        if((temp[i] > 9) || (temp[i] < 0))
        {
            outfile <<" "<< temp[i]<<",";
        }
        else
        {
            outfile <<"  "<< temp[i]<<",";
        }
        
        // Now increment so far
        soFar++;
        // Check to see if the next element is the last element
        if((i + 1) == (int)temp.size() - 1)
        {
            // Print out the last element. Then a curly brace, and break.
            if((temp[i+1] > 9) || (temp[i+1] < 0))
            {
                outfile <<" "<< temp[i + 1]<<"};\n";
            }
            else
            {
                outfile <<"  "<< temp[i + 1]<<"};\n";
            }
            break;
        }  
    }
    
    ofstream outfile2("temp.out");
    for(int i = 0; i < (int)temp.size(); i++)
    {
        outfile2 << temp[i] << " ";
    }   
}


/**
 * 
 * @param temp 
 * @param name 
 */
void SubMatrix::printVectorToFile(vector<short>& temp, const char* name)
{
    char nameOfFile[30];
    strcpy(nameOfFile, name);
    strcat(nameOfFile, ".out");

    ofstream outfile(nameOfFile);
    
    if(!outfile)
        cerr<<"oops failed to open !!!\n";
    
    for(int i = 0; i < (int)temp.size(); i++)
    {
        if((temp[i] > 9) || (temp[i] < 0))
        {
            outfile <<" "<< temp[i]<<",";
        }
        else
        {
            outfile <<"  "<< temp[i]<<",";
        }
    }
    outfile.close();
}


/**
 * 
 * @param alignResidueType 
 * @param alignType 
 * @return 
 */
Matrix* SubMatrix::getUserMatAddress(int alignResidueType, int alignType)
{
    if((alignResidueType == Protein) && (alignType == Pairwise))
    {
        return &pwUserMat;
    }
    else if((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        return &userMat; 
    }
    else if((alignResidueType == DNA) && (alignType == Pairwise))
    {
        return &pwUserDNAMat; 
    }
    else if((alignResidueType == DNA) && (alignType == MultipleAlign))
    {
        return &userDNAMat;  
    }
    return NULL;
}


/**
 * 
 * @param alignResidueType 
 * @param alignType 
 * @return 
 */
Xref* SubMatrix::getUserXrefAddress(int alignResidueType, int alignType)
{
    if((alignResidueType == Protein) && (alignType == Pairwise))
    {
        return &pwAAXref;
    }
    else if((alignResidueType == Protein) && (alignType == MultipleAlign))
    {
        return &AAXref; 
    }
    else if((alignResidueType == DNA) && (alignType == Pairwise))
    {
        return &pwDNAXref;  
    }
    else if((alignResidueType == DNA) && (alignType == MultipleAlign))
    {
        return &DNAXref;   
    }
    return NULL;
}

/**
 * This is an error handling routine. If an incorrect combination of values
 * is given, it will terminate the program.
 * @param alignResidueType 
 * @param alignType 
 */
void SubMatrix::checkResidueAndAlignType(int alignResidueType, int alignType)
{
    if(((alignResidueType != 0) && (alignResidueType != 1))
           || ((alignType != 0) && (alignType != 1)))
    {
        InvalidCombination ex(alignResidueType, alignType);
        ex.whatHappened();
        throw 1;
    }
}

/**
 * The function tempInterface is used to call the SubMatrix in the way it is
 * supposed to be used. It is for testing purposes.
 * @param alignResidueType 
 * @param alignType 
 */
void SubMatrix::tempInterface(int alignResidueType, int alignType)
{
    /*char userFile[FILENAMELEN + 1];

    userParameters->setDNAFlag(true);
    strcpy(userFile, "mat1");
    userParameters->setMenuFlag(false);   
    getUserMatFromFile(userFile, DNA, Pairwise);
    setCurrentNameAndNum(userFile, 4, 3, Pairwise);
    
    setCurrentNameAndNum("gonnet", 4, Protein, Pairwise);
    */
}

/**
 * A single matrix is used in scoring the alignment. This is Blosum45. This is the 
 * function to get it.
 * @param matrix[][] 
 * @return 
 */
int SubMatrix::getAlnScoreMatrix(int matrix[NUMRES][NUMRES])
{
    int _maxNumRes;
    /* 
       //_maxNumRes = getMatrix(blosum45mtVec, &defaultAAXref, matrix, true, 100);
       _maxNumRes = getMatrix(blosum45mtVec, &defaultAAXref, matrix, false, 1, true);
       //_maxNumRes = getMatrix(blosum62mt2Vec, &defaultAAXref, matrix, true, 100);
    */
    // 1.83 style
    _maxNumRes = getMatrix(blosum45mtVec, &defaultAAXref, matrix, true, 100);

    return _maxNumRes;
}

/**
 * This function is used to get the matrix that will be used for the calculation of 
 * the histogram. The histogram values are used in ClustalQt. 
 * @param matrix[][] 
 * @param matNum 
 * @param dnaMatNum 
 */
void SubMatrix::getQTMatrixForHistogram(int matrix[NUMRES][NUMRES])
{
    Matrix* _matPtrLocal;
    Xref* _matXrefLocal;
    int maxRes;
    if(userParameters->getDNAFlag())
    {
        if (QTDNAHistMatNum == DNAUSERDEFINED)
        {
            _matPtrLocal = &QTscoreUserDNAMatrix;
            _matXrefLocal = &QTscoreDNAXref;
        }
        else if (QTDNAHistMatNum == DNACLUSTALW)
        {
            _matPtrLocal = clustalvdnamtVec;
            _matXrefLocal = &defaultDNAXref;
        }
        else
        {
            _matPtrLocal = swgapdnamtVec;
            _matXrefLocal = &defaultDNAXref;            
        }
    }
    else 
    {
        if (QTAAHistMatNum == AAHISTIDENTITY)
        {
            _matPtrLocal = idmatVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTAAHistMatNum == AAHISTGONNETPAM80)
        {
            _matPtrLocal = gon80mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTAAHistMatNum == AAHISTGONNETPAM120)
        {
            _matPtrLocal = gon120mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTAAHistMatNum == AAHISTUSER)
        {
            _matPtrLocal = &QTscoreUserMatrix;
            _matXrefLocal = &QTscoreXref;
        }
        else if (QTAAHistMatNum == AAHISTGONNETPAM350)
        {
            _matPtrLocal = gon350mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else // Default
        {
            _matPtrLocal = gon250mtVec;
            _matXrefLocal = &defaultAAXref;            
        }
    }
    maxRes = getMatrix(_matPtrLocal, _matXrefLocal, matrix, false, 100);

}

void SubMatrix::getQTMatrixForLowScoreSeg(int matrix[NUMRES][NUMRES])
{
    Matrix* _matPtrLocal;
    Xref* _matXrefLocal;
    int maxRes;
    int _maxAA = userParameters->getMaxAA();
    int max = 0;
    int offset;
        
    if(userParameters->getDNAFlag())
    {
        if (QTsegmentDNAMatNum == DNAUSERDEFINED)
        {
            _matPtrLocal = &QTsegmentDNAMatrix;
            _matXrefLocal = &QTsegmentDNAXref;            
        }
        else if (QTsegmentDNAMatNum == DNACLUSTALW)
        {
            _matPtrLocal = clustalvdnamtVec;
            _matXrefLocal = &defaultDNAXref;
        }
        else
        {
            _matPtrLocal = swgapdnamtVec;
            _matXrefLocal = &defaultDNAXref;
        }
        /* get a positive matrix - then adjust it according to scale */
        maxRes = getMatrix(_matPtrLocal, _matXrefLocal, matrix, false, 100);
        /* find the maximum value */
        for(int i = 0; i <= _maxAA; i++)
        {
            for(int j = 0; j <= _maxAA; j++)
            {
                if(matrix[i][j] > max)
                { 
                    max = matrix[i][j];
                }
            }
        }
        /* subtract max * scale / 2 from each matrix value */
        offset = static_cast<int>(static_cast<float>(max *
                                  userParameters->getQTlowScoreDNAMarkingScale()) / 20.0);

        for(int i = 0; i <= _maxAA; i++)
        {
            for(int j = 0; j <= _maxAA; j++)
            {
                matrix[i][j] -= offset;
            }
        }
    }
    else
    {
        if (QTsegmentAAMatNum == QTAASEGGONNETPAM80)
        {
            _matPtrLocal = gon80mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTsegmentAAMatNum == QTAASEGGONNETPAM120)
        {
            _matPtrLocal = gon120mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else if (QTsegmentAAMatNum == QTAASEGUSER)
        {
            _matPtrLocal = &QTsegmentAAMatrix;
            _matXrefLocal = &QTsegmentAAXref;
        }
        else if (QTsegmentAAMatNum == QTAASEGGONNETPAM350)
        {
            _matPtrLocal = gon350mtVec;
            _matXrefLocal = &defaultAAXref;
        }
        else
        {
            _matPtrLocal = gon250mtVec;
            _matXrefLocal = &defaultAAXref;            
        }
        /* get a negative matrix */
        maxRes = getMatrix(_matPtrLocal, _matXrefLocal, matrix, true, 100);
    }
}

bool SubMatrix::getQTLowScoreMatFromFile(char* fileName, bool dna)
{
    int maxRes;

    FILE *infile;

    line2 = string(fileName);

    if(line2.size() == 0)
    {
        return false;
    }

    if((infile = fopen(line2.c_str(), "r")) == NULL)
    {
        utilityObject->error("Cannot find matrix file [%s]", line2.c_str());
        return false;
    }

    strcpy(fileName, line2.c_str());

    if(dna)
    {
        maxRes = readUserMatrix(fileName, QTsegmentDNAMatrix, QTsegmentDNAXref);
    }
    else
    {
        maxRes = readUserMatrix(fileName, QTsegmentAAMatrix, QTsegmentAAXref);
    }

    if (maxRes <= 0)
    {
        return false;
    }

    return true;
}

bool SubMatrix::getAAScoreMatFromFile(char *str)
{
    int maxRes;

    FILE *infile;

    line2 = string(str);

    if(line2.size() == 0)
    {
        return false;
    }

    if((infile = fopen(line2.c_str(), "r")) == NULL)
    {
        utilityObject->error("Cannot find matrix file [%s]", line2.c_str());
        return false;
    }

    strcpy(str, line2.c_str());

    maxRes = readUserMatrix(str, QTscoreUserMatrix, QTscoreXref);

    if (maxRes <= 0)
    {
        return false;
    }

    return true;
}

bool SubMatrix::getDNAScoreMatFromFile(char *str)
{
    int maxRes;

    FILE *infile;

    line2 = string(str);

    if(line2.size() == 0)
    {
        return false;
    }

    if((infile = fopen(line2.c_str(), "r")) == NULL)
    {
        utilityObject->error("Cannot find matrix file [%s]", line2.c_str());
        return false;
    }

    strcpy(str, line2.c_str());

    maxRes = readUserMatrix(str, QTscoreUserDNAMatrix, QTscoreDNAXref);

    if (maxRes <= 0)
    {
        return false;
    }

    return true;
}


bool SubMatrix::getQTLowScoreMatFromR(Rcpp::NumericMatrix substitutionMatrix, bool dna) {
	printf("getQTLowScoreMatFromFile\n");

	int maxRes;
    
    if(dna) {
        maxRes = readUserMatrixFromR(substitutionMatrix, QTsegmentDNAMatrix, QTsegmentDNAXref);
    }
    else {
        maxRes = readUserMatrixFromR(substitutionMatrix, QTsegmentAAMatrix, QTsegmentAAXref);
    }
    
    if (maxRes <= 0) {
        return false;
    }

    return true;
}

bool SubMatrix::getAAScoreMatFromR(Rcpp::NumericMatrix substitutionMatrix) {
	printf("getAAScoreMatFromFile\n");

	int maxRes = readUserMatrixFromR(substitutionMatrix, QTscoreUserMatrix, QTscoreXref);
    
    if (maxRes <= 0) {
        return false;
    }
    return true;
}

bool SubMatrix::getDNAScoreMatFromR(Rcpp::NumericMatrix substitutionMatrix) {
	printf("getDNAScoreMatFromFile\n");

	int maxRes = readUserMatrixFromR(substitutionMatrix, QTscoreUserDNAMatrix, QTscoreDNAXref);
    
    if (maxRes <= 0) {
        return false;
    }
    return true;
}

}
