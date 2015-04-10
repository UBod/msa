/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * @author Mark Larkin, Conway Institute, UCD. mark.larkin@ucd.ie
 * Changes:
 *
 * Mark: 23-01-2007: There was a problem with running an alignment with only 
 * one sequence in it. I needed to make a change to the multiSeqAlign function.
 ****************************************************************************/
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "MSA.h"
#include "MyersMillerProfileAlign.h"
//#include "../phylogeneticTree/ClusterTree.h"
#include "../general/debuglogObject.h"

namespace clustalw
{

/**
 * 
 * @param alnPtr 
 * @param distMat 
 * @param iStart 
 * @param phylipName 
 * @return 
 */
int MSA::multiSeqAlign(Alignment* alnPtr, DistMatrix* distMat, vector<int>* seqWeight, AlignmentSteps* progSteps, int iStart)
{
        
    if(!progSteps)
    {
        return 0;
    }
        
    int* aligned;
    vector<int> group;
    int ix;
    
    int* maxid;
    int max = 0, sum = 0;
    vector<int> treeWeight;
    int i = 0, j = 0, set = 0, iseq = 0;
    int entries = 0;
    int score = 0;
    int _numSteps = 0;
    utilityObject->info("Start of Multiple Alignment\n"); 

    int _numSeqs = alnPtr->getNumSeqs();
    
    vector<int> newOutputIndex(_numSeqs);   
   
    alnPtr->addSeqWeight(seqWeight);
    
    ProfileAlignAlgorithm* alignAlgorithm = new MyersMillerProfileAlign;
    _numSteps = progSteps->getNumSteps();
    // for each sequence, find the most closely related sequence

    maxid = new int[_numSeqs + 1];

    for (i = 1; i <= _numSeqs; i++)
    {
        maxid[i] =  -1;
        for (j = 1; j <= _numSeqs; j++)
            if (j != i && maxid[i] < (*distMat)(i, j))
            {
                maxid[i] = static_cast<int>((*distMat)(i, j));
            }
    }

    // group the sequences according to their relative divergence

    if (iStart == 0)
    {

        // start the multiple alignments.........

        utilityObject->info("Aligning...");
        // first pass, align closely related sequences first....
        
        
        ix = 0;
        aligned = new int[_numSeqs + 1];
        
        for (i = 0; i <= _numSeqs; i++)
        {
            aligned[i] = 0;
        }
        
        const vector<vector<int> >* ptrToSets = progSteps->getSteps();
        
        
        for (set = 1; set <= _numSteps; ++set)
        {
            entries = 0;
            for (i = 1; i <= _numSeqs; i++)
            {
                if (((*ptrToSets)[set][i] != 0) && 
                    (maxid[i] > userParameters->getDivergenceCutoff()))
                {
                    entries++;
                    if (aligned[i] == 0)
                    {
                        if (userParameters->getOutputOrder() == INPUT)
                        {
                            ++ix;
                            newOutputIndex[i - 1] = i;
                        }
                        else
                        {
                            if(ix >= (int)newOutputIndex.size())
                            {
                                cerr << "ERROR: size = " << newOutputIndex.size() 
                                     << "ix = " << ix << "\n";
                                throw 1;
                            }
                            else
                            {
                                newOutputIndex[ix] = i;
                                ++ix;
                            }
                        }

                        aligned[i] = 1;
                    }
                }
            }

            if (entries > 0)
            {
                #if DEBUGFULL
                    if(logObject && DEBUGLOG)
                    {    
                        logObject->logMsg("Doing profile align");
                    }
                #endif                 
                score = alignAlgorithm->profileAlign(alnPtr, distMat, progSteps->getStep(set),
                                                     aligned);
            }
            else
            {
                score = 0;
            }


            // negative score means fatal error... exit now! 

            if (score < 0)
            {
                return (-1);
            }
            if(userParameters->getDisplayInfo())
            {
                if ((entries > 0) && (score > 0))
                {
                    utilityObject->info("Group %d: Sequences:%4d      Score:%d",
                                        set, entries, score);                     
                }
                else
                {
                    utilityObject->info("Group %d:                     Delayed", set);
                }
            }
        }
    }

    else
    {
        aligned = new int[_numSeqs + 1];
        ix = 0;
        for (i = 1; i <= iStart + 1; i++)
        {
            aligned[i] = 1;
            ++ix;
            newOutputIndex[i - 1] = i;
        }
        for (i = iStart + 2; i <= _numSeqs; i++)
        {
            aligned[i] = 0;
        }
    }

    // second pass - align remaining, more divergent sequences..... 

    // if not all sequences were aligned, for each unaligned sequence,
    // find it's closest pair amongst the aligned sequences.

    group.resize(_numSeqs + 1); 
    treeWeight.resize(_numSeqs); 

    for (i = 0; i < _numSeqs; i++)
    {
        treeWeight[i] = (*seqWeight)[i];
    }

    // if we haven't aligned any sequences, in the first pass - align the
    // two most closely related sequences now
    if (ix == 0)
    {
        max =  -1;
        iseq = 0;

        for (i = 1; i <= _numSeqs; i++)
        {
            for (j = i + 1; j <= _numSeqs; j++)
            {
                if (max < (*distMat)(i, j))
                {
                    max = static_cast<int>((*distMat)(i, j)); // Mark change 17-5-07
                    iseq = i;
                }
            }
        }
        aligned[iseq] = 1;
        if (userParameters->getOutputOrder() == INPUT)
        {
            ++ix;
            newOutputIndex[iseq - 1] = iseq;
        }
        else
        {
            newOutputIndex[ix] = iseq;
            ++ix;
        }
    }

    while (ix < _numSeqs)
    {
        for (i = 1; i <= _numSeqs; i++)
        {
            if (aligned[i] == 0)
            {
                maxid[i] =  - 1;
                for (j = 1; j <= _numSeqs; j++)
                    if ((maxid[i] < (*distMat)(i, j)) && (aligned[j] != 0))
                    {
                        maxid[i] = static_cast<int>((*distMat)(i, j));// Mark change 17-5-07
                    }

            }
        }
        // find the most closely related sequence to those already aligned

        max =  - 1;
        iseq = 0;
        for (i = 1; i <= _numSeqs; i++)
        {
            if ((aligned[i] == 0) && (maxid[i] > max))
            {
                max = maxid[i];
                iseq = i;
            }
        }


        // align this sequence to the existing alignment
        // weight sequences with percent identity with profile
        // OR...., multiply sequence weights from tree by percent identity with new sequence 
        if (userParameters->getNoWeights() == false)
        {
            for (j = 0; j < _numSeqs; j++)
                if (aligned[j + 1] != 0)
                {
                    (*seqWeight)[j] = static_cast<int>(treeWeight[j] * (*distMat)(j + 1, iseq));
                }

            // Normalise the weights, such that the sum of the weights = INT_SCALE_FACTOR

            sum = 0;
            for (j = 0; j < _numSeqs; j++)
            {
                if (aligned[j + 1] != 0)
                {
                    sum += (*seqWeight)[j];
                }
            }

            if (sum == 0)
            {
                for (j = 0; j < _numSeqs; j++)
                {
                    (*seqWeight)[j] = 1;
                }

                sum = j;
            }
            for (j = 0; j < _numSeqs; j++)
            {
                if (aligned[j + 1] != 0)
                {
                    (*seqWeight)[j] = ((*seqWeight)[j] * INT_SCALE_FACTOR) / sum;
                    if ((*seqWeight)[j] < 1)
                    {
                        (*seqWeight)[j] = 1;
                    }
                }
            }
        }

        entries = 0;
        for (j = 1; j <= _numSeqs; j++)
        {
            if (aligned[j] != 0)
            {
                group[j] = 1;
                entries++;
            }
            else if (iseq == j)
            {
                group[j] = 2;
                entries++;
            }
        }
        
        alnPtr->addSeqWeight(seqWeight);
        aligned[iseq] = 1;
        
        score = alignAlgorithm->profileAlign(alnPtr, distMat, &group, aligned);
         
        if (userParameters->getOutputOrder() == INPUT)
        {
            ++ix;
            newOutputIndex[iseq - 1] = iseq;
        }
        else
        {
            newOutputIndex[ix] = iseq;
            ++ix;
        }
    }
        
    alnPtr->addOutputIndex(&newOutputIndex);
    
    if(userParameters->getDisplayInfo())
    {
      int alignmentScore = alnPtr->alignScore(); // ?? check, FS, 2009-05-18
    }
    
    delete alignAlgorithm;
    delete [] aligned;
    delete [] maxid;
    return (_numSeqs);
}


/**
 * 
 * @param alnPtr 
 * @param distMat 
 * @param iStart 
 * @param phylipName 
 * @return 
 */
int MSA::seqsAlignToProfile(Alignment* alnPtr, DistMatrix* distMat, vector<int>* seqWeight, int iStart, 
                           string phylipName)
{
    int *aligned;  
    vector<int> treeWeight;
    vector<int> group;
    int ix;

    int *maxid;
    int max = 0;
    int i = 0, j = 0, iseq = 0;
    int sum = 0, entries = 0;
    int score = 0;
    int _numSeqs = alnPtr->getNumSeqs();
    
    utilityObject->info("Start of Multiple Alignment\n"); 
       
    ProfileAlignAlgorithm* alignAlgorithm = new MyersMillerProfileAlign;
    
    // calculate sequence weights according to branch lengths of the tree -
     // weights in global variable seq_weight normalised to sum to 100
    vector<int> newOutputIndex(_numSeqs);    
    //groupTree.calcSeqWeights(0, _numSeqs, &seqWeight);

    treeWeight.resize(_numSeqs);
    for (i = 0; i < _numSeqs; i++)
    {
        treeWeight[i] = (*seqWeight)[i];
    }

    // for each sequence, find the most closely related sequence

    maxid = new int[_numSeqs + 1];
    
    for (i = 1; i <= _numSeqs; i++)
    {
        maxid[i] =  - 1;
        for (j = 1; j <= _numSeqs; j++)
        {
            if (maxid[i] < (*distMat)(i, j))
            {
                maxid[i] = static_cast<int>((*distMat)(i, j)); // Mark change 17-5-07
            }
        }
    }

    aligned = new int[_numSeqs + 1];
    ix = 0;
    for (i = 1; i <= iStart + 1; i++)
    {
        aligned[i] = 1;
        ++ix;
        newOutputIndex[i - 1] = i;
    }
    for (i = iStart + 2; i <= _numSeqs; i++)
    {
        aligned[i] = 0;
    }

    // for each unaligned sequence, find it's closest pair amongst the
    // aligned sequences. 

    group.resize(_numSeqs + 1);

    while (ix < _numSeqs)
    {
        if (ix > 0)
        {
            for (i = 1; i <= _numSeqs; i++)
            {
                if (aligned[i] == 0)
                {
                    maxid[i] =  - 1;
                    for (j = 1; j <= _numSeqs; j++)
                    {
                        if ((maxid[i] < (*distMat)(i, j)) && (aligned[j] != 0))
                        {
                            maxid[i] = static_cast<int>((*distMat)(i, j));
                        }
                    }
                }
            }
        }

        // find the most closely related sequence to those already aligned

        max =  -1;
        for (i = 1; i <= _numSeqs; i++)
        {
            if ((aligned[i] == 0) && (maxid[i] > max))
            {
                max = maxid[i];
                iseq = i;
            }
        }

        // align this sequence to the existing alignment 
        entries = 0;
        for (j = 1; j <= _numSeqs; j++)
        {
            if (aligned[j] != 0)
            {
                group[j] = 1;
                entries++;
            }
            else if (iseq == j)
            {
                group[j] = 2;
                entries++;
            }
        }

        aligned[iseq] = 1;

        // multiply sequence weights from tree by percent
        // identity with new sequence 

        for (j = 0; j < _numSeqs; j++)
        {
            (*seqWeight)[j] = static_cast<int>(treeWeight[j] * (*distMat)(j + 1, iseq));
        }

        //
        // Normalise the weights, such that the sum of the weights = INT_SCALE_FACTOR
        //

        sum = 0;
        for (j = 0; j < _numSeqs; j++)
        {
            if (group[j + 1] == 1)
            {
                sum += (*seqWeight)[j];
            }
        }

        if (sum == 0)
        {
            for (j = 0; j < _numSeqs; j++)
            {
                (*seqWeight)[j] = 1;
            }

            sum = j;
        }
        for (j = 0; j < _numSeqs; j++)
        {
            (*seqWeight)[j] = ((*seqWeight)[j] * INT_SCALE_FACTOR) / sum;
            if ((*seqWeight)[j] < 1)
            {
                (*seqWeight)[j] = 1;
            }
        }
        // Add the seqWeights to the Alignment object!!!!
        alnPtr->addSeqWeight(seqWeight);
        
        score = alignAlgorithm->profileAlign(alnPtr, distMat, &group, aligned);
        utilityObject->info("Sequence:%d     Score:%d", iseq, score);
        if (userParameters->getOutputOrder() == INPUT)
        {
            ++ix;
            newOutputIndex[iseq - 1] = iseq;
        }
        else
        {
            newOutputIndex[ix] = iseq;
            ++ix;
        }
    }

    delete [] aligned;
    delete [] maxid;
    delete alignAlgorithm;
    alnPtr->addOutputIndex(&newOutputIndex);
    
    if(userParameters->getDisplayInfo())
    {
        alnPtr->alignScore();
    }

    return (_numSeqs);
}


/**
 * 
 * @param alnPtr 
 * @param distMat 
 * @return 
 */
int MSA::calcPairwiseForProfileAlign(Alignment* alnPtr, DistMatrix* distMat)
{
    //Tree groupTree;
    int i, j, temp;
    int entries;
    int* aligned;  
    vector<int> group;
    vector<int> seqWeight;
    float dscore;
    int score;
    int _numSeqs = alnPtr->getNumSeqs();
    
    seqWeight.resize(_numSeqs);
    ProfileAlignAlgorithm* alignAlg = new MyersMillerProfileAlign;
    
    utilityObject->info("Start of Initial Alignment");
    /* calculate sequence weights according to branch lengths of the tree -
     * weights in global variable seq_weight normalised to sum to INT_SCALE_FACTOR */

    temp = INT_SCALE_FACTOR / _numSeqs;
    for (i = 0; i < _numSeqs; i++)
    {
        seqWeight[i] = temp;
    }

    userParameters->setDistanceTree(false);

    /* do the initial alignment.........  */

    group.resize(_numSeqs + 1);

    for (i = 1; i <= alnPtr->getProfile1NumSeqs(); ++i)
    {
        group[i] = 1;
    }

    for (i = alnPtr->getProfile1NumSeqs() + 1; i <= _numSeqs; ++i)
    {
        group[i] = 2;
    }

    entries = _numSeqs;

    aligned = new int[_numSeqs + 1];

    for (i = 1; i <= _numSeqs; i++)
    {
        aligned[i] = 1;
    }

    alnPtr->addSeqWeight(&seqWeight);
    
    score = alignAlg->profileAlign(alnPtr, distMat, &group, aligned);
    utilityObject->info("Sequences:%d      Score:%d", entries, score);
    delete [] aligned;
    
    for (i = 1; i <= _numSeqs; i++)
    {
        for (j = i + 1; j <= _numSeqs; j++)
        {
            dscore = alnPtr->countid(i, j);
            (*distMat)(i, j) = ((double)100.0 - (double)dscore) / (double)100.0;
            (*distMat)(j, i) = (*distMat)(i, j);
        }
    }
    delete alignAlg;
    return (_numSeqs);

}


/**
 * 
 * @param alnPtr 
 * @param distMat 
 * @param p1TreeName 
 * @param p2TreeName 
 * @return 
 */
int MSA::doProfileAlign(Alignment* alnPtr, DistMatrix* distMat, vector<int>* prof1Weight, vector<int>* prof2Weight)
{
    //Tree groupTree1, groupTree2;
    int i, j, sum, entries;
    int score;
    int *aligned;
    vector<int> group;
    int *maxid;
    //vector<int> prof1Weight, prof2Weight;
    int _profile1NumSeqs = alnPtr->getProfile1NumSeqs();
    int _numSeqs = alnPtr->getNumSeqs();
    vector<int> _seqWeight, _outputIndex;    
    
    utilityObject->info("Start of Multiple Alignment\n");
        
    _seqWeight.resize(_numSeqs + 1);
    _outputIndex.resize(_numSeqs);
    ProfileAlignAlgorithm* alignAlgorithm = new MyersMillerProfileAlign;

    // weight sequences with max percent identity with other profile

    maxid = new int[_numSeqs + 1];
    for (i = 0; i < _profile1NumSeqs; i++)
    {
        maxid[i] = 0;
        for (j = _profile1NumSeqs + 1; j <= _numSeqs; j++)
        {
            if (maxid[i] < (*distMat)(i + 1, j))
            {
                maxid[i] = static_cast<int>((*distMat)(i + 1, j)); // Mark change 17-5-07
            } 
        }
        _seqWeight[i] = maxid[i] * (*prof1Weight)[i];
    }

    for (i = _profile1NumSeqs; i < _numSeqs; i++)
    {
        maxid[i] =  - 1;
        for (j = 1; j <= _profile1NumSeqs; j++)
        {
            if (maxid[i] < (*distMat)(i + 1, j))
            {
                maxid[i] = static_cast<int>((*distMat)(i + 1, j));// Mark change 17-5-07
            }
        }
        _seqWeight[i] = maxid[i] * (*prof2Weight)[i];
    }
    //
    // Normalise the weights, such that the sum of the weights = INT_SCALE_FACTOR
    //

    sum = 0;
    for (j = 0; j < _numSeqs; j++)
    {
        sum += _seqWeight[j];
    }

    if (sum == 0)
    {
        for (j = 0; j < _numSeqs; j++)
        {
            _seqWeight[j] = 1;
        }

        sum = j;
    }
    for (j = 0; j < _numSeqs; j++)
    {
        _seqWeight[j] = (_seqWeight[j] * INT_SCALE_FACTOR) / sum;
        if (_seqWeight[j] < 1)
        {
            _seqWeight[j] = 1;
        }
    }

    // do the alignment.........  /

    utilityObject->info("Aligning...");
    group.resize(_numSeqs + 1);

    for (i = 1; i <= _profile1NumSeqs; ++i)
    {
        group[i] = 1;
    }

    for (i = _profile1NumSeqs + 1; i <= _numSeqs; ++i)
    {
        group[i] = 2;
    }

    entries = _numSeqs;

    aligned = new int[_numSeqs + 1];
    for (i = 1; i <= _numSeqs; i++)
    {
        aligned[i] = 1;
    }

    alnPtr->addSeqWeight(&_seqWeight);

    score = alignAlgorithm->profileAlign(alnPtr, distMat, &group, aligned);
   
    utilityObject->info("Sequences:%d      Score:%d", entries, score);
    
    for (i = 1; i <= _numSeqs; i++)
    {
        _outputIndex[i - 1] = i;
    }

    alnPtr->addOutputIndex(&_outputIndex);
    
    delete alignAlgorithm;
    delete [] aligned;
    delete [] maxid;
    return (_numSeqs);
    return 1;
}
 
}
