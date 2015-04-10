/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Changes:
 * Mark 10-5-2007: Bug fix # 42. Added getWeightsForQtLowScore function.
 * Mark 22-5-2007: Made a change to getWeightsForQtLowScore
 * Mark 23-5-2007: Made a change to getWeightsForQtLowScore
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "TreeInterface.h"
#include <sstream>
#include "UnRootedClusterTree.h"
#include "Tree.h"
#include "../multipleAlign/MSA.h"
#include "../general/userparams.h"
#include "../general/debuglogObject.h"
#include "UPGMA/RootedGuideTree.h"
#include "UPGMA/RootedClusterTree.h"
namespace clustalw
{

/**
 * This will be called by align!
 *
 */
auto_ptr<AlignmentSteps>
TreeInterface::getWeightsAndStepsFromDistMat(vector<int>* seqWeights, DistMatrix* distMat, 
                                             Alignment *alignPtr, int seq1, int nSeqs, 
                                             string* phylipName, bool* success)
{
    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("In getWeightsAndStepsFromDistMat\n");
        }
    #endif    
    if(userParameters->getClusterAlgorithm() == UPGMA)
    {
        return getWeightsAndStepsFromDistMatUPGMA(seqWeights, distMat, alignPtr, 
                                                seq1, nSeqs, phylipName, success);
    }
    else
    {
        return getWeightsAndStepsFromDistMatNJ(seqWeights, distMat, alignPtr, 
                                                seq1, nSeqs, phylipName, success);
    }
}

                                   
/**
 * This will be called by sequencesAlignToProfile. This function will put the distMat into
 * a similarity matrix.
 */
void
TreeInterface::getWeightsFromDistMat(vector<int>* seqWeights, DistMatrix* distMat, 
                                     Alignment *alignPtr, int seq1, int nSeqs, 
                                     string* phylipName, bool* success)
{
    if(userParameters->getClusterAlgorithm() == UPGMA)
    {    
        getWeightsFromDistMatUPGMA(seqWeights, distMat, alignPtr, seq1, nSeqs, phylipName,
                                   success);
    }
    else
    {
        getWeightsFromDistMatNJ(seqWeights, distMat, alignPtr, seq1, nSeqs, phylipName,
                                success);    
    }
}

/**
 * This function will be called by profileAlign
 */
void
TreeInterface::getWeightsForProfileAlign(Alignment* alignPtr, DistMatrix* distMat, 
                                         string* p1TreeName, vector<int>* p1Weights, string* p2TreeName, 
                                         vector<int>* p2Weights, int numSeqs, int profile1NumSeqs, bool useTree1, 
                                         bool useTree2, bool* success)
{    
    if(userParameters->getClusterAlgorithm() == UPGMA)
    {    
        getWeightsForProfileAlignUPGMA(alignPtr, distMat, p1TreeName, p1Weights, p2TreeName,
                                p2Weights, numSeqs, profile1NumSeqs, useTree1, useTree2,
                                success);
    }
    else
    {
        getWeightsForProfileAlignNJ(alignPtr, distMat, p1TreeName, p1Weights, p2TreeName,
                                p2Weights, numSeqs, profile1NumSeqs, useTree1, useTree2,
                                success);    
    }                               
}                    

/**
 * This function will be called by doAlignUseOldGuideTree
 *
 */                                   
auto_ptr<AlignmentSteps>
TreeInterface::getWeightsAndStepsFromTree(Alignment* alignPtr, 
                                          DistMatrix* distMat, string* treeName,
                                          vector<int>* seqWeights, int fSeq, 
                                          int numSeqs, bool* success)         
{
    /**
     * This will only use the NJ. It will not use UPGMA
     */
    return getWeightsAndStepsFromTreeNJ(alignPtr, distMat, treeName, seqWeights, 
                                        fSeq, numSeqs, success);
}

/**
 * Called by sequencesAlignToProfile
 */
int
TreeInterface::getWeightsFromGuideTree(Alignment* alignPtr, DistMatrix* distMat,
                                       string* treeName, vector<int>* seqWeights, int fSeq,
                                       int nSeqs, bool* success)
{
    /**
     * This will only use the NJ. It will not use UPGMA
     */    
    return getWeightsFromGuideTreeNJ(alignPtr, distMat, treeName, seqWeights, fSeq, nSeqs, 
                                      success);
}
 
/**
 * This will be called by doGuideTreeOnly. This does not put the distMat into similarity.
 *
 */                                     
void
TreeInterface::generateTreeFromDistMat(DistMatrix* distMat, Alignment *alignPtr, 
                                       int seq1, int nSeqs, 
                                       string* phylipName, bool* success)
{
    /**
     * This function does not put distMat into similarities
     */
    if(userParameters->getClusterAlgorithm() == UPGMA)
    { 
        RootedGuideTree guideTree;
        generateTreeFromDistMatUPGMA(&guideTree, distMat, alignPtr, seq1, nSeqs, 
                                     phylipName, success);
    }
    else
    {                                            
        generateTreeFromDistMatNJ(distMat, alignPtr, seq1, nSeqs, phylipName, success);
    }
}

void
TreeInterface::treeFromAlignment(TreeNames* treeNames, Alignment *alignPtr)
{
    if(userParameters->getClusterAlgorithm() == UPGMA)
    {
        RootedClusterTree clusterTree;
        clusterTree.treeFromAlignment(treeNames, alignPtr);
    }
    else
    {
        UnRootedClusterTree clusterTree;
        clusterTree.treeFromAlignment(treeNames, alignPtr);
    }
}                   
 
void
TreeInterface::bootstrapTree(TreeNames* treeNames, Alignment *alignPtr)
{
    /**
     * NOTE We only do bootstrapping on the NJ tree.
     */
    UnRootedClusterTree clusterTree;
    clusterTree.bootstrapTree(treeNames, alignPtr);
}


/**
 * Private functions!!!!
 *
 */

 
/** *******************
 *      
 *      Neighbour joining functions
 */ 
/**
 * Note: After this function has been called the distance matrix will all be in similarities.
 *
 */   
auto_ptr<AlignmentSteps>
TreeInterface::getWeightsAndStepsFromDistMatNJ(vector<int>* seqWeights, DistMatrix* distMat, 
                                               Alignment *alignPtr, int seq1, int nSeqs, 
                                               string* phylipName, bool* success)
{   
    auto_ptr<AlignmentSteps> progSteps;   
    generateTreeFromDistMatNJ(distMat, alignPtr, seq1, nSeqs, phylipName, success);
        
    progSteps = getWeightsAndStepsUseOldGuideTreeNJ(distMat, alignPtr, phylipName,
                                                    seqWeights, seq1, nSeqs, success);
    return progSteps;
}

auto_ptr<AlignmentSteps>
TreeInterface::getWeightsAndStepsUseOldGuideTreeNJ(DistMatrix* distMat, Alignment *alignPtr,
                                                   string* treeName, vector<int>* seqWeights, 
                                                   int fSeq, int nSeqs, bool* success)
{   
    Tree groupTree; 
    auto_ptr<AlignmentSteps> progSteps;
    
    if(nSeqs == 1)
    {
        utilityObject->info("Only 1 sequence, cannot do multiple alignment\n");
        *success = false;
        return progSteps;
    }
    
    int status = 0;
    
    status = readTreeAndCalcWeightsNJ(&groupTree, alignPtr, distMat, treeName,
                                          seqWeights, fSeq, nSeqs);
                                     
    if(status == 0)
    {
        *success = false;
        return progSteps;
    }
    
    progSteps = groupTree.createSets(0, nSeqs);
    int _numSteps = progSteps->getNumSteps();
    
    utilityObject->info("There are %d groups", _numSteps);        
    // clear the memory used for the phylogenetic tree

    if (nSeqs >= 2)
    {
        groupTree.clearTree(NULL);
    }
    *success = true;
    return progSteps;   
}
                                                                                              

                                                                                
/**
 * The function readTreeAndCalcWeightsNJ is used to read in the tree given the treeName.
 * It then calls the appropriate functions to calc the seqWeights and make sure the matrix 
 * is in similarity mode. 
 */                                                                                
int
TreeInterface::readTreeAndCalcWeightsNJ(Tree* groupTree, Alignment* alignPtr, 
                                        DistMatrix* distMat, string* treeName,
                                        vector<int>* seqWeights, int fSeq, int nSeqs)
{
    int status = 0;
    if (nSeqs >= 2)
    {
        status = groupTree->readTree(alignPtr, treeName->c_str(), fSeq - 1, nSeqs);
        if (status == 0)
        {
            return status;
        }
    }

    groupTree->calcSeqWeights(fSeq - 1, nSeqs, seqWeights);
    
    status = groupTree->calcSimilarities(alignPtr, distMat);
    
    return status;                                                            
}

int
TreeInterface::getWeightsFromGuideTreeNJ(Alignment* alignPtr, DistMatrix* distMat,
                                         string* treeName, vector<int>* seqWeights, int fSeq, int nSeqs, 
                                         bool* success)
{
    Tree groupTree;
    int status = readTreeAndCalcWeightsNJ(&groupTree, alignPtr, distMat, treeName,
                                          seqWeights, fSeq, nSeqs);
    if(status == 0)
    {
        *success = false;
    }
    else
    {
        *success = true;
    }
    return status;
}                                                                                             

void
TreeInterface::getWeightsFromDistMatNJ(vector<int>* seqWeights, DistMatrix* distMat, 
                                       Alignment *alignPtr, int seq1, int _nSeqs, 
                                       string* phylipName, bool* success)
{
    string copyOfPhylipName = string(*phylipName);     
    int nSeqs = _nSeqs;
    int status = 0;
    
    generateTreeFromDistMatNJ(distMat, alignPtr, seq1, nSeqs, phylipName, success);
    
    Tree groupTree;
    status = readTreeAndCalcWeightsNJ(&groupTree, alignPtr, distMat, phylipName,
                                          seqWeights, seq1, nSeqs);
    if(status == 0)
    {
        *success = false;
    }
    else
    {
        *success = true;
    }

}

void
TreeInterface::getWeightsForProfileAlignNJ(Alignment* alignPtr, DistMatrix* distMat, 
                                           string* p1TreeName, vector<int>* p1Weights, string* p2TreeName, 
                                           vector<int>* p2Weights, int numSeqs, int profile1NumSeqs, bool useTree1, 
                                           bool useTree2, bool* success)
{
    if(!useTree1)
    {
        if (profile1NumSeqs >= 2) 
        {        
            generateTreeFromDistMatNJ(distMat, alignPtr, 1, profile1NumSeqs, p1TreeName,
                                        success);
        }    
    }
    
    if(!useTree2)
    {
        if(numSeqs - profile1NumSeqs >= 2) 
        {        
            generateTreeFromDistMatNJ(distMat, alignPtr, profile1NumSeqs + 1, 
                                      numSeqs - profile1NumSeqs, p2TreeName, success);
        }  
    }
    
    if (userParameters->getNewTree1File() || userParameters->getNewTree2File()) 
    {
        *success = false;
        return;
    }    
    // MSA->CALCPAIRWISE
    
    MSA* msaObj = new MSA();
        
    int count = msaObj->calcPairwiseForProfileAlign(alignPtr, distMat);
    
    if (count == 0) 
    {
        *success = false;
        return;
    }
    
    Tree groupTree1, groupTree2;
    int status = 0;
        
    if (profile1NumSeqs >= 2)
    {
        status = groupTree1.readTree(alignPtr, p1TreeName->c_str(), 0, profile1NumSeqs);
        if (status == 0)
        {
            *success = false;
            return;
        }
    }
    
    groupTree1.calcSeqWeights(0, profile1NumSeqs, p1Weights);
    
    if (profile1NumSeqs >= 2)
    {
        groupTree1.clearTree(NULL);
    }
    
    if (numSeqs - profile1NumSeqs >= 2)
    {
        status = groupTree2.readTree(alignPtr, p2TreeName->c_str(), profile1NumSeqs, numSeqs);
        if (status == 0)
        {
            *success = false;
            return;
        }
    }
    
    groupTree2.calcSeqWeights(profile1NumSeqs, numSeqs, p2Weights);


    /* clear the memory for the phylogenetic tree */

    if (numSeqs - profile1NumSeqs >= 2)
    {
        groupTree2.clearTree(NULL);
    }
    
    /**
     * Convert distances to similarities!!!!!!
     */
    for (int i = 1; i < numSeqs; i++)
    {
        for (int j = i + 1; j <= numSeqs; j++)
        {
            (*distMat)(i, j) = 100.0 - (*distMat)(i, j) * 100.0;
            (*distMat)(j, i) = (*distMat)(i, j);
        }
    }    
      
    *success = true;
            
}

void
TreeInterface::generateTreeFromDistMatNJ(DistMatrix* distMat, Alignment *alignPtr, 
                                         int seq1, int nSeqs, 
                                         string* phylipName, bool* success)
{
    string copyOfPhylipName = string(*phylipName);     
    
    if (nSeqs >= 2) 
    {
        UnRootedClusterTree* clusterTree = new UnRootedClusterTree;   

        clusterTree->treeFromDistMatrix(distMat, alignPtr, seq1, nSeqs, copyOfPhylipName);

        *phylipName = copyOfPhylipName;
        // AW: message outputted by OutputFile function
        // utilityObject->info("Guide tree        file created:   [%s]",
        //                              phylipName->c_str());
        delete clusterTree;
    }
    *success = true;
}

auto_ptr<AlignmentSteps>
TreeInterface::getWeightsAndStepsFromTreeNJ(Alignment* alignPtr, 
                                            DistMatrix* distMat, string* treeName,
                                            vector<int>* seqWeights, int fSeq, int numSeqs, 
                                            bool* success)
{
    auto_ptr<AlignmentSteps> progSteps;
    Tree groupTree; 
    if(numSeqs == 1)
    {
        utilityObject->info("Only 1 sequence, cannot do multiple alignment\n");
        *success = false;
        return progSteps;
    }
    int status;
    status = readTreeAndCalcWeightsNJ(&groupTree, alignPtr, distMat, treeName, 
                                      seqWeights, fSeq, numSeqs);
    
    if (status == 0)
    {
        *success = false;
        return progSteps;
    }
    
    progSteps = groupTree.createSets(0, numSeqs);
    int _numSteps = progSteps->getNumSteps();
    utilityObject->info("There are %d groups", _numSteps);        
    // clear the memory used for the phylogenetic tree

    if (numSeqs >= 2)
    {
        groupTree.clearTree(NULL);
    }
    
    *success = true;            
    
    return progSteps;

}

/**
 * UPGMA functions
 *
 */
 
 
auto_ptr<AlignmentSteps>
TreeInterface::getWeightsAndStepsFromDistMatUPGMA(vector<int>* seqWeights, 
                                                  DistMatrix* distMat, Alignment *alignPtr, 
                                                  int seq1, int nSeqs, string* phylipName, bool* success)
{
    auto_ptr<AlignmentSteps> progSteps;
    RootedGuideTree guideTree;   

    progSteps = generateTreeFromDistMatUPGMA(&guideTree, distMat, alignPtr, seq1, nSeqs,
                                             phylipName, success);

    guideTree.calcSeqWeights(0, nSeqs, seqWeights);
    
    distMat->makeSimilarityMatrix();
                                                    
    return progSteps;
}                                                                                             

auto_ptr<AlignmentSteps> TreeInterface::generateTreeFromDistMatUPGMA(RootedGuideTree* guideTree, DistMatrix* distMat, Alignment *alignPtr, int seq1, int nSeqs, 
                                            string* phylipName, bool* success)
{
    auto_ptr<AlignmentSteps> progSteps;
    string copyOfPhylipName = string(*phylipName);     
    
    if (nSeqs >= 2) 
    {
        RootedClusterTree clusterTree;   
        progSteps = clusterTree.treeFromDistMatrix(guideTree, distMat, alignPtr, seq1,
                                                   nSeqs, copyOfPhylipName);
                                        
        *phylipName = copyOfPhylipName;
        // AW: message outputted by OutputFile function
        // utilityObject->info("Guide tree        file created:   [%s]",
        //                              phylipName->c_str());        
    }
    *success = true;
    return progSteps;
}
                                  
void TreeInterface::getWeightsFromDistMatUPGMA(vector<int>* seqWeights, DistMatrix* distMat, 
                                   Alignment *alignPtr, int seq1, int nSeqs, 
                                   string* phylipName, bool* success)
{
    getWeightsAndStepsFromDistMatUPGMA(seqWeights, distMat, alignPtr, 
                                 seq1, nSeqs, phylipName, success);
}

/**
 * The function getWeightsForProfileAlignUPGMA is used to generate the sequence weights
 * that will be used in a profile alignment. It also recalculates the distance matrix, and 
 * then returns it as a similarity matrix. This function uses the NJ code to read in a tree
 * from a file if we are using a previous tree. This is because that part of the code is able
 * to do the finding of the root etc.
 */
void TreeInterface::getWeightsForProfileAlignUPGMA(Alignment* alignPtr, DistMatrix* distMat, 
                             string* p1TreeName, vector<int>* p1Weights, string* p2TreeName, 
                    vector<int>* p2Weights, int numSeqs, int profile1NumSeqs, bool useTree1, 
                    bool useTree2, bool* success)
{
    int status = 0;
    
    if(useTree1)
    {
        // Use the code to read in the tree and get the seqWeights
        Tree groupTree1;
        if (profile1NumSeqs >= 2)
        {
            status = groupTree1.readTree(alignPtr, p1TreeName->c_str(), 0, profile1NumSeqs);
            if (status == 0)
            {
                *success = false;
                return;
            }
        }
        groupTree1.calcSeqWeights(0, profile1NumSeqs, p1Weights);
    
        if (profile1NumSeqs >= 2)
        {
            groupTree1.clearTree(NULL);
        }                
    }
    else
    {
        if (profile1NumSeqs >= 2) 
        {        
            RootedGuideTree guideTree;
            generateTreeFromDistMatUPGMA(&guideTree, distMat, alignPtr, 1, profile1NumSeqs, 
                                         p1TreeName, success);
            guideTree.calcSeqWeights(0, profile1NumSeqs, p1Weights);
        }                                      
    }
    
    if(useTree2)
    {
        Tree groupTree2;
        if (numSeqs - profile1NumSeqs >= 2)
        {
            status = groupTree2.readTree(alignPtr, p2TreeName->c_str(), 
                                         profile1NumSeqs, numSeqs);
            if (status == 0)
            {
                *success = false;
                return;
            }
        }
        groupTree2.calcSeqWeights(profile1NumSeqs, numSeqs, p2Weights);
    
        if (numSeqs - profile1NumSeqs >= 2)
        {
            groupTree2.clearTree(NULL);
        }    
    }
    else
    {
        if(numSeqs - profile1NumSeqs >= 2) 
        {        
            RootedGuideTree guideTree;
            generateTreeFromDistMatUPGMA(&guideTree, distMat, alignPtr, profile1NumSeqs + 1,
                                         numSeqs - profile1NumSeqs, p2TreeName, success);
            guideTree.calcSeqWeights(profile1NumSeqs, numSeqs, p2Weights);
        }    
    }
    
    if (userParameters->getNewTree1File() || userParameters->getNewTree2File()) 
    {
        *success = false;
        return;
    }    
    
    MSA* msaObj = new MSA();
        
    int count = msaObj->calcPairwiseForProfileAlign(alignPtr, distMat);
    
    delete msaObj;
    
    if (count == 0) 
    {
        *success = false;
        return;
    }
    
    /**
     * Convert distances to similarities!!!!!!
     */   
    distMat->makeSimilarityMatrix();
      
    *success = true;

}

/**
 * Mark 10-5-2007: Bug fix # 42
 * The following function is to be used only for calculating the weights for the Qt
 * low scoring segments.
 */
void TreeInterface::getWeightsForQtLowScore(vector<int>* seqWeights, DistMatrix* distMat, 
                                   Alignment *alignPtr, int seq1, int nSeqs, 
                                   string* phylipName, bool* success)
{
    string copyOfPhylipName = string(*phylipName);     
    int _nSeqs = nSeqs;
    int status = 0;
    
    generateTreeFromDistMatNJ(distMat, alignPtr, seq1, _nSeqs, phylipName, success);
    
    Tree groupTree;

    status = 0;
    if (_nSeqs >= 2)
    {
        status = groupTree.readTree(alignPtr, phylipName->c_str(), seq1 - 1,
                                    seq1 + _nSeqs-1); // mark 22-5-07
        if(status == 0)
        {
            *success = false;
            return;
        }
        else
        {
            *success = true;
        }
    }

    groupTree.calcSeqWeights(seq1 - 1, seq1 + _nSeqs - 1, seqWeights); // mark 23-5-07
}
                                                           
}
