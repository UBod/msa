/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Changes:
 * Mark 10-5-2007: Bug fix # 42. Added getWeightsForQtLowScore function.
 *
 *
 */
#ifndef TREEINTERFACE_H
#define TREEINTERFACE_H

#include <vector>
#include <string>
#include <memory>
#include "AlignmentSteps.h"
#include "../alignment/Alignment.h"
#include "../general/clustalw.h"

namespace clustalw
{

using namespace std;
class Tree;
class RootedGuideTree;

class TreeInterface
{
    public:
        /**
         * This function will be used to generate the phylogenetic tree from the distMat using
         * either UPGMA or NJ. It will then calculate the seqWeights and the steps. 
         * Note the Distmat
         * will be in similarity form after this. It will no longer be distances!!!!!
         */
        auto_ptr<AlignmentSteps> getWeightsAndStepsFromDistMat(vector<int>* seqWeights, 
                                                          DistMatrix* distMat, 
                                                          Alignment *alignPtr, 
                                                          int seq1, int nSeqs, 
                                                          string* phylipName, bool* success);
        /**
         * This will be called by sequencesAlignToProfile and QTcalcWeightsForLowScoreSeg
         */        
        void getWeightsFromDistMat(vector<int>* seqWeights, DistMatrix* distMat, 
                                   Alignment *alignPtr, int seq1, int nSeqs, 
                                   string* phylipName, bool* success);
                                   
        void getWeightsForQtLowScore(vector<int>* seqWeights, DistMatrix* distMat, 
                                   Alignment *alignPtr, int seq1, int nSeqs, 
                                   string* phylipName, bool* success);
                                   
        /**
         * This function will be called from doAlignUseOldTree
         */                                                  
        auto_ptr<AlignmentSteps> getWeightsAndStepsFromTree(Alignment* alignPtr, 
                                                DistMatrix* distMat, string* treeName,
                                                vector<int>* seqWeights, int fSeq, 
                                                int numSeqs, bool* success);
                                                             
        /**
         * This function will be called from sequencesAlignToProfile, it doesnt calc the
         * steps.
         */                                                     
        int getWeightsFromGuideTree(Alignment* alignPtr, DistMatrix* distMat,
                                    string* treeName, vector<int>* seqWeights, int fSeq,
                                    int nSeqs, bool* success);
        /**
         * This function is used to generate 2 guide trees for the profile align part.
         * IT must be done on its own. It calls calcPairwiseForProfileAlign in MSA.
         * It is called by profileAlign and removeFirstIterate.
         */
        void getWeightsForProfileAlign(Alignment* alignPtr, DistMatrix* distMat, 
                             string* p1TreeName, vector<int>* p1Weights, string* p2TreeName, 
                    vector<int>* p2Weights, int numSeqs, int profile1NumSeqs, bool useTree1, 
                    bool useTree2, bool* success); 
        /**
         * This function is used to generate the guide tree from the distance matrix. 
         * It doesnt return
         * any seqWeights or AlignmentSteps. It will be used by doGuideTreeOnly
         */
        void generateTreeFromDistMat(DistMatrix* distMat, Alignment *alignPtr, 
                                                          int seq1, int nSeqs, 
                                                          string* phylipName, bool* success);
        
        /**
         * This function is to simply to call either the UPGMA or NJ version of 
         * this function, and print out all the trees. It will be called 
         * from phylogeneticTree in the clustal class.
         */                                                  
        void treeFromAlignment(TreeNames* treeNames, Alignment *alignPtr);
        
        /**
         * This function will be used to bootstrap the tree and output the results. 
         * It will use either
         * UPGMA or NJ for the bootstrapping, depending on what is selected.
         */
        void bootstrapTree(TreeNames* treeNames, Alignment *alignPtr);
         
    private:
        auto_ptr<AlignmentSteps> getWeightsAndStepsFromDistMatNJ(vector<int>* seqWeights, 
                                                          DistMatrix* distMat, 
                                                          Alignment *alignPtr, 
                                                          int seq1, int nSeqs, 
                                                          string* phylipName, bool* success);
        
        auto_ptr<AlignmentSteps> getWeightsAndStepsUseOldGuideTreeNJ(DistMatrix* distMat, 
                                                   Alignment *alignPtr,  string* treeName,
                                                   vector<int>* seqWeights, 
                                                   int fSeq, int nSeqs, bool* success);
                                                             
        int readTreeAndCalcWeightsNJ(Tree* groupTree, Alignment* alignPtr, 
                            DistMatrix* distMat, string* treeName, vector<int>* seqWeights,
                            int fSeq, int nSeqs);
        
        int getWeightsFromGuideTreeNJ(Alignment* alignPtr, DistMatrix* distMat,
                                      string* treeName, vector<int>* seqWeights, int fSeq,
                                      int nSeqs, bool* success);
        
        void getWeightsFromDistMatNJ(vector<int>* seqWeights, DistMatrix* distMat, 
                                   Alignment *alignPtr, int seq1, int nSeqs, 
                                   string* phylipName, bool* success); 
        
        void getWeightsForProfileAlignNJ(Alignment* alignPtr, DistMatrix* distMat, 
                             string* p1TreeName, vector<int>* p1Weights, string* p2TreeName, 
                    vector<int>* p2Weights, int numSeqs, int profile1NumSeqs, bool useTree1, 
                    bool useTree2, bool* success);
                    
        void generateTreeFromDistMatNJ(DistMatrix* distMat, Alignment *alignPtr, 
                                 int seq1, int nSeqs, string* phylipName, bool* success);
        
        auto_ptr<AlignmentSteps> getWeightsAndStepsFromTreeNJ(Alignment* alignPtr, 
                               DistMatrix* distMat, string* treeName,
                               vector<int>* seqWeights, int fSeq, int numSeqs, bool* success);
        
        /** UPGMA functions */
        auto_ptr<AlignmentSteps> getWeightsAndStepsFromDistMatUPGMA(vector<int>* seqWeights, 
                                 DistMatrix* distMat, Alignment *alignPtr, 
                                 int seq1, int nSeqs, string* phylipName, bool* success);
                                 
        auto_ptr<AlignmentSteps> generateTreeFromDistMatUPGMA(RootedGuideTree* guideTree,
                             DistMatrix* distMat, Alignment *alignPtr, int seq1, int nSeqs, 
                                            string* phylipName, bool* success);
                                 
        void getWeightsFromDistMatUPGMA(vector<int>* seqWeights, DistMatrix* distMat, 
                                   Alignment *alignPtr, int seq1, int nSeqs, 
                                   string* phylipName, bool* success); 
                                   
        void getWeightsForProfileAlignUPGMA(Alignment* alignPtr, DistMatrix* distMat, 
                             string* p1TreeName, vector<int>* p1Weights, string* p2TreeName, 
                    vector<int>* p2Weights, int numSeqs, int profile1NumSeqs, bool useTree1, 
                    bool useTree2, bool* success);                                
};
}
#endif
