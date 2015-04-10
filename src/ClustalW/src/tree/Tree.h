/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef TREE_H
#define TREE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "../alignment/Alignment.h"
#include "AlignmentSteps.h"
#include <memory>
namespace clustalw
{

class Tree
{
    public:
        /* Functions */
        void calcSeqWeights(int firstSeq, int lastSeq, vector<int>* sweight);
        int readTree(Alignment* alignPtr, const string& treeFileName, int firstSeq, 
                     int lastSeq);
        auto_ptr<AlignmentSteps> createSets(int firstSeq, int lastSeq);
        int calcSimilarities(Alignment* alignPtr, DistMatrix* distMat);
        void clearTree(TreeNode* p);
        /* Attributes */

    private:
        /* Functions */
        void createTree(TreeNode* ptree, TreeNode* parent, ifstream* file);
        void createNode(TreeNode* pptr, TreeNode* parent);
        TreeNode* insertNode(TreeNode* pptr);
        void clearTreeNodes(TreeNode* p);
        TreeNode* reRoot(TreeNode* ptree, int nseqs);
        TreeNode* insertRoot(TreeNode* p, float diff);
        float calcRootMean(TreeNode* root, float *maxDist);
        float calcMean(TreeNode* nptr, float *maxDist, int nSeqs);
        void orderNodes();
        int calcWeight(int leaf);
        void skipSpace(ifstream* file);
        void groupSeqs(TreeNode* p, int *nextGroups, int nSeqs, AlignmentSteps* stepsPtr);
        void markGroup1(TreeNode* p, int *groups, int n);
        void markGroup2(TreeNode* p, int *groups, int n);
        TreeNode* avail();
        void setInfo(TreeNode* p, TreeNode* parent, int pleaf, string pname, float
                     pdist);
        void debugPrintAllNodes(int nseqs);
            
        /* Attributes */
        AlignmentSteps progSteps;
        char charFromFile;
        ifstream file;
        TreeNode** lptr;
        TreeNode** olptr;
        TreeNode** nptr;
        TreeNode** ptrs;
        int nnodes;
        int ntotal;
        bool rootedTree;
        TreeNode* seqTree;
        TreeNode* root;
        int* groups;
        int numSeq;
        int numSets;
        const static int MAXERRS = 10;
};

}
#endif
