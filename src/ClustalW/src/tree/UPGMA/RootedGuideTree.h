/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Changes:
 * mark May 8th 2007: removed makeEmpty(Node* t function), changed makeEmpty(). I also
 * made changes to any functions that accessed Node's data members 
 */
#ifndef ROOTEDGUIDETREE_H
#define ROOTEDGUIDETREE_H

#include "Node.h"
#include <fstream>

#include "../../alignment/Alignment.h"

namespace clustalw
{

using namespace std;

class RootedGuideTree
{
    public:
        RootedGuideTree();
        RootedGuideTree(Node* root);
        ~RootedGuideTree();
        void setRoot(Node* r);
        void makeEmpty();
        void calcSeqWeights(int firstSeq, int lastSeq, vector<int>* seqWeights);
        Node* getRoot(){return root;}
    private:
        void orderNodes();
        int calcOrderNode(Node* node);
        void calcWeights(vector<int>* seqWeights);
        void doWeightCalc(float weightSoFar, vector<float>* seqWeights, Node* t);
        Node* root;
};

}
#endif
