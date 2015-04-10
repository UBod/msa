/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <iostream>
#include "RootedGuideTree.h"
#include "../../general/userparams.h"
using namespace std;

namespace clustalw
{

RootedGuideTree::RootedGuideTree()
 : root(0){ }

RootedGuideTree::RootedGuideTree(Node* r)
 : root(r) { }

RootedGuideTree::~RootedGuideTree()
{
    root->makeEmpty();
}

void RootedGuideTree::setRoot(Node* r)
{
    makeEmpty();
    
    root = r;
}
 
void RootedGuideTree::makeEmpty()
{
    root->makeEmpty();
}

void RootedGuideTree::orderNodes()
{
    calcOrderNode(root);
}

int RootedGuideTree::calcOrderNode(Node* node)
{
    if(node != 0)
    {
        if(node->getLeft() == 0 && node->getRight() == 0) // Leaf Node
        {
            node->setOrder(1);
            return 1;
        }
        else // Internal node
        {
            node->setOrder(calcOrderNode(node->getLeft()) + calcOrderNode(node->getRight()));
            return node->getOrder();            
        }
    }
    return 0;
}

void RootedGuideTree::calcSeqWeights(int firstSeq, int lastSeq, vector<int>* seqWeights)
{
    if((int)seqWeights->size() < lastSeq - 1)
    {
        seqWeights->resize(lastSeq - 1);
    }
    
    int i = 0, _nSeqs = 0;
    int temp = 0, sum = 0;
    //
    // If there are more than three sequences....
    //
    _nSeqs = lastSeq - firstSeq;
    if ((_nSeqs >= 2) && (userParameters->getDistanceTree() == true) && 
        (userParameters->getNoWeights() == false))
    {
        //
        // Calculate sequence weights based on Phylip tree.
        //
        orderNodes();
        calcWeights(seqWeights);

        //
        // Normalise the weights, such that the sum of the weights = INT_SCALE_FACTOR
        //

        sum = 0;
        for (i = firstSeq; i < lastSeq; i++)
        {
            sum += (*seqWeights)[i];
        }

        if (sum == 0)
        {
            for (i = firstSeq; i < lastSeq; i++)
            {
                (*seqWeights)[i] = 1;
            }
            sum = i;
        }

        for (i = firstSeq; i < lastSeq; i++)
        {
            (*seqWeights)[i] = ((*seqWeights)[i] * INT_SCALE_FACTOR) / sum;
            if ((*seqWeights)[i] < 1)
            {
                (*seqWeights)[i] = 1;
            }
        }
    }
    else
    {
        //
        // Otherwise, use identity weights.
        //
        temp = INT_SCALE_FACTOR / _nSeqs;
        // AW 2009-07-09: goes wrong if we have more than
        // INT_SCALE_FACTOR seqs. if so, set to 1, just as above
        // same as in Tree.cpp
        if (temp < 1)
            temp = 1;

        
        for (i = firstSeq; i < lastSeq; i++)
        {
            (*seqWeights)[i] = temp;
        }
    }
}

void RootedGuideTree::calcWeights(vector<int>* seqWeights)
{
    vector<float> weights;
    int sizeSeqWeights = seqWeights->size();
    weights.resize(sizeSeqWeights, 0.0);
    
    doWeightCalc(0.0, &weights, root);

    for(int i = 0; i < sizeSeqWeights; i++)
    {
        (*seqWeights)[i] = static_cast<int>(weights[i] * 100);        
    }
}

void RootedGuideTree::doWeightCalc(float weightSoFar, vector<float>* weights, Node* t)
{
    if(t != 0)
    {
        if(t->getLeft() == 0 && t->getRight() == 0) // Leaf Node
        {
            (*weights)[t->getSeqNum() - 1] = weightSoFar;
        }
        else // Internal node
        {
            float w = weightSoFar + (t->getHeight() / t->getOrder());
            doWeightCalc(w, weights, t->getLeft());
            doWeightCalc(w, weights, t->getRight());            
        }
    }
}

}
