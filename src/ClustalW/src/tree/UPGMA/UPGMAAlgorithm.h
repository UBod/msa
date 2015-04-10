#ifndef __UPGMAALGORITHM_H
#define __UPGMAALGORITHM_H

#include "Node.h" 
#include <limits>
#include "../../general/clustalw.h"
#include "../AlignmentSteps.h"
#include "RootedGuideTree.h"
#include <list>
#include <memory>

namespace clustalw
{

using namespace std;

class UPGMAAlgorithm
{
    public:
        bool overwriteMatrix; 

        UPGMAAlgorithm();
        auto_ptr<AlignmentSteps> generateTree(RootedGuideTree* phyTree, 
                                    DistMatrix* distMat, SeqInfo* seqInfo,
                                    bool overwrite, ofstream* tree = 0);
        void setVerbose(bool _verbose){verbose = _verbose;}                            
                  
    private:
        Node **initialiseNodes(double *distanceMatrix, int firstSeq);

        Node *doUPGMA(Node **nodes, ofstream* tree);
        void printAllNodes(Node** nodes);
        void addAlignmentStep(vector<int>* group1, vector<int>* group2);
        Node** getNodeWithMinDist(Node** clusters);
        void recomputeNodeToJoin1DistMatRow(Node* nodeToJoin1, double** nodeToJoin2DistIter);
        void computeAllOtherDistsToNewNode(Node* nodeToJoin1, Node* nodeToJoin2,
                                                   double** nodeToJoin2DistIter);
        void computeDistsUpToNodeToJoin2(Node* nToJoin1, Node* nToJoin2, 
                                         double** nodeToJoin2DistIter);
        void computeDistsForNodesAfterNode2(Node* nToJoin2);
        void movePtrPastUnusedDistances(double** ptrToDist)
        {            
            while(**ptrToDist < 0)
            {
                (*ptrToDist)++;
            }
        }
        double calcNewDist(double dist1, double dist2);                             
        //vector<Node*>::iterator getIterToNodeWithMinDist(vector<Node*>* nodesLeft);
        
        auto_ptr<AlignmentSteps> progSteps;
        int numSeqs; 
        bool verbose;
        int orderNode1, orderNode2, orderNewNode;
};

}

#endif
