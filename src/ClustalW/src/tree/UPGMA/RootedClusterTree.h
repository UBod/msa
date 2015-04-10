/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef ROOTEDCLUSTERTREE_H
#define ROOTEDCLUSTERTREE_H
#include "../ClusterTree.h"
#include "RootedGuideTree.h"
#include "../AlignmentSteps.h"
namespace clustalw
{

class RootedClusterTree : private ClusterTree
{
    public:    
        //RootedClusterTree();
        /**
         * NOTE these will have different signatures!!!!
         */ 
        void treeFromAlignment(TreeNames* treeNames, Alignment *alignPtr);
        auto_ptr<AlignmentSteps> treeFromDistMatrix(RootedGuideTree* phyloTree, 
                                        DistMatrix* distMat, Alignment *alignPtr, int seq1, 
                                        int nSeqs, string& phylipName);
};

}
#endif
