/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef UNROOTEDCLUSTERTREE_H
#define UNROOTEDCLUSTERTREE_H
#include "ClusterTree.h"

namespace clustalw
{

class UnRootedClusterTree : private ClusterTree
{
    public:    
        UnRootedClusterTree();
        void treeFromAlignment(TreeNames* treeNames, Alignment *alignPtr);
        void treeFromDistMatrix(DistMatrix* distMat, Alignment *alignPtr, int seq1, 
                                int nSeqs, string& phylipName);
        void bootstrapTree(TreeNames* treeNames, Alignment *alignPtr);
    private:
        PhyloTree* phyloTree;
        ClusterTreeOutput* outputTree;    
};

}
#endif
