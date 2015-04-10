/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Changes:
 * mark 8-5-2007: I removed the isLeafNode function. Changed both the traversal functions
 * to access Node's data members via functions.
 */
#ifndef ROOTEDTREEOUTPUT_H
#define ROOTEDTREEOUTPUT_H

#include <fstream>
#include "RootedGuideTree.h"
#include "../../general/clustalw.h"


/** this is only used for upgma?!
 *
 */


namespace clustalw
{

class RootedTreeOutput
{
    public:
        RootedTreeOutput(SeqInfo* seqInfo);
        void printPhylipTree(RootedGuideTree* tree, ofstream* ptrToFile, Alignment *alignPtr,
                             DistMatrix* distMat);
        void printNexusTree(RootedGuideTree* tree, ofstream* ptrToFile, Alignment *alignPtr, 
                            DistMatrix* distMat);                   
                
    private:
        void phylipTraverse(ofstream* ptrToFile, Alignment *alignPtr, Node* tree);
        void nexusTraverse(ofstream* ptrToFile, Alignment *alignPtr, Node* tree);
        int firstSeq;
        int lastSeq;
        int numSeqs;                  
};

}

#endif
