/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef CLUSTERTREEOUTPUT_H
#define CLUSTERTREEOUTPUT_H

#include <memory>
#include <fstream>
#include "../alignment/Alignment.h"
#include "../general/clustalw.h"

namespace clustalw
{

class ClusterTreeOutput
{
    public:
        ClusterTreeOutput(clustalw::SeqInfo* seqInfo, int boot);
        void printNexusTree(clustalw::PhyloTree* phyloTree, ofstream* tree,
                   clustalw::Alignment *alignPtr, clustalw::DistMatrix* distMat, vector<int>* bootTotals);
        void printTree(clustalw::PhyloTree* phyloTree, ofstream* tree, vector<int>* totals);
        void printPhylipTree(clustalw::PhyloTree* phyloTree, ofstream* tree,
                   clustalw::Alignment *alignPtr, clustalw::DistMatrix* distMat, vector<int>* bootTotals);
        void printTreeDesc(clustalw::PhyloTree* phyloTree);
          
    private:
        ClusterTreeOutput(); // Dont allow contruction with default!!!!
        int twoWaySplit(clustalw::PhyloTree* phyloTree, ofstream* tree, int startRow, 
                    int flag, clustalw::Alignment *alignPtr, vector<int>* bootTotals);
        int twoWaySplitNexus(clustalw::PhyloTree* phyloTree, ofstream* tree, int startRow,
                    int flag, clustalw::Alignment *alignPtr, vector<int>* bootTotals);
        /* Attributes! */
        int firstSeq;
        int lastSeq;
        int numSeqs;
        int bootstrap;
};

}

#endif
