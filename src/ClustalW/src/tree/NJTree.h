/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef NJTREE_H
#define NJTREE_H

#include <vector>
#include <iomanip>
#include "ClusterTreeAlgorithm.h"
#include "../general/userparams.h"

namespace clustalw
{

class NJTree : public ClusterTreeAlgorithm
{
    public:
        NJTree(): verbose(false){};
	virtual ~NJTree(){};

        /** calculate an NJ tree
         *
         * @param phyTree the tree structure
         * @param distMat distance matrix
         * @param seqInfo holding sequence number info
         * @param log ofstream to log info to (used by -outputtree)
         *
         */
        virtual void generateTree(clustalw::PhyloTree* phyTree,
                                  clustalw::DistMatrix* distMat,
                                  clustalw::SeqInfo* seqInfo,
                                  ofstream* log = 0);
        /** be verbose during tree generation
         *
         *  if set to true, generateTree will need a log ofstream
         */
        virtual void setVerbose(bool choice){verbose = choice;};
    private:
        vector<double> av;
        vector<int> tkill;
        bool verbose;
};

}
#endif
