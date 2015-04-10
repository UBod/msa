/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef CLUSTERTREEALGORITHM_H
#define CLUSTERTREEALGORITHM_H

#include <memory>
#include <fstream>
#include <iostream>
#include "../general/clustalw.h"

namespace clustalw
{

class ClusterTreeAlgorithm
{
    public:
  virtual ~ClusterTreeAlgorithm(){};

        virtual void generateTree(clustalw::PhyloTree* phyTree, clustalw::DistMatrix* distMat, clustalw::SeqInfo* seqInfo,
                                  ofstream* tree = 0) = 0;
        virtual void setVerbose(bool choice) = 0;
};

}
#endif
