/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef PROFILEALIGNALGORITHM_H
#define PROFILEALIGNALGORITHM_H

#include <vector>
#include "../alignment/Alignment.h"
namespace clustalw
{

class ProfileAlignAlgorithm
{
    public:
  virtual ~ProfileAlignAlgorithm(){};

    /* Functions */
    virtual int profileAlign(Alignment* alnPtr, DistMatrix* distMat, vector<int>* group, 
                             int* aligned) = 0;
    /* Attributes */

    protected:
    /* Attributes */
        int prfLength1;
        int prfLength2;
        SeqArray seqArray;
        vector<int> alnWeight;
        int nseqs1;
        int nseqs2;     
    private:
    /* Functions */
    
    /* Attributes */
       
};

}
#endif
