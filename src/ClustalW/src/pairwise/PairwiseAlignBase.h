/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef PAIRWISEALIGNBASE_H
#define PAIRWISEALIGNBASE_H

#include <vector>
#include "../alignment/Alignment.h"

namespace clustalw
{

class PairwiseAlignBase
{
    public:
  virtual ~PairwiseAlignBase(){};
        /* Functions */
        virtual void pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, 
                                   int iEnd, int jStart, int jEnd) = 0; 
        /* Attributes */

    private:
        /* Functions */

        /* Attributes */
};

}
#endif
