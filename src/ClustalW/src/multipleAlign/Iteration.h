/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef ITERATION_H
#define ITERATION_H
#include "../alignment/Alignment.h"
#include "../general/clustalw.h"

namespace clustalw
{

class Iteration
{
    public:
        //Iteration();
        bool iterationOnTreeNode(int numSeqsProf1, int numSeqsProf2, int& prfLength1,
                                    int& prfLength2, SeqArray* seqArray);
        bool removeFirstIterate(Alignment* alnPtr);
    private:
        void printSeqArray(SeqArray* arrayToPrint);
};

}

#endif

