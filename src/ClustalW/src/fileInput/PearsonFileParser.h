/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/** 
 * This file is for parsing pearson format files.
 * CHANGE: 
 */
#ifndef PEARSONFILEPARSER_H
#define PEARSONFILEPARSER_H

#include <string>
#include "FileParser.h"

namespace clustalw
{
using namespace std;

class PearsonFileParser : public FileParser
{
    public:
        /* Functions */
        PearsonFileParser(string filePath);
        virtual vector<Sequence> getSeqRange(int firstSeq, int num, string *offendingSeq=NULL);
        virtual Sequence getSeq(int seqNum, string *offendingSeq=NULL);
        virtual int countSeqs();
        virtual void getSecStructure(vector<char>& gapPenaltyMask, 
                                     vector<char>& secStructMask, string& secStructName, 
                                     int &structPenalties, int length); 

        /* Attributes */

    private:
        /* Functions */

        /* Attributes */
        string fileName;
};

}
#endif

