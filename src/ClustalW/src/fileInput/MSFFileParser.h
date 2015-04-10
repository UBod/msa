/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef MSFFILEPARSER_H
#define MSFFILEPARSER_H

#include <string>
#include "FileParser.h"

namespace clustalw
{

class MSFFileParser : public FileParser
{
    public:
        /* Functions */
        MSFFileParser(string filePath);
        virtual Sequence getSeq(int seqNum, string *offendingSeq=NULL);
        virtual vector<Sequence> getSeqRange(int firstSeq, int num, string *offendingSeq=NULL);
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


