/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef EMBLFILEPARSER_H
#define EMBLFILEPARSER_H

#include <string>
#include "FileParser.h"

namespace clustalw
{

class EMBLFileParser : public FileParser
{
    public:
        /* Functions */
        EMBLFileParser(string filePath);
        ~EMBLFileParser();
        virtual Sequence getSeq(int seqNum, string *offendingSeq=NULL);
        virtual vector<Sequence> getSeqRange(int firstSeq, int num, string *offendingSeq=NULL);
        virtual int countSeqs();
        virtual void getSecStructure(vector<char>& gapPenaltyMask, 
                vector<char>& secStructMask,
               string& secStructName, int &structPenalties, int length); 

        /* Attributes */

    private:
        /* Functions */
        void getSwissFeature(char* line, vector<char>& secStructMask, int length);
        void getSwissMask(char* line, vector<char>& gapPenaltyMask, int length);
        /* Attributes */
        string fileName;
};

}
#endif

