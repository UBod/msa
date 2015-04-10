/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef GDEFILEPARSER_H
#define GDEFILEPARSER_H

#include <string>
#include "FileParser.h"

namespace clustalw
{

class GDEFileParser : public FileParser
{
    public:
        /* Functions */
        GDEFileParser(string filePath);
        ~GDEFileParser();
        virtual vector<Sequence> getSeqRange(int firstSeq, int num, string *offendingSeq=NULL);
        virtual Sequence getSeq(int seqNum, string *offendingSeq=NULL);
        virtual int countSeqs();
        virtual void getSecStructure(vector<char>& gapPenaltyMask, 
                vector<char>& secStructMask,
               string& secStructName, int &structPenalties, int length); 

        /* Attributes */

    private:
        /* Functions */

        /* Attributes */
        string fileName;
};

}
#endif

