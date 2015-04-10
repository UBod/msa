/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef CLUSTALFILEPARSER_H
#define CLUSTALFILEPARSER_H

#include <string>
#include "FileParser.h"

namespace clustalw
{

class ClustalFileParser : public FileParser
{
    public:
        /* Functions */
        ClustalFileParser(string filePath);
        ~ClustalFileParser();
        virtual Sequence getSeq(int seqNum, string *offendingSeq=NULL);
        virtual vector<Sequence> getSeqRange(int firstSeq, int num, string *offendingSeq=NULL);
        virtual int countSeqs();
        virtual void getSecStructure(vector<char>& gapPenaltyMask, 
                vector<char>& secStructMask,
               string& secStructName, int &structPenalties, int length); 

        /* Attributes */

    private:
        /* Functions */
        bool clustalBlankline(char* line); // Only used in this class!
        string fileName;
        /* Attributes */
        
};

}
#endif
