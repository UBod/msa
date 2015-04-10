/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/** The aim of this class is to return one sequence at a time.
 * Note that the file must be open when it is passed to the FileParser.
 * The parser does not know the name of the file to open. Only the filereader knows.
 *
 * Changes: 
 *
 * Mark 24-1-2007. I added the function findDelimiter to determine if '\r' or '\n' 
 * will be used as the line delimiter when parsing the file.
 *
 * 10-02-07,Nigel Brown(EMBL): Removed delimiter and findDelimiter()
 * members, as functionality now handled by the stream class.
 */
#ifndef FILEPARSER_H
#define FILEPARSER_H

#include <ctype.h>
#include "../alignment/Sequence.h"
#include "../general/userparams.h"
#include <iostream>
#include "InFileStream.h"
#include "../RClustalW.h"

namespace clustalw
{
 
class FileParser
{
    public:
        /* Functions */
        FileParser();
        virtual ~FileParser();
        virtual vector<Sequence> getSeqRange(int firstSeq, int num, string *offendingSeq) = 0;
        virtual Sequence getSeq(int seqNum, string *offendingSeq=NULL) = 0;
        virtual int countSeqs() = 0; // VIRTUAL 
        virtual void getSecStructure(vector<char>& gapPenaltyMask, 
                                     vector<char>& secStructMask,
                                     string& secStructName, int &structPenalties, int length) = 0;

        vector<Sequence> getSeqRangeR(int firstSeq, int num, string *offendingSeq, ClustalWInput *input);

        void fillCharTab(void);
        char getDelimiter(string filename);
        /* Attributes */
        char chartab[128];
        int getParseExitCode() { return parseExitCode; };
        
    protected:
        void freeFileResources(InFileStream* filePtr);
        InFileStream* _fileIn;
        int parseExitCode; // reason for returning empty sequence
                           // vector; same as used in FileReader

    private:
        /* Functions */

        /* Attributes */

};

}
#endif

