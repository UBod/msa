/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef FILEREADER_H
#define FILEREADER_H

#include <vector>
#include <string>
#include <memory>
#include "../alignment/Alignment.h"
#include "../alignment/Sequence.h"
#include "../general/userparams.h"
#include "../general/utils.h"
#include "FileParser.h"
#include "ClustalFileParser.h"
#include "PearsonFileParser.h"
#include "PIRFileParser.h"
#include "GDEFileParser.h"
#include "MSFFileParser.h"
#include "RSFFileParser.h"
#include "EMBLFileParser.h"

namespace clustalw
{

class FileReader
{
    public:
        /* Functions */
        FileReader();
        ~FileReader();
        int seqInput(Alignment* alignPtr, bool append, string *offendingSeq);
        int readSeqs(Alignment* alignPtr, int firstSeq, string *offendingSeq);
        int readCharacterSeqs(Alignment* alignPtr, int firstSeq, string *offendingSeq, ClustalWInput *input);
        int profileInput(Alignment* alignPtr);

        /* Attributes */

    private:
        /* Functions */
        void checkInfile(int* nseqs, auto_ptr<FileParser>& fileParser);
        /* Attributes */
        string sequenceFileName;
        bool noEmptySequence(vector<Sequence> seqRangeVector, string *offendingSeq);
            
        InFileStream* fileIn;
        int structPenalties;
        string secStructName;
        vector<char> secStructMask; // Will need to be cleared out before every reading!
        vector<char> gapPenaltyMask;
        vector<string> formatNames; 
};
}
#endif

