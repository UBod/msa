/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * 10-02-07,Nigel Brown(EMBL): Removed delimiter and findDelimiter()
 * members, as functionality now handled by the stream class.
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "FileParser.h"
const char LF = 0x0a;  //linefeed
const char CR = 0x0d;  //carriage return





namespace clustalw
{
FileParser::FileParser()
{
    parseExitCode = OK;
}

FileParser::~FileParser()
{
    // Dont do anything!!!!
}

vector<Sequence> FileParser::getSeqRangeR(int firstSeq, int num, string *offendingSeq, ClustalWInput *input) {
	vector<Sequence> seqRangeVector;
	    int i;

		int n = (*input).inputSeqs.size();

		string title("msa with RClustalW");

	    for (i=0; i < num; i++) {
	        if ((int)(*input).inputSeqs[i].length() > userParameters->getMaxAllowedSeqLength()) {
	        	seqRangeVector.clear();
	        	return seqRangeVector;
	        }

	        Sequence tempSeq((*input).inputSeqs[i], (*input).seqNames[i], title);
	        seqRangeVector.push_back(tempSeq);
	    }
	    return seqRangeVector;
}

void FileParser::fillCharTab(void)
{
    register int i;
    register char c;

    for (i = 0; i < 128; chartab[i++] = 0)
        ;
    for (i = 0; i <= userParameters->getMaxAA() + 1; i++)
    {
        c = userParameters->getAminoAcidCode(i);
        chartab[(int)c] = chartab[tolower(c)] = c;
    }
}

void FileParser::freeFileResources(InFileStream* filePtr)
{
    if(filePtr != 0)
    {
        filePtr->close();
        delete filePtr;
        filePtr = 0;    
    }
}


char FileParser::getDelimiter(string filename)
{
    ifstream in;
    int type = 0;
    char delim;

    in.open(filename.c_str(), ios::in);
    in.seekg(0, ios::beg);

    //look for CR or LF or CRLF (or LFCR)
    if (in.is_open()) {
        char c;
        while (in.get(c)) {
            if (c == CR)
                type |= 1;
            else if (c == LF)
                type |= 2;
            else if (type)
                break;
        }
    }
    in.close();

    switch (type) {
        case 1:
            //cout << "file is Mac System 9" << endl;
            delim = '\r';
            break;
        case 2:
            //cout << "file is UNIX" << endl;
            delim = '\n';
            break;
        case 3:
            //cout << "file is DOS" << endl;
            delim = '\n';
            break;
        default: //short or empty file
            //cout << "file is UNIX (default)" << endl;
            delim = '\n';
    }
    return delim;
}

}
