/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Changes:
 *
 * 10-02-07,Nigel Brown(EMBL): changed ifstream to InFileStream to handle
 * cross-platform end-of-lines.
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "MSFFileParser.h"

namespace clustalw
{

/**
 * MSFFileParser contructor sets up the chartab array.
 * @param filePath 
 * @return 
 */
MSFFileParser::MSFFileParser(string filePath)
{
    fileName = filePath; 
    fillCharTab();
}


    
    vector<Sequence> MSFFileParser::getSeqRange(int firstSeq, int no, string *offendingSeq)
{
    vector<Sequence> seqRangeVector;
    int i;

    for (i=0; i<no; i++)
    { 
        Sequence tempSeq = getSeq(firstSeq + i);
        if (parseExitCode!=OK) {
            seqRangeVector.clear();
            return seqRangeVector;
        }
        seqRangeVector.push_back(tempSeq);
    }
    return seqRangeVector;
}



/**
 * The function getSeq finds the sequence seqNum in the file and returns it.
 * @param seqNum The number of the sequence in the file to get.
 * @return A sequence object containing the seqNum'th sequence from the file.
 */
Sequence MSFFileParser::getSeq(int seqNum, string *offendingSeq)
{
    char _line[MAXLINE + 1];
    char _sname[MAXNAMES + 1];
    string characterSeq = "";
    string name = "";
    string title = "";
    string blank = "";
    
    _line[0] = EOS;
    int i, j, k;
    unsigned char c;

    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg);       
        
        for (i = 0;; i++)
        {
            if (!_fileIn->getline(_line, MAXLINE + 1))
            {
                _fileIn->close();
                return Sequence(blank, blank, blank);
            }
            // read the title
            if (utilityObject->lineType(_line, "//"))
            {
                break;
            }
            // lines...ignore
        }

        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            if (!utilityObject->blankLine(_line))
            {
                for (i = 1; i < seqNum; i++)
                {
                    _fileIn->getline(_line, MAXLINE + 1);
                }
                for (j = 0; j <= (int)strlen(_line); j++)
                {
                    if (_line[j] != ' ')
                    {
                        break;
                    }
                }
                for (k = j; k <= (int)strlen(_line); k++)
                {
                    if (_line[k] == ' ')
                    {
                        break;
                    }
                }
                
                // Get the name of the sequence
                strncpy(_sname, _line + j, utilityObject->MIN(MAXNAMES, k - j));
                _sname[utilityObject->MIN(MAXNAMES, k - j)] = EOS;
                utilityObject->rTrim(_sname);
                utilityObject->blankToUnderscore(_sname);
                name = string(_sname);

                for (i = k; i <= MAXLINE; i++)
                {
                    c = _line[i];
                    if (c == '.' || c == '~')
                    {
                        c = '-';
                    }
                    if (c == '*')
                    {
                        c = 'X';
                    }
                    if (c == '\n' || c == EOS)
                    {
                        break;
                    }
                    // EOL 
                    c = chartab[c];
                    if (c)
                    {
                        characterSeq += c;
                    }
                }

                for (i = 0;; i++)
                {
                    if (!_fileIn->getline(_line, MAXLINE + 1))
                    {
                        _fileIn->close();
                        return Sequence(characterSeq, name, title);
                    }
                    if (utilityObject->blankLine(_line))
                    {
                        break;
                    }
                }
            }
        }
        _fileIn->close();
        
        if ((int)characterSeq.length() > userParameters->getMaxAllowedSeqLength())
        {
            parseExitCode=SEQUENCETOOBIG;
            if (offendingSeq!=NULL)
                offendingSeq->assign(name);
            // return empty seq
            return Sequence(blank, blank, blank);
        }
        return Sequence(characterSeq, name, title);;
    }
    catch(...)
    {
        _fileIn->close();
        cerr << "An exception has occured in the function MSFFileParser::getSeq()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}

/**
 * The function countSeqs counts the number of sequences in the file.
 * @return The number of sequences in the file.
 */
int MSFFileParser::countSeqs()
{
    char _line[MAXLINE + 1];
    int _numSeqs;
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
    
        if(!_fileIn->is_open())
        {
            return 0; // No sequences found!
        }
    
        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            if (utilityObject->lineType(_line, "//"))
            {
                break;
            }
        }

        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            if (!utilityObject->blankLine(_line))
            {
                break;
            }
            // Look for next non- blank line
        } 
        _numSeqs = 1;

        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            if (utilityObject->blankLineNumericLabel(_line))
            {
                _fileIn->close();
                return _numSeqs;
            }
            _numSeqs++;
        }

        return 0; // if you got to here-funny format/no seqs.
    }
    catch(...)
    {
        _fileIn->close();
        cerr << "An exception has occured in the function MSFFileParser::countSeqs()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}

/**
 * There is no secondary structure information in MSF files. Set structPenalties to NONE.
 * @param gapPenaltyMask 
 * @param secStructMask 
 * @param secStructName 
 * @param structPenalties 
 * @param length 
 */
void MSFFileParser::getSecStructure(vector<char>& gapPenaltyMask, vector<char>& secStructMask,
                                    string& secStructName, int &structPenalties, int length)
{
    structPenalties = NONE;
}


}

