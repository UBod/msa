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
#include "PIRFileParser.h"

namespace clustalw
{

/**
 * PIRFileParser contructor sets up the chartab array.
 * @param filePath 
 */
PIRFileParser::PIRFileParser(string filePath)
{
    fileName = filePath; 
    fillCharTab();
}


/*
 * get range of sequences
 */
    vector<Sequence> PIRFileParser::getSeqRange(int firstSeq, int no, string *offendingSeq)
{
    vector<Sequence> seqRangeVector;
    int i;

    for (i=0; i<no; i++)
    { 
        Sequence tempSeq = getSeq(firstSeq + i, offendingSeq);
        if (parseExitCode!=OK) {
            seqRangeVector.clear();
            return seqRangeVector;
        }
        seqRangeVector.push_back(tempSeq);
    }
    return seqRangeVector;
}



/**
 * The function getSeq finds the sequence 'seqNum' in the file and returns it.
 * @param seqNum The number of the sequence to get from the file.
 * @return The 'seqNum' sequence from the file.
 */
    Sequence PIRFileParser::getSeq(int seqNum, string *offendingSeq)
{
    char _line[MAXLINE + 1];
    char _sname[MAXNAMES + 1];
    char _title[MAXTITLES + 1];
    string characterSeq = "";
    string name = "";
    string title = "";
    string blank = "";
    
    _line[0] = EOS;
    int i;
    unsigned char c;
    int _currentSeqNum = 0;
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg);
        
        // Read in lines until we get to the begining of sequence seqNum.
        while (_currentSeqNum != seqNum)
        {
            while(*_line != '>')
            {
                if(!_fileIn->getline(_line, MAXLINE + 1)) // If we cannot get anymore!
                {
                    _fileIn->close();
                    return Sequence(blank, blank, blank);
                }
            }
            ++_currentSeqNum;
            if(_currentSeqNum == seqNum) // Found the sequence
            {
                break;
            }
            // Get next line so that we are past the '>' line
            _fileIn->getline(_line, MAXLINE + 1);
        }        
        
        // line contains the name of the sequence
        for (i = 4; i <= (int)strlen(_line); i++)
        {
            if (_line[i] != ' ')
            {
                break;
            }
        }
        
        strncpy(_sname, _line + i, MAXNAMES); // remember entryname 
        _sname[MAXNAMES] = EOS;
        utilityObject->rTrim(_sname);
        utilityObject->blankToUnderscore(_sname); // replace blanks with '_'
        name = string(_sname);
        
        _fileIn->getline(_line, MAXLINE + 1);
        strncpy(_title, _line, MAXTITLES);
        _title[MAXTITLES] = EOS;
        i = strlen(_title);
        if (_title[i - 1] == '\n')
        {
            _title[i - 1] = EOS;
        }
        title = string(_title);
        
        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            for (i = 0; i <= MAXLINE; i++)
            {
                c = _line[i];
                if (c == '\n' || c == EOS || c == '*')
                {
                    break;
                }

                c = chartab[c];
                if (c)
                {
                    characterSeq += c;
                }
            }
            if (c == '*')
            {
                break;
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
        return Sequence(characterSeq, name, title);
    }
    catch(...)
    {
        _fileIn->close();
        cerr << "There was an exception in the PIRFileParser::getSeq function.\n"
             << "Need to end program\n";
        throw 1;
    }
}

/**
 * The function countSeqs finds the number of sequences in the file and returns it.
 * @return The number of sequences in the file.
 */
int PIRFileParser::countSeqs()
{
    char line[MAXLINE + 1], c;
    line[0] = EOS;
    int numSeqs, i;
    bool seqOk;
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
    
        if(!_fileIn->is_open())
        {
            return 0; // No sequences found!
        }
    
        // Get to begining of sequences!
        while (_fileIn->getline(line, MAXLINE + 1))
        {
            if (!utilityObject->blankLine(line))
            {
                break;
            }
        }
    
        // Now check the 1st sequence to make sure it ends with *
        seqOk = false;
        while (_fileIn->getline(line, MAXLINE + 1))
        {
             // Look for end of first seq
            if (*line == '>')
            {
                break;
            }
            for (i = 0; seqOk == false; i++)
            {
                c = line[i];
                if (c == '*')
                {
                    seqOk = true; // ok - end of sequence found
                    break;
                } // EOL 
                if (c == '\n' || c == EOS)
                {
                    break;
                }
                // EOL
            }
            if (seqOk == true)
            {
                break;
            }
        }
        if (seqOk == false)
        {
            _fileIn->close();
            utilityObject->error("PIR format sequence end marker '*'\nmissing for one or more sequences.\n");     
            return 0; // funny format
        }

        numSeqs = 1;
    
        while (_fileIn->getline(line, MAXLINE + 1))
        {
            if (*line == '>')
            {
                // Look for start of next seq 
                seqOk = false;
                while (_fileIn->getline(line, MAXLINE + 1))
                {
                    // Look for end of seq
                    if (*line == '>')
                    {
                        _fileIn->close();
                        utilityObject->error("PIR format sequence end marker '*'\nmissing for one or more sequences.\n");     
                        return 0; // funny format
                    }
                    for (i = 0; seqOk == false; i++)
                    {
                        c = line[i];
                        if (c == '*')
                        {
                            seqOk = true; // ok - sequence found
                            break;
                        }
                        if (c == '\n' || c == EOS)
                        {
                            break;
                        }
                    }
                    if (seqOk == true)
                    {
                        numSeqs++;
                        break;
                    }
                }
            }
        }
    
        _fileIn->close();
    
        return numSeqs;
    }
    catch(...)
    {
        _fileIn->close();
        cerr << "An exception has occured in the function PIRFileParser::countSeqs()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}

/**
 * There is no secondary structure information in PIR files!
 * @param gapPenaltyMask 
 * @param secStructMask 
 * @param secStructName 
 * @param structPenalties 
 * @param length 
 */
void PIRFileParser::getSecStructure(vector<char>& gapPenaltyMask, vector<char>& secStructMask,
                                    string& secStructName, int &structPenalties, int length)
{
    structPenalties = NONE;
}

}

