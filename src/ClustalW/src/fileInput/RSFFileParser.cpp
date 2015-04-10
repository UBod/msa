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
#include "RSFFileParser.h"

namespace clustalw
{

/**
 * Constructor sets up the chartab array.
 * @param filePath 
 */
RSFFileParser::RSFFileParser(string filePath)
{
    fileName = filePath; 
    fillCharTab();
}

    vector<Sequence> RSFFileParser::getSeqRange(int firstSeq, int no, string *offendingSeq)
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
 * get the sequence seqNum from the file.
 * @param seqNum The number of the sequence to get.
 * @return The sequence seqNum.
 */
    Sequence RSFFileParser::getSeq(int seqNum, string *offendingSeq)
{
    char _line[MAXLINE + 1];
    char _sname[MAXNAMES + 1];
    string characterSeq = "";
    string name = "";
    string title = "";
    string blank = "";
    _line[0] = EOS;
    
    int i;
    unsigned char c;
    int _currentSeqNum = 0; // Not at any sequence yet!
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg); // start at the beginning
        
        // Need to get the cursor to the begining of the correct sequence.
        // This will be the case when we get to the seqNum {
        while (_currentSeqNum != seqNum)
        {
            while(*_line != '{')
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
            // Get next line so that we are past the '{' line
            _fileIn->getline(_line, MAXLINE + 1);
        }

        while (!keyword(_line, "name"))
        {
            if (!_fileIn->getline(_line, MAXLINE + 1))
            {
                _fileIn->close();
                return Sequence(blank, blank, blank);
            }
        }
        for (i = 5; i <= (int)strlen(_line); i++)
        {
            if (_line[i] != ' ')
            {
                break;
            }
        }
        strncpy(_sname, _line + i, MAXNAMES); // remember entryname
        for (i = 0; i <= (int)strlen(_sname); i++)
        {
            if (_sname[i] == ' ')
            {
                _sname[i] = EOS;
                break;
            }
        }

        _sname[MAXNAMES] = EOS;
        utilityObject->rTrim(_sname);
        utilityObject->blankToUnderscore(_sname); // replace blanks with '_'
        name = string(_sname);


        while (!keyword(_line, "sequence"))
        {
            if (!_fileIn->getline(_line, MAXLINE + 1))
            {
                _fileIn->close();
                return Sequence(blank, blank, blank);
            }
        }
            
        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            for (i = 0; i <= MAXLINE; i++)
            {
                c = _line[i];
                if (c == EOS || c == '}')
                {
                    break;
                }
                 // EOL
                if (c == '.')
                {
                    characterSeq += '-';
                }
                c = chartab[c];
                if (c)
                {
                    characterSeq += c;
                }
            }
            if (c == '}')
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
        cerr << "There was an exception in the RSFFileParser::getSeq function.\n"
             << "Need to end program\n";
        throw 1;
    }            
}

/**
 * count the number of sequences in a GCG RSF alignment file 
 * @return The number of sequences in the file.
 */
int RSFFileParser::countSeqs()
{
    char _line[MAXLINE + 1];
    int numSeqs;

    try
    {
        numSeqs = 0;
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg); // start at the beginning
                
        if(!_fileIn->is_open())
        {
            return 0; // No sequences found!
        }
                        
        // skip the comments 
        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            // NOTE needed to change to -1 and -2 (it was -2 and -3)
            // This is because getline does not put the \n in!
            if (_line[strlen(_line) - 1] == '.' && _line[strlen(_line) - 2] == '.')
            {
                break;
            }
        }

        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            if (*_line == '{')
            {
                numSeqs++;
            }
        }
        _fileIn->close();
        return numSeqs;
    }
    catch(...)
    {
        _fileIn->close();
        cerr << "An exception has occured in the function RSFFileParser::countSeqs()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }    
}

/**
 * Get the secondary structure information from the file.
 * @param gapPenaltyMask 
 * @param secStructMask 
 * @param secStructName 
 * @param structPenalties 
 * @param length 
 */
void RSFFileParser::getSecStructure(vector<char>& gapPenaltyMask, vector<char>& secStructMask,
                     string& secStructName, int &structPenalties, int length)
{
    bool guigetss = false;
    if(userParameters->getProfileNum() == 1 && userParameters->getStructPenalties1())
         guigetss = true;
    if(userParameters->getProfileNum() == 2 && userParameters->getStructPenalties2())
         guigetss = true;

    char _title[MAXLINE + 1];
    char _line[MAXLINE + 1];
    char _lin2[MAXLINE + 1];
    char _sname[MAXNAMES + 1];
    int i;
    _line[0] = EOS;
    
    try
    {
        secStructMask.clear();
        secStructMask.assign(length, '.');        
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg); // Need to start at begining
                
        // skip the comments 
        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            if (_line[strlen(_line) - 1] == '.' && _line[strlen(_line) - 2] == '.')
            {
                break;
            }
        }

        // find the start of the sequence entry 
        for (;;)
        {
            while (_fileIn->getline(_line, MAXLINE + 1))
                if (*_line == '{')
                {
                    break;
                }

            while (!keyword(_line, "name"))
            {
                if (!_fileIn->getline(_line, MAXLINE + 1))
                {
                    _fileIn->close();
                    return;
                }
            }
        
            for (i = 5; i <= (int)strlen(_line); i++)
            {
                if (_line[i] != ' ')
                {
                    break;
                }
            }
            strncpy(_sname, _line + i, MAXNAMES); // remember entryname
            for (i = 0; i <= (int)strlen(_sname); i++)
            {
                if (_sname[i] == ' ')
                {
                    _sname[i] = EOS;
                    break;
                }
            }
            _sname[MAXNAMES] = EOS;
            utilityObject->rTrim(_sname);
            utilityObject->blankToUnderscore(_sname); // replace blanks with '_'

            // look for secondary structure feature table / gap penalty mask
            while (_fileIn->getline(_line, MAXLINE + 1))
            {
                if (keyword(_line, "feature"))
                {
                    if (userParameters->getInteractive() && !userParameters->getGui())
                    {
                        strcpy(_title, "Found secondary structure in alignment file: ");
                        strcat(_title, _sname);
                        (*_lin2) = utilityObject->promptForYesNo(_title,
                            "Use it to set local gap penalties ");
                    }
                    else
                    {
                        (*_lin2) = 'y';
                    }
                    if (guigetss || ((*_lin2 != 'n') && (*_lin2 != 'N')))
                    {
                        structPenalties = SECST;
                        secStructMask.assign(length, '.');
                        do
                        {
                            if (keyword(_line, "feature"))
                            {
                                getRSFFeature(&_line[7], secStructMask, length);
                            }
                            _fileIn->getline(_line, MAXLINE + 1);
                        }
                        while (!keyword(_line, "sequence"));
                    }
                    else
                    {
                        do
                        {
                            _fileIn->getline(_line, MAXLINE + 1);
                        }
                        while (!keyword(_line, "sequence"));
                    }
                    secStructName = string(_sname);
                }
                else if (keyword(_line, "sequence"))
                {
                    break;
                }

                if (structPenalties != NONE)
                {
                    break;
                }
            }
        }
        _fileIn->close();
    }
    catch(...)
    {
        _fileIn->close();
        cerr << "An exception has occured in the function RSFFileParser::getSecStructure()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }    
}

/**
 * get a feature from the file. Called by getSecStructure
 * @param line 
 * @param secStructMask 
 * @param length 
 */
void RSFFileParser::getRSFFeature(char* line, vector<char>& secStructMask, int length)
{
    char c, s;
    char str1[MAXLINE + 1], str2[MAXLINE + 1], feature[MAXLINE + 1];
    int i, tmp, startPos, endPos;

    try
    {
        if (sscanf(line, "%d%d%d%s%s%s", &startPos, &endPos, &tmp, str1, str2,
            feature) != 6)
        {
            return;
        }

        if (strcmp(feature, "HELIX") == 0)
        {
            c = 'A';
            s = '$';
        }
        else if (strcmp(feature, "STRAND") == 0)
        {
            c = 'B';
            s = '%';
        }
        else
        {
            return ;
        }

        if (startPos >= length || endPos >= length)
        {
            return ;
        }
        secStructMask[startPos - 1] = s;
        for (i = startPos; i < endPos - 1; i++)
        {
            secStructMask[i] = c;
        }
        secStructMask[endPos - 1] = s;
    }
    catch(...)
    {
        cerr << "An exception has occured in the function RSFFileParser::getRSFFeature()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}

/**
 * keyword checks if code is on the line!
 * @param line 
 * @param code 
 * @return 
 */
bool RSFFileParser::keyword(char *line, const char *code)
{
    int i;
    char key[MAXLINE];

    for (i = 0; !isspace(line[i]) && line[i] != EOS; i++)
    {
        key[i] = line[i];
    }
    key[i] = EOS;
    return (strcmp(key, code) == 0);
}

}


