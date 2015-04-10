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
#include "EMBLFileParser.h"

namespace clustalw
{

/*
 * constructor sets up chartab array.
 */
EMBLFileParser::EMBLFileParser(string filePath)
{
    fileName = filePath; 
    fillCharTab();
}

/*
 * dont need to destruct anything.
 */
EMBLFileParser::~EMBLFileParser()
{

}


/*
 * get range of sequences
 */
    vector<Sequence> EMBLFileParser::getSeqRange(int firstSeq, int no, string *offendingSeq)
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


/*
 * get the sequence seqNum in the file.
 */
    Sequence EMBLFileParser::getSeq(int seqNum, string *offendingSeq)
{
    char _line[MAXLINE + 1];
    //char _tseq[MAXLINE + 1];
    char _sname[MAXNAMES + 1];
    //char _title[MAXTITLES + 1];
    string characterSeq = "";
    string name = "";
    string title = "";
    
    _line[0] = EOS;
    int i;
    //int j;
    bool gotSeq = false;
    unsigned char c;
    int _currentSeqNum = 0;
    string blank = "";

    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg); // start at the beginning

        // Read in lines until we get to the begining of sequence seqNum.
        while (_currentSeqNum != seqNum)
        {
            while(!utilityObject->lineType(_line, "ID"))
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
        utilityObject->blankToUnderscore(_sname);
        name = string(_sname);
        // Andreas Wilm (UCD): why cout here? cout << name << "\n"; 
        
        while (!utilityObject->lineType(_line, "SQ"))
        {
            if(!_fileIn->getline(_line, MAXLINE + 1)) // If we cannot get anymore!
            {
                _fileIn->close();
                // FIXME AW: why return with a name but otherwise empty seq?
                return Sequence(blank, name, blank);
            }
        }

        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            if (gotSeq && utilityObject->blankLine(_line))
            {
                break;
            }
            
            // NOTE I changed this to -1 and -2 because the getline doesnt return the \n
            if (strlen(_line) > 2 && _line[strlen(_line) - 1] == '.' &&
                    _line[strlen(_line) - 2] == '.')
            {
                continue;
            }

            for (i = 0; i <= MAXLINE; i++)
            {
                c = _line[i];
                if (c == '\n' || c == EOS || c == '/')
                {
                    break;
                }
                // EOL     
                c = chartab[c];
                if (c)
                {
                    gotSeq = true;
                    characterSeq += c;
                }
            }
            if (c == '/')
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
        cerr << "There was an exception in the EMBLFileParser::getSeq function.\n"
             << "Need to end program\n";
        throw 1;
    }            
}

/*
 * count the number of sequences in the file and return the number
 */
int EMBLFileParser::countSeqs()
{
    char line[MAXLINE + 1];
    // char c;
    line[0] = EOS;
    int numSeqs;
    // int i;
    //bool seqOk;
    numSeqs = 0;
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
    
        if(!_fileIn->is_open())
        {
            return 0; // No sequences found!
        }       
                 
        while (_fileIn->getline(line, MAXLINE + 1))
        {
            if (utilityObject->lineType(line, "ID"))
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
        cerr << "An exception has occured in the function EMBLFileParser::countSeqs()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}

/*
 * get secondary structure information from the file.
 */
void EMBLFileParser::getSecStructure(vector<char>& gapPenaltyMask, vector<char>& secStructMask, string& secStructName, int &structPenalties, int length)
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
    char _feature[MAXLINE + 1];
    int i;
    _line[0] = '\0';
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg);

        // clear out the masks
        gapPenaltyMask.clear();
        secStructMask.clear();
            
        // find the start of the sequence entry 
        for (;;)
        {
            while (!utilityObject->lineType(_line, "ID"))
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
            utilityObject->blankToUnderscore(_sname);
            
            // look for secondary structure feature table / gap penalty mask 
            while (_fileIn->getline(_line, MAXLINE + 1))
            {
                if (utilityObject->lineType(_line, "FT"))
                {
                    sscanf(_line + 2, "%s", _feature);
                    if (strcmp(_feature, "HELIX") == 0 || strcmp(_feature, "STRAND") == 0)
                    {
                        if (userParameters->getInteractive() && !userParameters->getGui())
                        {
                            strcpy(_title,
                                "Found secondary structure in alignment file: ");
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
                            for (i = 0; i < length; i++)
                            {
                                secStructMask.push_back('.');
                            }
                            do
                            {
                                getSwissFeature(&_line[2], secStructMask, length);
                                _fileIn->getline(_line, MAXLINE + 1);
                            }
                            while (utilityObject->lineType(_line, "FT"));
                        }
                        else
                        {
                            do
                            {
                                _fileIn->getline(_line, MAXLINE + 1);
                            }
                            while (utilityObject->lineType(_line, "FT"));
                        }
                        secStructName = string(_sname);
                    }
                }
                else if (utilityObject->lineType(_line, "GM"))
                {
                    if (userParameters->getInteractive())
                    {
                        strcpy(_title, "Found gap penalty mask in alignment file: ");
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
                        structPenalties = GMASK;
                        for (i = 0; i < length; i++)
                        {
                            gapPenaltyMask.push_back('1');
                        }
                        do
                        {
                            getSwissMask(&_line[2], gapPenaltyMask, length);
                            _fileIn->getline(_line, MAXLINE + 1);
                        }
                        while (utilityObject->lineType(_line, "GM"));
                    }
                    else
                    {
                        do
                        {
                            _fileIn->getline(_line, MAXLINE + 1);
                        }
                        while (utilityObject->lineType(_line, "GM"));
                    }
                    secStructName = string(_sname);
                }
                if (utilityObject->lineType(_line, "SQ"))
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
        cerr << "An exception has occured in the function EMBLFileParser::getSecStructure()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
    
}

/*
 * get the sec structure mask
 */
void EMBLFileParser::getSwissFeature(char* line, vector<char>& secStructMask, int length)
{
    char c, s, feature[MAXLINE + 1];
    int i, startPos, endPos;

    try
    {
        if (sscanf(line, "%s%d%d", feature, &startPos, &endPos) != 3)
        {
            return ;
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
        cerr << "An exception has occured in the function EMBLFileParser::getSwissFeature()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}

/*
 * get the gap penalty mask.
 */
void EMBLFileParser::getSwissMask(char* line, vector<char>& gapPenaltyMask, int length)
{
    int i, value, startPos, endPos;

    try
    {
        if (sscanf(line, "%d%d%d", &value, &startPos, &endPos) != 3)
        {
            return;
        }

        if (value < 1 || value > 9)
        {
            return ;
        }

        if (startPos >= length || endPos >= length)
        {
            return ;
        }
        for (i = startPos - 1; i < endPos; i++)
        {
            gapPenaltyMask[i] = value + '0';
        }
    }
    catch(...)
    {
        cerr << "An exception has occured in the function EMBLFileParser::getSwissMask()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }    
}


}
