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
#include "GDEFileParser.h"

namespace clustalw
{

/*
 * Constructor sets up the chartab array.
 *
 */
GDEFileParser::GDEFileParser(string filePath)
{
    fileName = filePath; 
    fillCharTab();
}

/*
 * Nothing to do in destruction of object.
 */
GDEFileParser::~GDEFileParser()
{

}


    vector<Sequence> GDEFileParser::getSeqRange(int firstSeq, int no, string *offendingSeq)
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
 * The getSeq function is used to get sequence number seqNum from the file.
 */
    Sequence GDEFileParser::getSeq(int seqNum, string *offendingSeq)
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
        _fileIn->seekg(0, std::ios::beg);
                
        bool dnaFlagSet = userParameters->getDNAFlag();
        while (_currentSeqNum != seqNum)
        {
            while((*_line != '#' && dnaFlagSet) ||
                  (*_line != '%' && !dnaFlagSet))
            {
                if(!_fileIn->getline(_line, MAXLINE + 1)) 
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
            // Get next line so that we are past the '#' or '%' line
            _fileIn->getline(_line, MAXLINE + 1);  //nige
        }
        
        for (i = 1; i <= MAXNAMES; i++)
        {
            if (_line[i] == '(' || _line[i] == '\n' || _line[i] == '\r')
            {
                break;
            }
            _sname[i - 1] = _line[i];
        }
        i--;
        _sname[i] = EOS;

        for (i--; i > 0; i--)
        {
            if (isspace(_sname[i]))
            {
                _sname[i] = EOS;
            }
            else
            {
                break;
            }
        }
        utilityObject->blankToUnderscore(_sname);
        name = string(_sname);
        title = "";
        
        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            if (*_line == '%' ||  *_line == '#' ||  *_line == '"')
            {
               break;
            }
            for (i = 0; i <= MAXLINE; i++)
            {
                c = _line[i];
                if (c == '\n' || c == EOS)
                {
                    break;
                }

                c = chartab[c];
                if (c)
                {
                    characterSeq += c;
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
        return Sequence(characterSeq, name, title);
    }
    catch(...)
    {
        _fileIn->close();
        cerr << "There was an exception in the GDEnFileParser::getSeq function.\n"
             << "Need to end program\n";
        throw 1;
    }

}

/*
 * The countSeqs function returns the number of sequences in the file. 
 */
int GDEFileParser::countSeqs()
{
    char line[MAXLINE + 1];
    int _nseqs = 0;
    
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
            if ((*line == '%') && (userParameters->getDNAFlag() == false))
            {
                _nseqs++;
            }
            else if ((*line == '#') && (userParameters->getDNAFlag() == true))
            {
                _nseqs++;
            }
        }
        _fileIn->close();

        return _nseqs;
    }
    catch(...)
    {
        _fileIn->close();
        cerr << "An exception has occured in the function GDEFileParser::countSeqs()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }    
}

/*
 * getSecStructure gets the secondary structure from the file.
 */
void GDEFileParser::getSecStructure(vector<char>& gapPenaltyMask, vector<char>& secStructMask,
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
    int i, len, offset = 0;
    unsigned char c;
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg);
    
        // NOTE I think I should empty the masks before pushing onto them!
        gapPenaltyMask.clear();
        secStructMask.clear();
    
        for (;;)
        {
            _line[0] = '\0';
            // search for the next comment line
            while (*_line != '"')
            {
                if (!_fileIn->getline(_line, MAXLINE + 1))
                {
                    _fileIn->close();
                    return;
                }
            }

            // is it a secondary structure entry? 
            if (strncmp(&_line[1], "SS_", 3) == 0)
            {
                for (i = 1; i <= MAXNAMES - 3; i++)
                {
                    if (_line[i + 3] == '(' || _line[i + 3] == '\n' || _line[i + 3] == '\r')
                    {
                        break;
                    }
                    _sname[i - 1] = _line[i + 3];
                }
                i--;
                _sname[i] = EOS;
            
                // NOTE NOTE NOTE
                // Is it possible for this to be executed????????????????
                // if _line contains ( then we break and dont put it into _sname
                // So how can sname have it???????
                if (_sname[i - 1] == '(')
                {
                    sscanf(&_line[i + 3], "%d", &offset);
                }
                else
                {
                    offset = 0;
                }
                for (i--; i > 0; i--)
                {
                    if (isspace(_sname[i]))
                    {
                        _sname[i] = EOS;
                    }
                    else
                    {
                        break;
                    }
                }

                utilityObject->blankToUnderscore(_sname);
                secStructName = string(_sname);
            
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
                    for (i = 0; i < length; i++)
                    {
                        secStructMask.push_back('.');
                    }
                    len = 0;
                    while (_fileIn->getline(_line, MAXLINE + 1))
                    {
                        if (*_line == '%' ||  *_line == '#' ||  *_line == '"')
                        {
                            break;
                        }
                        for (i = offset; i < length; i++)
                        {
                            c = _line[i];
                            if (c == '\n' || c == EOS)
                            {
                                break;
                            }
                            // EOL
                            secStructMask[len++] = c;
                        }
                        if (len >= length) // NOTE i put in >=
                        {
                            break;
                        }
                    }
                }
            }
        
            // or is it a gap penalty mask entry?
            else if (strncmp(&_line[1], "GM_", 3) == 0)
            {
                for (i = 1; i <= MAXNAMES - 3; i++)
                {
                    if (_line[i + 3] == '(' || _line[i + 3] == '\n')
                    {
                        break;
                    }
                    _sname[i - 1] = _line[i + 3];
                }
                i--;
                _sname[i] = EOS;
            
                // NOTE NOTE
                // Again I dont think it is possible for _sname to have ( !!!!
                if (_sname[i - 1] == '(')
                {
                    sscanf(&_line[i + 3], "%d", &offset);
                }
                else
                {
                    offset = 0;
                }
                for (i--; i > 0; i--)
                {
                    if (isspace(_sname[i]))
                    {
                        _sname[i] = EOS;
                    }
                    else
                    {
                        break;
                    }
                }
            
                utilityObject->blankToUnderscore(_sname);
                secStructName = string(_sname);

                if (userParameters->getInteractive() && !userParameters->getGui())
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
                    len = 0;
                    while (_fileIn->getline(_line, MAXLINE + 1))
                    {
                        if (*_line == '%' ||  *_line == '#' ||  *_line == '"')
                        {
                            break;
                        }
                        for (i = offset; i < length; i++)
                        {
                            c = _line[i];
                            if (c == '\n' || c == EOS)
                            {
                                break;
                            }
                            // EOL
                            gapPenaltyMask[len++] = c;
                        }
                        if (len >= length) // NOTE I put in >=
                        {
                            break;
                        }
                    }
                }
            }
            if (structPenalties != NONE)
            {
                break;
            }
        }
        _fileIn->close();
    }
    catch(...)
    {
        _fileIn->close();
        cerr << "An exception has occured in the function GDEFileParser::getSecStructure()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}

}

