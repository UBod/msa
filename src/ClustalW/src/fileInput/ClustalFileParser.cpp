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
 *
 * 27-4-2007, Mark Larkin (UCD): Made 2 small changes to getSecStructure function. There 
 * was a problem with the secondary structure info in windows.
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "ClustalFileParser.h"

namespace clustalw
{

ClustalFileParser::ClustalFileParser(string filePath)
{
    fileName = filePath; 
    fillCharTab();
    _fileIn = 0;
}
        
ClustalFileParser::~ClustalFileParser()
{
    
}


    vector<Sequence> ClustalFileParser::getSeqRange(int firstSeq, int no, string *offendingSeq)
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


Sequence ClustalFileParser::getSeq(int seqNum, string *offendingSeq)
{
    char line[MAXLINE + 1];
    line[0] = EOS;
    char tseq[MAXLINE + 1];
    tseq[0] = '\0';
    char sname[MAXNAMES + 150];

    for(int i = 1; i < MAXNAMES + 1; i++)
    {
        line[i] = tseq[i] = sname[i] = '0';
    }
    string characterSeq = "";
    string name = "";
    string title = ""; // Nothing happens it here!!!

    int i, j;
    unsigned char c;
    string blank = "";
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg); // start at the beginning
             
        _fileIn->getline(line, MAXLINE + 1); // read the title line...ignore it
    
        while (_fileIn->getline(line, MAXLINE + 1))  //nige
        {
            if (!clustalBlankline(line))
            {

                for (i = 1; i < seqNum; i++)
                {
                    _fileIn->getline(line, MAXLINE + 1);  //nige
                }
                for (j = 0; j <= (int)strlen(line); j++)
                    if (line[j] != ' ')
                    {
                        break;
                    }
                string _tempStr = string(line);

                sscanf(_tempStr.c_str(), "%s%s", sname, tseq);
                for (j = 0; j < MAXNAMES; j++)
                {
                    if (sname[j] == ' ')
                    {
                        break;
                    }
                }
                sname[j] = EOS;
                utilityObject->rTrim(sname);
                utilityObject->blankToUnderscore(sname); // replace blanks with '_'
                name = string(sname);           
               
                for (i = 0; i <= MAXLINE; i++)
                {
                    c = tseq[i];
                    if (isspace(c) || c == EOS)
                    {
                        break;
                    }
                    // EOL 
                    c = chartab[c];
                    if (c)
                    {
                        characterSeq += c; // Add the character to the sequence
                    }
                }

                for (i = 0;; i++)
                {
                    if (!_fileIn->getline(line, MAXLINE + 1)) // If we cant get another line!
                    {
                        freeFileResources(_fileIn);
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
                    if (clustalBlankline(line))
                    {
                        break;
                    }
                }
            }
        }
        freeFileResources(_fileIn);

	// getSecStructure(vector<char>& gapPenaltyMask, 
	//      vector<char>& secStructMask, string& secStructName, int &structPenalties, int length)

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
        freeFileResources(_fileIn);
        cerr << "An exception has occured in the function ClustalFileParser::getSeq()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}
        
/*
 * The function countSeqs tells us how many sequences are in a clustal format file.
 * Need to check if the file is open!
 */
int ClustalFileParser::countSeqs()
{
    char line[MAXLINE + 1];
    int _nseqs;
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
    
        if(!_fileIn->is_open())
        {
            freeFileResources(_fileIn);            
            return 0; // No sequences found!
        }


        while (_fileIn->getline(line, MAXLINE + 1))
        {
            if (!utilityObject->blankLine(line))
            {
                break;
            }
        }
    
        // This gets us to the begining of the sequence lines!
        while (_fileIn->getline(line, MAXLINE + 1))
        {
            if (!clustalBlankline(line))
            {
                break;
            }
        }
        _nseqs = 1;

        while (_fileIn->getline(line, MAXLINE + 1))
        {
            if (clustalBlankline(line))
            {
                freeFileResources(_fileIn);
                return _nseqs;
            }
            _nseqs++;
        }
        freeFileResources(_fileIn);
        return (int)0; // if you got to here-funny format/no seqs.
    }
    catch(...)
    {
        freeFileResources(_fileIn);
        cerr << "An exception has occured in the function ClustalFileParser::countSeqs()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }    
}

/*
 * This function is to get the secondary structure for a Clustal format file.
 * I am aware that I am using some C and some C++ here, and this may seem like
 * bad style, but I think it is better to use the C functions for processing the
 * strings as they are already working.
 */        
void ClustalFileParser::getSecStructure(vector<char>& gapPenaltyMask, 
      vector<char>& secStructMask, string& secStructName, int &structPenalties, int length)
{
    bool guigetss = false;
    if(userParameters->getProfileNum() == 1 && userParameters->getStructPenalties1())
         guigetss = true;
    if(userParameters->getProfileNum() == 2 && userParameters->getStructPenalties2())
         guigetss = true;

    char title[MAXLINE + 1];
    title[0] = '\0';
    char line[MAXLINE + 1];
    line[0] = '\0';
    char lin2[MAXLINE + 1];
    lin2[0] = '\0';
    char tseq[MAXLINE + 1];
    tseq[0] = '\0';
    char sname[MAXNAMES + 1];
    sname[0] = '\0';
    
    for(int i = 1; i < MAXNAMES + 1; i++)
    {
        title[i] = line[i] = lin2[i] = tseq[i] = sname[i] ='0';
    }
    
    int i, j, len, ix, struct_index = 0;
    char c;
    
    try
    {
        _fileIn = new InFileStream;  //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg);

        // NOTE clear out the masks
        gapPenaltyMask.clear();
        secStructMask.clear();
        
        len = 0; // initialise length to zero 
    
        if (!_fileIn->getline(line, MAXLINE + 1))
        {
            freeFileResources(_fileIn);
            return ;
        }
        // read the title line...ignore it 

        if (!_fileIn->getline(line, MAXLINE + 1))
        {
            freeFileResources(_fileIn);
            return ;
        }
        // read the next line...
        // skip any blank lines 
        for (;;)
        {
            if (!_fileIn->getline(line, MAXLINE + 1))
            {
                freeFileResources(_fileIn);   
                return ;
            }
            if (!utilityObject->blankLine(line))
            {
                break;
            }
        }

        // look for structure table lines
        ix =  -1;
        for (;;)
        {
            if (line[0] != '!')
            {
                break;
            }
            if (strncmp(line, "!SS", 3) == 0)
            {
                ix++;
                sscanf(line + 4, "%s%s", sname, tseq);
                for (j = 0; j < MAXNAMES; j++)
                {
                    if (sname[j] == ' ')
                    {
                        break;
                    }
                }
                sname[j] = EOS;
                utilityObject->rTrim(sname);
                utilityObject->blankToUnderscore(sname);
            
                if (userParameters->getInteractive() && !userParameters->getGui())
                {
                    strcpy(title, "Found secondary structure in alignment file: ");
                    strcat(title, sname);
                    (*lin2) = utilityObject->promptForYesNo(title,
                        "Use it to set local gap penalties ");
                }
                else
                {
                    (*lin2) = 'y';
                }
                if (guigetss || ((*lin2 != 'n') && (*lin2 != 'N')))
                {
                    structPenalties = SECST;
                    struct_index = ix;
                    for (i = 0; i < length; i++)
                    {
                        secStructMask.push_back('.');
                        gapPenaltyMask.push_back('.');
                    }
                
                    secStructName = string(sname);
                
                    for (i = 0; len < length; i++)
                    {
                        c = tseq[i];
                        if (c == '\n' || c == EOS)
                        {
                            break;
                        }
                        // EOL
                        if (!isspace(c))
                        {
                            if(len < (int)secStructMask.size())
                            {
                                secStructMask[len++] = c; // NOTE array notation = BAD
                            }
                        }
                    }
                }
            }
            else if (strncmp(line, "!GM", 3) == 0)
            {
                ix++;
                sscanf(line + 4, "%s%s", sname, tseq);
                for (j = 0; j < MAXNAMES; j++)
                {
                    if (sname[j] == ' ')
                    {
                        break;
                    }
                }
                sname[j] = EOS;
                utilityObject->rTrim(sname);
                utilityObject->blankToUnderscore(sname);
            
                if (userParameters->getInteractive() && !userParameters->getGui())
                {
                    strcpy(title, "Found gap penalty mask in alignment file: ");
                    strcat(title, sname);
                    (*lin2) = utilityObject->promptForYesNo(title,
                        "Use it to set local gap penalties ");
                }
                else
                {
                    (*lin2) = 'y';
                }
                if (guigetss || ((*lin2 != 'n') && (*lin2 != 'N')))
                {
                    structPenalties = GMASK;
                    struct_index = ix;
                    for (i = 0; i < length; i++)
                    {
                        gapPenaltyMask.push_back('1');
                    }
                
                    secStructName = string(sname);
                
                    for (i = 0; len < length; i++)
                    {
                        c = tseq[i];
                        if (c == '\n' || c == EOS)
                        {
                            break;
                        }
                        // EOL
                        if (!isspace(c))
                        {
                            if(len < (int)gapPenaltyMask.size())
                            {                            
                                gapPenaltyMask[len++] = c;
                            }
                        }
                    }
                }
            }
            if (structPenalties != NONE)
            {
                break;
            }
            if (!_fileIn->getline(line, MAXLINE + 1))
            {
                freeFileResources(_fileIn);   
                return ;
            }
        }
        if (structPenalties == NONE)
        {
            freeFileResources(_fileIn);    
            return ;
        }

        // skip any more comment lines
        while (line[0] == '!')
        {
            if (!_fileIn->getline(line, MAXLINE + 1))
            {
                freeFileResources(_fileIn);   
                return ;
            }
        }

        // skip the sequence lines and any comments after the alignment 
        for (;;)
        {
            if (isspace(line[0]) || line[0] == '\0') // Mark change 27-4-2007
            {
                break;
            }
            if (!_fileIn->getline(line, MAXLINE + 1))
            {
                freeFileResources(_fileIn);   
                return ;
            }
        }


        // read the rest of the alignment

        for (;;)
        {
            // skip any blank lines
            for (;;)
            {
                if (!utilityObject->blankLine(line))
                {
                    break;
                }
                if (!_fileIn->getline(line, MAXLINE + 1))
                {
                    freeFileResources(_fileIn);    
                    return ;
                }
            }
            // get structure table line 
            for (ix = 0; ix < struct_index; ix++)
            {
                if (line[0] != '!')
                {
                    if (structPenalties == SECST)
                    {
                        utilityObject->error("bad secondary structure format\n");
                    }
                    else
                    {
                        utilityObject->error("bad gap penalty mask format\n");
                    }
                    structPenalties = NONE;
                    freeFileResources(_fileIn);    
                    return ;
                }
                if (!_fileIn->getline(line, MAXLINE + 1))
                {
                    freeFileResources(_fileIn);    
                    return ;
                }
            }
            if (structPenalties == SECST)
            {
                if (strncmp(line, "!SS", 3) != 0)
                {
                    utilityObject->error("bad secondary structure format\n");
                    structPenalties = NONE;
                    freeFileResources(_fileIn);    
                    return ;
                }
                sscanf(line + 4, "%s%s", sname, tseq);
                for (i = 0; len < length; i++)
                {
                    c = tseq[i];
                    if (c == '\n' || c == EOS)
                    {
                        break;
                    }
                    // EOL
                    if (!isspace(c))
                    {
                        secStructMask[len++] = c;
                    }
                }
            }
            else if (structPenalties == GMASK)
            {
                if (strncmp(line, "!GM", 3) != 0)
                {
                    utilityObject->error("bad gap penalty mask format\n");
                    structPenalties = NONE;
                    freeFileResources(_fileIn);    
                    return ;
                }
                sscanf(line + 4, "%s%s", sname, tseq);
                for (i = 0; len < length; i++)
                {
                    c = tseq[i];
                    if (c == '\n' || c == EOS)
                {
                    break;
                }
                // EOL
                if (!isspace(c))
                {
                    gapPenaltyMask[len++] = c;
                }
            }
        }

        // skip any more comment lines
        while (line[0] == '!')
        {
            if (!_fileIn->getline(line, MAXLINE + 1))
            {
                freeFileResources(_fileIn);   
                return ;
            }
        }

            // skip the sequence lines
            for (;;)
            {
                if (isspace(line[0]) || line[0] == '\0') // Mark change 27-4-2007
                {
                    break;
                }
                if (!_fileIn->getline(line, MAXLINE + 1))
                {
                    freeFileResources(_fileIn);
                    return ;
                }
            }
        }
        freeFileResources(_fileIn);
    }
    catch(...)
    {
        freeFileResources(_fileIn);
        cerr << "An exception has occured in the function ClustalFileParser::getSecStructure()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}

bool ClustalFileParser::clustalBlankline(char* line)
{
    int i;

    if (line[0] == '!')
    {
        return true;
    }

    for (i = 0; line[i] != '\n' && line[i] != EOS; i++)
    {
        if (isdigit(line[i]) || isspace(line[i]) || (line[i] == '*') ||
            (line[i] == ':') || (line[i] == '.'))
            ;
        else
        {
            return false;
        }
    }
    return true;
}

}
