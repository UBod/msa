/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/** 
 * Changes: 
 *
 * Mark 24-1-2007. I am now using the char "delimiter" for the delimiter when using
 * getline. This is to get around the problem of some files having '\r'.
 *
 * 10-02-07,Nigel Brown(EMBL): Changed ifstream to InFileStream to handle
 * cross-platform end-of-lines and removed delimiter member.
 *
 * 28-12-07,Paul McGettigan : replaced array processing with string processing this fixes bug #72
 *
 * 9-2-2008, Paul McGettigan : fixed problem where space after '>' but before sequence name was causing
 *                             alignment to fail due to no sequence name being read in
 * 15-2-2008, Paul McGettigan : fixed bug 91 where Pseudo -FASTA format files were not being processed as 
 *                              previously in v1.83
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "PearsonFileParser.h"

namespace clustalw
{

/**
 * Constructor for the Pearson file parser.
 * @param filePath 
 * @return 
 */
PearsonFileParser::PearsonFileParser(string filePath)
{
    fileName = filePath;
    fillCharTab();
}


/**
 * reads fasta/pearson file in one go instead of calling getSeq for
 * each single sequence.
 *
 * FIXME AW: only PearsonFileParser::getSeqRange is special, rest is the
 * same. should be defined in FileParser and then overloaded in special
 * cases like here
 */
vector<Sequence>
PearsonFileParser::getSeqRange(int firstSeq, int nSeqsToRead, string *offendingSeq)
{
    string characterSeq = "";
    string name = "";
    string title = "";
    string blank = "";
    string greater = ">";
    //_line[0] = EOS;
    vector<Sequence> seqRangeVector;
    
    string line;
    
    
    //int i, j;
    int nSeqsRead = 0;
    unsigned char c;
    char delim;
    int _currentSeqNum = 0; // Not at any sequence yet!
    
    try
    {
       delim=FileParser::getDelimiter(fileName);
       //cout << "delim = " << delim << endl;
       ifstream _fileIn;
       _fileIn.open(fileName.c_str(),ios::in);

        // Read in lines until we get to the begining of sequence firstSeq.
        string line="";

        do {
          std::getline(_fileIn,line,delim);
          if(line.substr(0,1) == greater){
              _currentSeqNum++;
          }
        } while(_currentSeqNum <firstSeq);
        
        
        while (nSeqsRead < nSeqsToRead)
        {
            // get sequence name from current line (excluded '>' and read up to first ' ' or MAXNAMES
            // remove the first char i.e. '>'
            name=line.substr(1,MAXNAMES);
            //if(name.find(">") != string::npos){
            //  andreas wilm: exit if angle bracket within header?
            //}
            
            while(name.substr(0,1)==" "){
                name=name.substr(1,MAXNAMES);
            }
            //int i;
            //i = name.find(" ");
            if(name.find(" ") != string::npos){
                name=name.substr(0,name.find(" "));
            }
            utilityObject->rTrim(&name); // also replaces linef

            name=utilityObject->blankToUnderscore(name); // replace blanks with '_'
            
            
            // Read in lines until we get to the begining of next sequence.
            
            title = ""; // No title information
            
            while(std::getline(_fileIn,line,delim) ){
                 
               string::const_iterator iterator1 = line.begin();
                while(iterator1 != line.end()){

                    // Andreas Wilm (UCD): exit if angle brackets within sequence
                    if(*iterator1=='>' && iterator1!=line.begin()) {
                        /* error output handled in Clustal.cpp
                        cerr << "\nMultiple angle brackets inside sequence found:"
                             << " invalid format.\n"
                             << "Maybe you forgot a linebreak between sequences?\n";
                        */
                        
                        parseExitCode=BADFORMAT;
                        _fileIn.close();
                        seqRangeVector.clear();
                        return seqRangeVector;
                    }
                       
                    if(*iterator1 =='\n' || *iterator1 =='\r' || *iterator1 == EOS || *iterator1 =='>'){
                        break;
                    }
                    c = *iterator1;

                    c = chartab[c];
                    if(c){
                        characterSeq.append(1,c);
                    }
                    iterator1++;
                }
                if(*iterator1 == '>'){
                    break;
               } 
            }
            
            // check sequence
            if ((int)characterSeq.length() > userParameters->getMaxAllowedSeqLength())
            {
                /* error output handled in Clustal.cpp */
                parseExitCode=SEQUENCETOOBIG;
                if (offendingSeq!=NULL)
                    offendingSeq->assign(name);
                _fileIn.close();
                seqRangeVector.clear();
                return seqRangeVector;
            }
            else if (characterSeq.length() == 0)
            {
                parseExitCode=EMPTYSEQUENCE;
                if (offendingSeq!=NULL)
                    offendingSeq->assign(name);
                _fileIn.close();
                seqRangeVector.clear();
                return seqRangeVector;
            }

            seqRangeVector.push_back(Sequence(characterSeq, name, title));
            characterSeq = "";
            nSeqsRead++;
        } // while (nSeqsRead < nSeqsToRead)

        _fileIn.close();

        return seqRangeVector;
    }

    catch(...)
    {
        cerr << "There was an exception in the PearsonFileParser::getSeqRange function.\n"
             << "Need to end program\n";
        throw 1;
    }
}



/**
 * The function getSeq is used to get the sequence 'seqNum' in the file. It returns a
 * sequence object containing the sequence.
 * Deprecated: where possible use faster getSeqRange which reads
 * sequences in one go
 * @param seqNum The number of the sequence to get.
 * @return 
 */
    Sequence PearsonFileParser::getSeq(int seqNum, string *offendingSeq)
{
    //char _line[MAXLINE + 1];
    //char tseq[MAXLINE + 1];
    //char sname[MAXNAMES + 1];
    //sname [MAXNAMES] = '\0';
    string characterSeq = "";
    string name = "";
    string title = "";
    string blank = "";
    string greater = ">";
    //_line[0] = EOS;
    
    string line;
    
    cerr << "Use of PearsonFileParser::getSeq is deprecated!\n";
    //int i, j;
    unsigned char c;
    char delim;
    int _currentSeqNum = 0; // Not at any sequence yet!
    
    try
    {
        /*
        _fileIn = new InFileStream; //nige
        _fileIn->open(fileName.c_str());  //nige
        _fileIn->seekg(0, std::ios::beg); // start at the beginning
       */
       delim=FileParser::getDelimiter(fileName);
       //cout << "delim = " << delim << endl;
       ifstream _fileIn;
       _fileIn.open(fileName.c_str(),ios::in);

        //////////////////////////////////////////////////
        //PMcG replace char array with string processing
        //////////////////////////////////////////////////

        // Read in lines until we get to the begining of sequence seqNum.
        string line="";

        do {
          std::getline(_fileIn,line,delim);
          if(line.substr(0,1) == greater){
            _currentSeqNum++;
          }
        } while(_currentSeqNum <seqNum);
        
        
        // get sequence name from current line (excluded '>' and read up to first ' ' or MAXNAMES
        // remove the first char i.e. '>'
        name=line.substr(1,MAXNAMES);
           
        //////////////////////////////////////
        // PMcG 9-2-2008 need to handle spaces at start of sequence name to conform to 1.83 handling
        //////////////////////////////////////
        while(name.substr(0,1)==" "){
          name=name.substr(1,MAXNAMES);
        }
        //int i;
        //i = name.find(" ");
        if(name.find(" ") != string::npos){
          name=name.substr(0,name.find(" "));
        } 
        name=utilityObject->blankToUnderscore(name); // replace blanks with '_'


        // Read in lines until we get to the begining of sequence seqNum.
          
        /* PMcG replace char array with string processing
        while (_currentSeqNum != seqNum)
        {
            while(*_line != '>')
            {
                if(!_fileIn->getline(_line, MAXLINE + 1))
                {
                    freeFileResources(_fileIn);
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
        for (i = 1; i <= strlen(_line); i++)
        {
            if (_line[i] != ' ')
            {
                break;
            }
        }
        strncpy(sname, _line + i, MAXNAMES); // remember entryname 
        for (i = 1; i <= strlen(sname); i++)
        {
            if (sname[i] == ' ')
            {
                break;
            }
        }
        sname[i] = EOS;
        utilityObject->rTrim(sname);
        utilityObject->blankToUnderscore(sname); // replace blanks with '_'
        name = string(sname);        
        */
        
        title = ""; // No title information

        string seqLine = "";
        while(std::getline(_fileIn,seqLine,delim) ){
          string::const_iterator iterator1 = seqLine.begin();
          while(iterator1 != seqLine.end()){
            if(*iterator1 =='\n' || *iterator1 =='\r' || *iterator1 == EOS || *iterator1 =='>'){
              break;
            }
            c = *iterator1;
            c = chartab[c];
            // PMcG 15-02-2008 bug 91
            // strip out spaces and numbers from pseudo_fasta files
            // but need to maintain gaps if present in sequence input
            // to replicate behaviour of v1.83
            //if(*iterator1 != ' ' && !isdigit(*iterator1)){
            if(c){
              characterSeq.append(1,c);
            }
            iterator1++;
          }
          if(*iterator1 == '>'){
            break;
          }
        }
 /*
        while (_fileIn->getline(_line, MAXLINE + 1))
        {
            for (i = 0; i <= MAXLINE; i++)
            {
                c = _line[i];
                if (c == '\n' || c == EOS || c == '>')
                {
                    break;
                }

                c = chartab[c];
                if (c)
                {
                    characterSeq += c;
                }
            }
            if (c == '>')
            {
                break;
            }
        }
*/

        _fileIn.close();

        if ((int)characterSeq.length() > userParameters->getMaxAllowedSeqLength())
        {
            parseExitCode=SEQUENCETOOBIG;
            // return empty seq
            return Sequence(blank, blank, blank);
        }
        else if (characterSeq.length() == 0)
        {
            parseExitCode=EMPTYSEQUENCE;
            // return empty seq
            return Sequence(blank, blank, blank);
        }
        
        return Sequence(characterSeq, name, title);
    }

    catch(...)
    {
        cerr << "There was an exception in the PearsonFileParser::getSeq function.\n"
             << "Need to end program\n";
        throw 1;
    }
}

/**
 * The function countSeqs, counts the number of sequences in a file.
 * @return The number of sequences in the file.
 */
int PearsonFileParser::countSeqs()
{
    //char line[1000 + 1];
    int _nseqs = 0;
    string line2;
    char delim;

    try
    {
        //_fileIn = new InFileStream;  //nige
        //_fileIn->open(fileName.c_str());  //nige
        delim=FileParser::getDelimiter(fileName);
        ifstream _fileIn;
        _fileIn.open(fileName.c_str(),ios::in);

    
        if(!_fileIn.is_open())
        {
            return 0; // No sequences found!
        }
    
        /* while ((*_fileIn) >> line2/@_fileIn->getline(line, 1000 + 1)@/)
           {
           /@if(_nseqs == 50)
           {
           cout << "\n\n" << line << "\n\n";
           throw 1;
           }@/
            */
        while (std::getline(_fileIn,line2,delim)) {                 
            if (line2[0] == '>')
                {
                    _nseqs++;
                }
        }
        _fileIn.close();
        return _nseqs;
    }
    catch(...)
    {
        freeFileResources(_fileIn);
        cerr << "An exception has occured in the function PearsonFileParser::countSeqs()\n"
             << "Program needs to terminate.\nPlease contact the Clustal developers\n";
        throw 1;
    }
}

/**
 * There is no secondary structure information in the Pearson file. This is here to 
 * set the structPenalties to NONE.
 * @param gapPenaltyMask 
 * @param secStructMask 
 * @param secStructName 
 * @param structPenalties 
 * @param length 
 */
void PearsonFileParser::getSecStructure(vector<char>& gapPenaltyMask, 
                         vector<char>& secStructMask, string& secStructName, 
                          int &structPenalties, int length)
{
    structPenalties = NONE;
}

}


