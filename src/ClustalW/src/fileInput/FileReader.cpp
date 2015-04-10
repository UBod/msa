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
 * 23-03-07,Nigel Brown(EMBL): added call to Alignment::testUniqueNames() to
 * test sequence names prior to appending; moved
 * Alignment::checkAllNamesDifferent() test into block handling loading of new
 * sequences.
 *
 * 02-04-07,Nigel Brown(EMBL): commented out the "No sequences in file..."
 * warnings in seqInput() and profileInput()and enabled the
 * higher-up-the-stack counterpart in MainWindow::errorHandler(), since the
 * former was causing a crash during the repaint after the user accepts the
 * message, because some state (what exactly?) is unstable at this depth after
 * a sequence load error, specifically when loading an empty file after
 * already loading some sequences.
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "FileReader.h"


namespace clustalw
{

FileReader::FileReader()
{
    fileIn = new InFileStream;
    structPenalties = NONE; // I think this should be ok.
}

FileReader::~FileReader()
{
    delete fileIn;
}


/* check if we've read sequences without any information, i.e. header
 *  but no sequence at all
 */
bool FileReader::noEmptySequence(vector<Sequence> seqVector, string *offendingSeq)
{
    vector<Sequence>::iterator si;
    for (si = seqVector.begin(); si != seqVector.end(); si++) {
        if (si->isEmpty()) {
            offendingSeq->assign(si->getName());
            return false;
        }
    }
    return true;
}



/*
 * The function seqInput is called from the interactive menu only. This is because it
 * allows us to append seqs. But this would not happen on command line.
 * It is called twice in the interactive menu, both times with append as false. It calls
 * the readseqs function to do the work.
 */
    int FileReader::seqInput(Alignment *alignPtr, bool append, string *offendingSeq)
{
    int code;
    
    if(userParameters->getMenuFlag()) 
    {
        cout << "\n\nSequences should all be in 1 file.\n"; 
        cout << "\n7 formats accepted: \n";
        cout << "NBRF/PIR, EMBL/SwissProt, Pearson (Fasta), GDE, Clustal, GCG/MSF, \
                 RSF.\n\n\n";
    }

    if (append)
    {
        int numSeqsAlready = alignPtr->getNumSeqs();
        code = readSeqs(alignPtr, numSeqsAlready + 1, offendingSeq);
    }
    else
    {
        code = readSeqs(alignPtr, 1, offendingSeq);  //  1 is the first seq to be read 
    }

    if(code == OK)
    {
        userParameters->setStructPenalties1(false);
        userParameters->setStructPenalties2(false);

        alignPtr->clearSecStruct1(); 
        alignPtr->clearSecStruct2(); 
        
        string typeOfAlign = userParameters->getDNAFlag() ? "DNA" : "PROTEIN";
        cout << "Sequences assumed to be "
             << typeOfAlign << endl;
             
        if (userParameters->getMenuFlag()) 
        {
            cout << "\n\n";
            alignPtr->printSequencesAddedInfo();  
        }    
        if(userParameters->getDNAFlag()) 
        {
            userParameters->setDNAMultiGap();
        }
        else 
        {
            userParameters->setProtMultiGap();
        }
        userParameters->setEmpty(false);
    }
    else if(code == NOSEQUENCESINFILE) // no sequences
    {
        /* 02-04-07,nige: this causes a fatal repaint: let
         *  MainWindow::errorHandler deal with it.
         */
        //utilityObject->error("No sequences in file!  Bad format?\n");
        return code;
    }
    else
    {
        return code;
    }
    return code;  
}

/*
 * The function readSeqs will create the FileParser that is needed and
 * it will read in all of the sequences. Depending on the filetype, it
 * will also check for secondary structure information.
 *
 * If there is a problem with one of the sequences its name will be
 * assigned to offendingSeq
 */
    int FileReader::readSeqs(Alignment *alignPtr, int firstSeq, string *offendingSeq)
{
    string line;
    string fileName;
    string linuxFilePath = "file://";
    //    static char *seq1, sname1[MAXNAMES + 1], title[MAXTITLES + 1];
    //int i, j;
    int noSeqs;
    //static int l1;
    static bool DNAFlag1;
    vector<Sequence> seqVector;
    vector<Sequence> seqRangeVector;
    
    auto_ptr<FileParser> fileParser; // Means we dont need to delete it!
    
    vector<int> _outputIndex; // Local version of outputIndex.
    

    if (userParameters->getMenuFlag())
    {
        utilityObject->getStr(string("Enter the name of the sequence file "), line);
    }
    else
    {
        line = userParameters->getSeqName();
    }
    if (line.length() == 0)
    {
        // We have no fileName, so return -1
        return -1;
    }
    
    // If we have file:// in the string, remove it!!
    string::size_type loc = line.find(linuxFilePath, 0);
    if(loc != string::npos)
    {
        fileName = line.substr(loc + linuxFilePath.length());
        line = fileName;
    }
    // Now we have the file name, we must open the file.
    fileIn->open(line.c_str());  //nige
    if (fileIn->fail())
        return CANNOTOPENFILE; 

    sequenceFileName = line;
    
    // Check if the file exists!
    if (!fileIn->is_open())
    {
        utilityObject->error("Cannot open input file (%s)\n", line.c_str());
        return CANNOTOPENFILE; 
    }
        
    userParameters->setSeqName(line);
   
    // NOTE I made a change here because the profile2Name was not stored!
    if(userParameters->getProfileNum() == 2 && userParameters->getProfile2Name().empty())
    {
        userParameters->setProfile2Name(line);
    }
    else if(userParameters->getProfileNum() == 1 && userParameters->getProfile1Name().empty())
    {
        userParameters->setProfile1Name(line);
    }
    
    noSeqs = 0;
    
    checkInfile(&noSeqs, fileParser);

    if (noSeqs == 0)
    {
        return NOSEQUENCESINFILE;
    }


    seqRangeVector = fileParser->getSeqRange(1, noSeqs, offendingSeq);
    if (seqRangeVector.size()==0)
        return fileParser->getParseExitCode();
    // FIXME Andreas Wilm (UCD): noEmptySequence check should be done
    // internally by FileParser instances
    if (noEmptySequence(seqRangeVector, offendingSeq) == false) {
        return EMPTYSEQUENCE; // Error there are same names.
    }
     
    for (int i=0; i < (int)seqRangeVector.size(); i++) {
        // Andreas Wilm (UCD): fixed wrong default output order
        // which is important when no alignment (convert!) is done
        // _outputIndex.push_back(i); // default output order
        _outputIndex.push_back(i+1); // default output order

        Sequence tempSeq = seqRangeVector[i];

        if (!userParameters->getExplicitDNAFlag())
        {
            DNAFlag1 = tempSeq.checkDNAFlag(); // check DNA/Prot
            if (i == 1)
            {
                userParameters->setDNAFlag(DNAFlag1);
            }
        } // type decided by first seq
        else
        {
            DNAFlag1 = userParameters->getDNAFlag();
        }
              
        seqVector.push_back(tempSeq);

    }
    
    if(firstSeq == 1) // New mulitple alignment, or else profile1
    {
        int prfNum = userParameters->getProfileNum();
        alignPtr->addSequences(&seqVector);
        //test names after saving
        if(alignPtr->checkAllNamesDifferent(offendingSeq) == false) {  //nige moved
            return ALLNAMESNOTDIFFERENT; // Error there are same names.
        }
        
        userParameters->setProfileNum(prfNum);
        if(!alignPtr->addOutputIndex(&_outputIndex)) 
        {
            return OTHERERROR;            
        }
    }
    else // profile2
    {
        //test names before appending
        if (! alignPtr->testUniqueNames(&seqVector, offendingSeq)) {  //nige added
            return ALLNAMESNOTDIFFERENT;
        }
        alignPtr->appendSequences(&seqVector);
        if(!alignPtr->appendOutputIndex(&_outputIndex))
        {
            return OTHERERROR;            
        }
    }
    
    // Then need to pass the length to the parser to get the secondary structure info.
    
    int maxAlnLength; // Get values from the Alignment object.
    maxAlnLength = alignPtr->getMaxAlnLength();    
    
    // look for a feature table / gap penalty mask (only if this is a profile) 
    if (userParameters->getProfileNum() > 0)
    {
        structPenalties = NONE;
        secStructMask.clear();
        gapPenaltyMask.clear();
        secStructName = "";
        
        fileParser->getSecStructure(gapPenaltyMask, secStructMask, secStructName,
                                  structPenalties, maxAlnLength);
        
        // Andreas Wilm (UCD): bug 114
        // after calling getSecStructure gapPenaltyMask needs to be
        // allocated (it's copied/linked to later)
        if (gapPenaltyMask.empty())
        {
            gapPenaltyMask.resize(secStructMask.size()); 
        }
    }    
  
    if (fileIn->is_open())
    {
        fileIn->close();
    }
    return OK; // return the number of seqs. read in this call
}

    int FileReader::readCharacterSeqs(Alignment* alignPtr, int firstSeq, string *offendingSeq, ClustalWInput *input) {
		string line;
	    int noSeqs;
	    //static int l1;
	    static bool DNAFlag1;
	    vector<Sequence> seqVector;
	    vector<Sequence> seqRangeVector;

	    vector<int> _outputIndex; // Local version of outputIndex.

	    auto_ptr<FileParser> fileParser; // Means we dont need to delete it!

	    if (userParameters->getMenuFlag())
	    {
	        utilityObject->getStr(string("Enter the name of the sequence file "), line);
	    }
	    else
	    {
	        line = userParameters->getSeqName();
	    }
	    if (line.length() == 0)
	    {
	        // We have no fileName, so return -1
	        return -1;
	    }

	    sequenceFileName = line;
	    userParameters->setSeqName(line);

	    // NOTE I made a change here because the profile2Name was not stored!
	    if(userParameters->getProfileNum() == 2 && userParameters->getProfile2Name().empty())
	    {
	        userParameters->setProfile2Name(line);
	    }
	    else if(userParameters->getProfileNum() == 1 && userParameters->getProfile1Name().empty())
	    {
	        userParameters->setProfile1Name(line);
	    }

	    noSeqs = (*input).inputSeqs.size();

	    if (noSeqs == 0) {
	        return NOSEQUENCESINFILE;
	    }


	    seqRangeVector = fileParser->getSeqRangeR(1, noSeqs, offendingSeq, input);
	    //if (seqRangeVector.size()==0)
	    //    return fileParser->getParseExitCode();
	    // FIXME Andreas Wilm (UCD): noEmptySequence check should be done
	    // internally by FileParser instances
	    if (noEmptySequence(seqRangeVector, offendingSeq) == false) {
	        return EMPTYSEQUENCE; // Error there are same names.
	    }

	    for (int i=0; i < (int)seqRangeVector.size(); i++) {
	        // Andreas Wilm (UCD): fixed wrong default output order
	        // which is important when no alignment (convert!) is done
	        // _outputIndex.push_back(i); // default output order
	        _outputIndex.push_back(i+1); // default output order

	        Sequence tempSeq = seqRangeVector[i];

	        if (!userParameters->getExplicitDNAFlag())
	        {
	            DNAFlag1 = tempSeq.checkDNAFlag(); // check DNA/Prot
	            if (i == 1)
	            {
	                userParameters->setDNAFlag(DNAFlag1);
	            }
	        } // type decided by first seq
	        else
	        {
	            DNAFlag1 = userParameters->getDNAFlag();
	        }

	        seqVector.push_back(tempSeq);

	    }

	    if(firstSeq == 1) // New mulitple alignment, or else profile1
	    {
	        int prfNum = userParameters->getProfileNum();
	        alignPtr->addSequences(&seqVector);
	        //test names after saving
	        if(alignPtr->checkAllNamesDifferent(offendingSeq) == false) {  //nige moved
	            return ALLNAMESNOTDIFFERENT; // Error there are same names.
	        }

	        userParameters->setProfileNum(prfNum);
	        if(!alignPtr->addOutputIndex(&_outputIndex))
	        {
	            return OTHERERROR;
	        }
	    }
	    else // profile2
	    {
	        //test names before appending
	        if (! alignPtr->testUniqueNames(&seqVector, offendingSeq)) {  //nige added
	            return ALLNAMESNOTDIFFERENT;
	        }
	        alignPtr->appendSequences(&seqVector);
	        if(!alignPtr->appendOutputIndex(&_outputIndex))
	        {
	            return OTHERERROR;
	        }
	    }

	    // Then need to pass the length to the parser to get the secondary structure info.

	    int maxAlnLength; // Get values from the Alignment object.
	    maxAlnLength = alignPtr->getMaxAlnLength();

	    // look for a feature table / gap penalty mask (only if this is a profile)
	    if (userParameters->getProfileNum() > 0)
	    {
	        structPenalties = NONE;
	        secStructMask.clear();
	        gapPenaltyMask.clear();
	        secStructName = "";

	        //FIXME TODO ... uncomment
	        //fileParser->getSecStructure(gapPenaltyMask, secStructMask, secStructName,
	        //                          structPenalties, maxAlnLength);

	        // Andreas Wilm (UCD): bug 114
	        // after calling getSecStructure gapPenaltyMask needs to be
	        // allocated (it's copied/linked to later)
	        if (gapPenaltyMask.empty())
	        {
	            gapPenaltyMask.resize(secStructMask.size());
	        }
	    }

	    return OK; // return the number of seqs. read in this call
}

/*
 * The function profileInput is used to read profiles into the Alignment. If the first
 * profile is already there it will read in the second profile. It returns the number
 * of seqs. If it returns -1, couldnt open file.
 */
int FileReader::profileInput(Alignment *alignPtr)
{
    int code;
    //int i, totalNumOfSeqs = 0;
    string offendingSeq;
    
    if(userParameters->getProfileNum() == 2 && userParameters->getProfile1Empty()) 
    {
        utilityObject->error("You must read in profile number 1 first\n");
        return MUSTREADINPROFILE1FIRST;
    }

    if(userParameters->getProfileNum() == 1)     // for the 1st profile 
    {
        code = readSeqs(alignPtr, 1, &offendingSeq);

        if (code == OK)
        {   
            // success; found some seqs.
            userParameters->setStructPenalties1(NONE);
            alignPtr->clearSecStruct1(); // Clear old info

            if (structPenalties != NONE) // feature table / mask in alignment
            {
                userParameters->setStructPenalties1(structPenalties);
                if (structPenalties == SECST) 
                {
                    
                    alignPtr->addSecStructMask1(&secStructMask);
                }
                alignPtr->addGapPenaltyMask1(&gapPenaltyMask);
                alignPtr->addSecStructName1(secStructName);

            }
            
            alignPtr->setProfile1NumSeqs(alignPtr->getNumSeqs());
            
            userParameters->setProfile1Empty(false);
            userParameters->setProfile2Empty(true);
            cout << "No. of seqs = " << alignPtr->getNumSeqs() << endl;
        }
        else
        {
            return code; // FIXME and offendingSeq?
        }
    }
    else
    {                    

        // first seq to be read = profile1_nseqs + 1 
        int profile1NumOfSeqs = alignPtr->getNumSeqs();
        code = readSeqs(alignPtr, profile1NumOfSeqs + 1, &offendingSeq); 
        if(code == OK)
        {
            userParameters->setStructPenalties2(NONE);
            alignPtr->clearSecStruct2(); // Clear old info
            
            if (structPenalties != NONE) // feature table / mask in alignment 
            {
                userParameters->setStructPenalties2(structPenalties);
                if (structPenalties == SECST) 
                {
                    alignPtr->addSecStructMask2(&secStructMask);
                }
                alignPtr->addGapPenaltyMask2(&gapPenaltyMask);
                alignPtr->addSecStructName2(secStructName);
            }
            cout << "No. of seqs in profile=" << alignPtr->getNumSeqs() - profile1NumOfSeqs 
                 << endl;
            cout << "Total no. of seqs     =" << alignPtr->getNumSeqs() << endl;
            userParameters->setProfile2Empty(false);
            userParameters->setEmpty(false); 
        }
        else
        {
            return code; // FIXME and offendingSeq?
        }
    }
    // Clear out the masks used in this class. This is important!!!!
    secStructMask.clear();
    gapPenaltyMask.clear();
    secStructName = "";
    
    string typeOfAlign = userParameters->getDNAFlag() ? "DNA" : "PROTEIN";
    cout << "Sequences assumed to be "
         << typeOfAlign << endl;
         
    if (userParameters->getMenuFlag())
    { 
        cout<< "\n\n";
    }
    
    alignPtr->printSequencesAddedInfo(); 
    
    if(userParameters->getDNAFlag()) 
    {
        userParameters->setDNAMultiGap();
    }
    else 
    {
        userParameters->setProtMultiGap();
    }
    return OK;
}


/*
 * The function checkInfile is used to find out which format the file is in, and it
 * also calls the appropriate function to count the seqences.
 */
void FileReader::checkInfile(int *nseqs, auto_ptr<FileParser>& fileParser)
{
    char _lineIn[MAXLINE + 1];
    int i;
    int lengthLine = 0;
    *nseqs = 0;

    while (fileIn->getline(_lineIn, MAXLINE + 1))
    {
        if (!utilityObject->blankLine(_lineIn))
        {
            break;
        }
    }
    lengthLine = strlen(_lineIn) - 1; // Mark Change 20-6-07
    
    for (i = lengthLine; i >= 0; i--)
    {
        if (isgraph(_lineIn[i]))
        {
            break;
        }
    }
    _lineIn[i + 1] = EOS;

    // Put first 7 characters into upper case! 
    for (i = 0; i <= 6 && i <= lengthLine; i++)
    {
        _lineIn[i] = toupper(_lineIn[i]);
    }

    // Create the parser to read the file.
    if (utilityObject->lineType(_lineIn, "ID"))
    {
         // EMBL/Swiss-Prot format ?
        fileParser.reset(new EMBLFileParser(sequenceFileName));
        if(userParameters->getDisplayInfo())
            cout << "Sequence format is EMBL\n";
    }
    else if (utilityObject->lineType(_lineIn, "CLUSTAL"))
    {
        fileParser.reset(new ClustalFileParser(sequenceFileName));
        if(userParameters->getDisplayInfo())
            cout << "Sequence format is CLUSTAL\n";
    }
    else if (utilityObject->lineType(_lineIn, "PILEUP")) // MSF
    {
        fileParser.reset(new MSFFileParser(sequenceFileName));
        if(userParameters->getDisplayInfo())
            cout << "Sequence format is MSF\n";
    }
    else if (utilityObject->lineType(_lineIn, "!!AA_MULTIPLE_ALIGNMENT")) // MSF
    {
        fileParser.reset(new MSFFileParser(sequenceFileName));
        userParameters->setDNAFlag(false);
        if(userParameters->getDisplayInfo())
            cout << "Sequence format is MSF\n";
    }
    else if (utilityObject->lineType(_lineIn, "!!NA_MULTIPLE_ALIGNMENT")) // MSF
    {
        fileParser.reset(new MSFFileParser(sequenceFileName));
        userParameters->setDNAFlag(true);
        if(userParameters->getDisplayInfo())
            cout << "Sequence format is MSF\n";
    }
    else if (strstr(_lineIn, "MSF") && _lineIn[strlen(_lineIn) - 1] == '.' &&
        _lineIn[strlen(_lineIn) - 2] == '.') // MSF
    {
        fileParser.reset(new MSFFileParser(sequenceFileName));
        if(userParameters->getDisplayInfo())
            cout << "Sequence format is MSF\n";
    }
    else if (utilityObject->lineType(_lineIn, "!!RICH_SEQUENCE")) // RSF
    {
        fileParser.reset(new RSFFileParser(sequenceFileName));
        if(userParameters->getDisplayInfo())
            cout << "Sequence format is RSF\n";
    }
    else if (utilityObject->lineType(_lineIn, "#NEXUS"))
    {
        //utilityObject->error("Cannot read NEXUS format\n");
        return;
    }
    else if (*_lineIn == '>')
    {
        if ((lengthLine>=3) && _lineIn[3] == ';') {
            //if(_lineIn[3] == ';') // distinguish PIR and Pearson
            { 
                // PIR 
                fileParser.reset(new PIRFileParser(sequenceFileName));
                if(userParameters->getDisplayInfo())
                    cout << "Sequence format is PIR\n";
            }
        }
        else
        {
            // PEARSON
            fileParser.reset(new PearsonFileParser(sequenceFileName));
            if(userParameters->getDisplayInfo())
                cout << "Sequence format is Pearson\n"; 
        }
    }
    else if ((*_lineIn == '"') || (*_lineIn == '%') || (*_lineIn == '#'))
    {
        // GDE format
        fileParser.reset(new GDEFileParser(sequenceFileName));
        if(userParameters->getDisplayInfo())
            cout << "Sequence format is GDE\n";

        if (*_lineIn == '%')
        {
            userParameters->setDNAFlag(false);
        }
        else if (*_lineIn == '#')
        {
            userParameters->setDNAFlag(true);
        }
    }
    else
    {
        return ;
    }
    
    try
    {
        // Get the number of sequences. This is passed back as a pointer!
        *nseqs = fileParser->countSeqs();
        // no output in 1.83: cout << "number of seqs is: " << *nseqs << "\n";
    }
    catch(exception ex)
    {
        *nseqs = 0;
    }
}


}

