/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Output routines ported from clustalx C code in 'interface.c'.
 *
 * Changes:
 *
 * 18-04-07,Nigel Brown(EMBL): clustalx code used (firsRes, length) convention
 * to define output columns, while the port uses (firsRes, lastRes). Fixed
 * problems in all 7 output routines where the conventions were mixed giving
 * over-long output blocks.
 *
 * 18-7-07: Mark (UCD), made changes to fastaOut, clustalOut, nbrfOut and gdeOut
 *
 * 9-2-08: Paul (UCD), changes to clustalOut to make output the same as 1.83
 *                     added const NAMESWIDTH which is calculated based on the
 *                     length of the longest sequence name and the upper and lower
 *                     limits specified in MAXNAMESTODISPLAY and MINNAMESTODISPLAY
 */

#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "AlignmentOutput.h"
#include <sstream>

namespace clustalw
{

AlignmentOutput::AlignmentOutput()
: clusSecStructOffset(9),
  clusSequenceOffset(15+1) // MARK made a change here for test case findRangeValues 10seqs
{
    // NOTE in the old clustal these arrays ended with NULL as a value.
    try
    {
        strongGroup.resize(9);
        strongGroup[0] = string("STA");
        strongGroup[1] = string("NEQK");
        strongGroup[2] = string("NHQK");
        strongGroup[3] = string("NDEQ");
        strongGroup[4] = string("QHRK");
        strongGroup[5] = string("MILV");
        strongGroup[6] = string("MILF");
        strongGroup[7] = string("HY");
        strongGroup[8] = string("FYW");
    
        weakGroup.resize(11);
        weakGroup[0] = string("CSA");
        weakGroup[1] = string("ATV");
        weakGroup[2] = string("SAG");
        weakGroup[3] = string("STNK");
        weakGroup[4] = string("STPA");
        weakGroup[5] = string("SGND");
        weakGroup[6] = string("SNDEQK");
        weakGroup[7] = string("NDEQHK");
        weakGroup[8] = string("NEQHRK");
        weakGroup[9] = string("FVLIM");
        weakGroup[10] = string("HFY");
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the contructor of AlignmentOutput.\n"
             << e.what() << "\n";
        throw 1;
    }
}

/*
 * The function createAlignmentOutput is used to call all the individual functions
 * for the different file types. It is possible to output the alignment in all file
 * types at the same time.
 */
void AlignmentOutput::createAlignmentOutput(Alignment* alignPtr, int firstSeq, int lastSeq, ClustalWOutput *output)
{
    int length;  
    int firstRes; // starting sequence range - Ramu
    int lastRes; // ending sequence range 
    bool rangeOK;
    
    if((firstSeq <= 0) || (lastSeq < firstSeq))
    {
        utilityObject->error("Cannot produce alignment output."
                " Incorrect call to createAlignmentOutput."
                " firstSeq = %d"
                " lastSeq = %d\n", firstSeq, lastSeq);
        return;
    }
    
    length = 0;
    firstRes = 1;
    lastRes = 0;
    rangeOK = false;
    
    try
    {
        length = alignPtr->getLengthLongestSequence();
        lastRes = length;

        if (userParameters->getRangeFromToSet()) 
        {
            firstRes = userParameters->getRangeFrom();
            lastRes = userParameters->getRangeTo();
            // Check if the numbers are ok.
            if ((firstRes > lastRes) || (firstRes == -1) || (lastRes == -1)) 
            {
                cerr << "seqrange numbers are not set properly, using default....\n";
                firstRes = 1;
                lastRes = length;
            }
            else
            {
                rangeOK = true;
            }
        }
        if (rangeOK && (lastRes > length)) 
        {
            lastRes = length;
            cout << "Seqrange " << lastRes << " is more than the " << length 
                 << "  setting it to " << length << " \n";
        }

        if (userParameters->getMenuFlag())
        {    
            cout << "Consensus length = " << lastRes << " \n";
        }

        outputRegion partToOutput;
        partToOutput._firstSeq = firstSeq;
        partToOutput._lastSeq = lastSeq;
        partToOutput._firstRes = firstRes;
        partToOutput._lastRes = lastRes;

        if(userParameters->getOutputClustal()) 
        {
            //if(clustalOutFile.get() != 0 && clustalOutFile->is_open())
            //{
                clustalOut(alignPtr, partToOutput, output);
                //if(clustalOutFile->is_open())
                //{
                    clustalOutFile->close();
                //}
                //utilityObject->info("CLUSTAL-Alignment result created"); //  [%s]\n",
                                   // clustalOutName.c_str());
            //}
        }
        if(userParameters->getOutputNbrf())  
        {
            if(nbrfOutFile.get() != 0 && nbrfOutFile->is_open())
            {            
                nbrfOut(alignPtr, partToOutput);
                if(nbrfOutFile->is_open())
                {
                    nbrfOutFile->close();
                }
                utilityObject->info("NBRF/PIR-Alignment file created  [%s]\n",
                                    nbrfOutName.c_str());                
            }
        }
        if(userParameters->getOutputGCG())  
        {
            if(gcgOutFile.get() != 0 && gcgOutFile->is_open())
            {
                gcgOut(alignPtr, partToOutput);
                if(gcgOutFile->is_open())
                {
                    gcgOutFile->close();
                }
                utilityObject->info("GCG-Alignment file created      [%s]\n",
                                    gcgOutName.c_str());                
            }
        }
        if(userParameters->getOutputPhylip())  
        {
            if(phylipOutFile.get() != 0 && phylipOutFile->is_open())
            {
                phylipOut(alignPtr, partToOutput);
                if(phylipOutFile->is_open())
                {
                    phylipOutFile->close();
                }
                utilityObject->info("PHYLIP-Alignment file created   [%s]\n",
                                    phylipOutName.c_str());                
            }
        }
        
        if(userParameters->getOutputGde())  
        {
            if(gdeOutFile.get() != 0 && gdeOutFile->is_open())
            {
                gdeOut(alignPtr, partToOutput);
                if(gdeOutFile->is_open())
                {
                    gdeOutFile->close();
                }
                utilityObject->info("GDE-Alignment file created      [%s]\n",
                                    gdeOutName.c_str());                
            }
        }
        
        if(userParameters->getOutputNexus())  
        {
            if(nexusOutFile.get() != 0 && nexusOutFile->is_open())
            {
                nexusOut(alignPtr, partToOutput);
                if(nexusOutFile->is_open())
                {
                    nexusOutFile->close();
                }
                utilityObject->info("NEXUS-Alignment file created    [%s]\n",
                                    nexusOutName.c_str());                
            }
        }

        if(userParameters->getOutputFasta())  
        {
            //if(fastaOutFile.get() != 0 && fastaOutFile->is_open())
            //{
                fastaOut(alignPtr, partToOutput, output);
                //if(fastaOutFile->is_open())
                //{
                    fastaOutFile->close();
                //}
                utilityObject->info("Fasta-Alignment result created    [%s]\n",
                                    fastaOutName.c_str());                
            //}
        }
        
        // Output the alignment to the screen if it is required.
        if (userParameters->getShowAlign() && userParameters->getMenuFlag()) 
        {
            showAlign();
        }
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the createAlignmentOutput function.\n"
             << e.what() << "\n";
        throw 1;
    }
    return;
}

bool AlignmentOutput::QTOpenFilesForOutput(AlignmentFileNames fileNames)
{
    if(!userParameters->getOutputClustal() && 
       !userParameters->getOutputNbrf() && !userParameters->getOutputGCG() &&
       !userParameters->getOutputPhylip() && !userParameters->getOutputGde() &&
       !userParameters->getOutputNexus() && !userParameters->getOutputFasta()) 
    {
         utilityObject->error("You must select an alignment output format\n");
         return false;
    }
        
    if(fileNames.clustalFile == "" && fileNames.fastaFile == "" && fileNames.gcgFile == ""
       && fileNames.gdeFile == "" && fileNames.nexusFile == "" && fileNames.nrbfFile == ""
       && fileNames.phylipFile == "")
    {
        utilityObject->error("No names for output files. Cannot output alignment.\n");
        return false;
    }
    
    if(fileNames.clustalFile != "")
    {
        clustalOutName = fileNames.clustalFile;
        if(!openExplicitFile(clustalOutFile, clustalOutName))
        { 
            return false;
        }
    }  
    if(fileNames.fastaFile != "")
    {
        fastaOutName = fileNames.fastaFile;
        if(!openExplicitFile(fastaOutFile, fastaOutName))
        { 
            return false;
        }
    }
    if(fileNames.gcgFile != "")
    {
        gcgOutName = fileNames.gcgFile;
        if(!openExplicitFile(gcgOutFile, gcgOutName))
        { 
            return false;
        }
    }
    if(fileNames.gdeFile != "")
    {
        gdeOutName = fileNames.gdeFile;
        if(!openExplicitFile(gdeOutFile, gdeOutName))
        { 
            return false;
        }
    }
    if(fileNames.nexusFile != "")
    {
        nexusOutName = fileNames.nexusFile;
        if(!openExplicitFile(nexusOutFile, nexusOutName))
        { 
            return false;
        }
    }
    if(fileNames.nrbfFile != "")
    {
        nbrfOutName = fileNames.nrbfFile;
        if(!openExplicitFile(nbrfOutFile, nbrfOutName))
        { 
            return false;
        }
    }
    if(fileNames.phylipFile != "")
    {
        phylipOutName = fileNames.phylipFile;
        if(!openExplicitFile(phylipOutFile, phylipOutName))
        { 
            return false;
        }
    } 
    return true;                         
}
/*
 * The function openAlignmentOutput opens a file for output. It returns true if it
 * has been successful and false if it has not been successful
 */
bool AlignmentOutput::openAlignmentOutput(string path)
{
    if(!userParameters->getOutputClustal() && 
       !userParameters->getOutputNbrf() && !userParameters->getOutputGCG() &&
       !userParameters->getOutputPhylip() && !userParameters->getOutputGde() &&
       !userParameters->getOutputNexus() && !userParameters->getOutputFasta()) 
    {
         utilityObject->error("You must select an alignment output format\n");
         return false;
    }
    string _fileNameToOutput = path;
    if(_fileNameToOutput == "")
    {
      /*_fileNameToOutput = userParameters->getSeqName();*/
      /* BUG 166 , file extension, FS, 2009-01-26 */
      utilityObject->getPath(userParameters->getSeqName(), &_fileNameToOutput);
    }
    
    if(userParameters->getOutputClustal()) 
    {
        if (userParameters->getOutfileName() != "") 
        {
            clustalOutName = userParameters->getOutfileName();
            if(!openExplicitFile(clustalOutFile, clustalOutName))
            { 
                return false;
            }
        }
        else 
        {
            clustalOutName = openOutputFile(clustalOutFile, 
                             "\nEnter a name for the CLUSTAL output file ",
                             _fileNameToOutput, "aln");
                             
            if(clustalOutName == "")
            {
                return false; // We have not been successful.
            }
        }
    }

    if(userParameters->getOutputNbrf()) 
    {
        if (userParameters->getOutfileName() != "") 
        {
            nbrfOutName = userParameters->getOutfileName();
            if(!openExplicitFile(nbrfOutFile, nbrfOutName))
            { 
                return false;
            }            
        }
        else
        {
            nbrfOutName = openOutputFile(nbrfOutFile,
                      "\nEnter a name for the NBRF/PIR output file", _fileNameToOutput,
                       "pir");
                      
            if(nbrfOutName == "")
            {
                return false; // We have not been successful.
            }
        }
    }

    if(userParameters->getOutputGCG()) 
    {
        if (userParameters->getOutfileName() != "") 
        {
            gcgOutName = userParameters->getOutfileName();
            if(!openExplicitFile(gcgOutFile, gcgOutName))
            { 
                return false;
            }   
        }
        else
        {
            gcgOutName = openOutputFile(gcgOutFile,
                      "\nEnter a name for the GCG output file", _fileNameToOutput, "msf");
                      
            if(gcgOutName == "")
            {
                return false; // We have not been successful.
            }
        }
    }

    if(userParameters->getOutputPhylip()) 
    {
        if (userParameters->getOutfileName() != "") 
        {
            phylipOutName = userParameters->getOutfileName();
            if(!openExplicitFile(phylipOutFile, phylipOutName))
            { 
                return false;
            }
        }
        else
        {
            phylipOutName = openOutputFile(phylipOutFile,
                      "\nEnter a name for the PHYLIP output file", _fileNameToOutput, "phy");
                      
            if(phylipOutName == "")
            {
                return false; // We have not been successful.
            }
        }
    }

    if(userParameters->getOutputGde()) 
    {
        if (userParameters->getOutfileName() != "") 
        {
            gdeOutName = userParameters->getOutfileName();
            if(!openExplicitFile(gdeOutFile, gdeOutName))
            { 
                return false;
            }
        }
        else
        {
            gdeOutName = openOutputFile(gdeOutFile, 
                             "\nEnter a name for the GDE output file     ",
                             _fileNameToOutput, "gde");
                 
            if(gdeOutName == "")
            {
                return false; // We have not been successful.
            }                       
        }
    }

    if(userParameters->getOutputNexus()) 
    {
        if (userParameters->getOutfileName() != "") 
        {
            nexusOutName = userParameters->getOutfileName();
            if(!openExplicitFile(nexusOutFile, nexusOutName))
            { 
                return false;
            }
        }
        else
        {
            nexusOutName = openOutputFile(nexusOutFile, 
                             "\nEnter a name for the NEXUS output file   ",
                             _fileNameToOutput, "nxs");
                 
            if(nexusOutName == "")
            {
                return false; // We have not been successful.
            }                       
        }
    }
   
    if(userParameters->getOutputFasta()) 
    {
        if (userParameters->getOutfileName() != "") 
        {
            fastaOutName = userParameters->getOutfileName();
            if(!openExplicitFile(fastaOutFile, fastaOutName))
            { 
                return false;
            }
        }
        else
        {
            fastaOutName = openOutputFile(fastaOutFile, 
                             "\nEnter a name for the Fasta output file ",
                             _fileNameToOutput, "fasta");
                 
            if(fastaOutName == "")
            {
                return false; // We have not been successful.
            }
        }
    }
  
    return true;

}


/*
 * The function openExplicitFile is used to open a file when we have the name already.
 * It returns true if the file has been opened correctly, false otherwise.
 */
bool AlignmentOutput::openExplicitFile(auto_ptr<ofstream>& outFile, string fileName)
{
    if (fileName == "") 
    {
        cerr << "Bad output file [" << fileName << "]\n";
        utilityObject->error("Bad output file [%s]\n", fileName.c_str());
        return false;
    }
    
    outFile.reset(new ofstream(fileName.c_str(), ofstream::trunc));
    if(!outFile->is_open()) 
    {
        utilityObject->error("Cannot open output file [%s]\n", fileName.c_str());
        return false;
    }
    return true;
}

/*
 * The function openOutputFile is used to open a file when we dont have the name yet.
 * This function returns the name if it has been successful. If it has not been 
 * successful it returns a blank string ""
 */
string AlignmentOutput::openOutputFile(auto_ptr<ofstream>& outFile, string prompt, 
                                       string path, string fileExtension)
{
    string temp;
    string _fileName; // Will return this name.
    string message;
    _fileName = path + fileExtension;

    if(_fileName.compare(userParameters->getSeqName()) == 0) 
    {
        cout << "Output file name is the same as input file.\n";
        if (userParameters->getMenuFlag()) 
        {
            message = "\n\nEnter new name to avoid overwriting  [" + _fileName + "]";

            utilityObject->getStr(message, temp);
            if(temp != "")
            {
                _fileName = temp;
            }
        }
    }
    else if (userParameters->getMenuFlag()) 
    {

        message = prompt + " [" + _fileName + "]";
        utilityObject->getStr(message, temp);
        if(temp != "")
        {
            _fileName = temp;
        }
    }

    outFile.reset(new ofstream(_fileName.c_str(), ofstream::trunc));
    if(!outFile->is_open()) 
    {
        utilityObject->error("Cannot open output file [%s]\n", _fileName.c_str());
        return "";
    }    
    return _fileName;
}

/*
 * fastaOut: output the alignment in FASTA format
 */
void AlignmentOutput::fastaOut(Alignment* alignPtr, outputRegion partToOutput, ClustalWOutput *output)
{
	stringstream line;
    char residue;
    int val;
    int i, ii;
    int j, slen;    
    int lineLength;
    int firstRes = partToOutput._firstRes;
    int lastRes = partToOutput._lastRes;
    int firstSeq = partToOutput._firstSeq;
    int lastSeq = partToOutput._lastSeq;
    rangeNum rnum;  
    cout << "firstres = " << firstRes << " lastres = " << lastRes << "\n";
    try
    {
        const SeqArray* alignment = alignPtr->getSeqArray(); //NOTE june29
        vector<char> sequence;
        sequence.assign(lastRes + 1, '0');     
    
        lineLength = PAGEWIDTH - alignPtr->getMaxNames();
                
        lineLength = lineLength - lineLength % 10; // round to a multiple of 10
        if (lineLength > LINELENGTH || lineLength <= 0) // Mark 18-7-07
        {    
            lineLength = LINELENGTH;
        }


        for(ii = firstSeq; ii <= lastSeq; ii++) 
        {
            i = alignPtr->getOutputIndex(ii - 1);
            slen = 0;
            //for(j = firstRes; j < firstRes + lastRes; j++) //- nige
            for(j = firstRes; j <= lastRes; j++) //- nige
            {
                if (j <= alignPtr->getSeqLength(i))
                {
                    //val = alignPtr.getResidue(i, j);
                    val = (*alignment)[i][j]; // NOTE june29
                }
                else
                { 
                    val = -3;
                }
                if((val == -3) || (val == 253))
                {
                    break;
                }
                else if((val < 0) || (val > userParameters->getMaxAA())) 
                {
                    residue = '-';
                }
                else 
                {
                    residue = userParameters->getAminoAcidCode(val);
                }
                if (userParameters->getLowercase())
                {    
                    sequence[j - firstRes] = residue;
                }
                else
                {
                    sequence[j - firstRes] = residue;
                }
                slen++;
            }

            //(*fastaOutFile) << ">" << nameonly(alignPtr->getName(i));
            line << ">" << nameonly(alignPtr->getName(i));
            if(userParameters->getSeqRange()) 
            {
                findRangeValues(alignPtr, &rnum, firstRes, lastRes, ii);
                //(*fastaOutFile) << "/" << rnum.start << "-" << rnum.end;
                line << "/" << rnum.start << "-" << rnum.end;
                
            }
            //(*fastaOutFile) << "\n";
            output->msa.push_back(line.str());
            for(j = 1; j <= slen; j++) 
            {
                //(*fastaOutFile) << sequence[j-1];
            	line.str("");
                line << sequence[j-1];
            	if((j % lineLength == 0) || (j == slen))
                { 
                    //(*fastaOutFile) << "\n";
            		output->msa.push_back(line.str());
                }
            }
        }
        //fastaOutFile->close();
        cout << "FASTA string created!\n";
    }
    catch(const bad_alloc& e)
    {
        //fastaOutFile->close();
        cerr << "A bad_alloc exception has occured in the fastaOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(VectorOutOfRange e)
    {
        //fastaOutFile->close();
        cerr << "An exception has occured in the fastaOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(...)
    {
        //fastaOutFile->close();
        cerr << "An exception has occured in the fastaOut function.\n";
        throw 1;
    }
    return;
}

/* 
 * gcgOut: output the alignment in gcg file format
 */
void AlignmentOutput::gcgOut(Alignment* alignPtr, outputRegion partToOutput)
{
    char residue;
    int val;
    int i, ii, chunks, block;
    int j, k, pos1, pos2;    
    long grandChecksum;
 
    try
    {  
        int firstRes = partToOutput._firstRes;
        int lastRes = partToOutput._lastRes;
        int firstSeq = partToOutput._firstSeq;
        int lastSeq = partToOutput._lastSeq;    
        rangeNum rnum;
        vector<char> sequence;
        vector<int> allChecks;
        int _maxLength = alignPtr->getMaxAlnLength();
        sequence.assign(_maxLength + 1, '0');
        allChecks.assign(lastSeq + 1, 0);
        const SeqArray* alignment = alignPtr->getSeqArray();
        
        for(i = firstSeq; i <= lastSeq; i++) 
        {
            //for(j = firstRes; j <= firstRes + lastRes - 1; j++) //- nige
            for(j = firstRes; j <= lastRes; j++) //- nige
            {
                if (j <= alignPtr->getSeqLength(i))
                {
                    //val = alignPtr.getResidue(i, j);
                    val = (*alignment)[i][j]; // NOTE june29
                }
                else
                { 
                    val = -3;
                }
                if((val == -3) || (val == 253)) 
                {
                    break;
                }
                else if((val < 0) || (val > userParameters->getMaxAA()))
                {
                    residue = '.';
                }
                else 
                {
                    residue = userParameters->getAminoAcidCode(val);
                }
                sequence[j - firstRes + 1] = residue;
            }
            // pad any short sequences with gaps, to make all sequences the same length 
            for(; j <= firstRes + lastRes - 1; j++)
            { 
                sequence[j - firstRes + 1] = '.';
            }
            allChecks[i] = SeqGCGCheckSum(&sequence, lastRes);
        }    
        int _index;
        grandChecksum = 0;
        for(i = 1; i <= alignPtr->getNumSeqs(); i++)
        {
            _index = alignPtr->getOutputIndex(i - 1);
            grandChecksum += allChecks[_index];
        }
        grandChecksum = grandChecksum % 10000;
        
        (*gcgOutFile) << "PileUp\n\n";
        (*gcgOutFile) << "\n\n   MSF:" << setw(5) << lastRes << "  Type: ";
                        
        if(userParameters->getDNAFlag())
        {
            (*gcgOutFile) << "N";
        }
        else
        {
            (*gcgOutFile) << "P";
        }

        (*gcgOutFile) << "    Check:" << setw(6) << grandChecksum << "   .. \n\n";
        float _seqWeight;

        int length = lastRes - firstRes + 1; //- nige

        for(ii = firstSeq; ii <= lastSeq; ii++)  
        {
            i = alignPtr->getOutputIndex(ii - 1);
            _seqWeight = (alignPtr->getSeqWeight(i - 1) * 100) / (float) INT_SCALE_FACTOR;

            (*gcgOutFile) << " Name: " << alignPtr->getName(i) << " oo  Len:"
                        //<< setw(5) << lastRes << "  Check:" << setw(6) << allChecks[i]  //- nige
                          << setw(5) << length << "  Check:" << setw(6) << allChecks[i] 
                          << "  Weight:  " << fixed << setprecision(1) << _seqWeight << "\n";
        }  
        (*gcgOutFile) << "\n//\n";
        
        chunks = length / GCG_LINELENGTH; //- nige
        if(length % GCG_LINELENGTH != 0)  //- nige
        {
            ++chunks;
        }
  
        for(block = 1; block <= chunks; block++) 
        {
            (*gcgOutFile) << "\n\n";
            pos1 = ((block - 1) * GCG_LINELENGTH) + 1;
            pos2 = (length < pos1 + GCG_LINELENGTH - 1)? length : pos1 + GCG_LINELENGTH - 1;//- nige
            
            for(ii = firstSeq; ii <= lastSeq; ii++) 
            {
                i = alignPtr->getOutputIndex(ii - 1);
                if (!userParameters->getSeqRange()) 
                {
                    (*gcgOutFile) << "\n" << setw(alignPtr->getMaxNames() + 5) << left
                                  << alignPtr->getName(i) << " ";
                }
                else 
                {
                    findRangeValues(alignPtr, &rnum, firstRes, lastRes, ii);
                    std::stringstream ss;
                    std::stringstream ss2;
                    string rangeStart;
                    string rangeEnd;
                    ss << rnum.start;
                    ss >> rangeStart;
                    ss2 << rnum.end;
                    ss2 >> rangeEnd;
                    string nameAndRange = nameonly(alignPtr->getName(i)) + "/" + rangeStart;
                    nameAndRange += "-" + rangeEnd;
                    (*gcgOutFile) << "\n" << setw(alignPtr->getMaxNames() + 15) << left
                                  << nameAndRange;
                }
                for(j = pos1, k = 1; j <= pos2; j++, k++) 
                {
            
                    //JULIE -
                    //check for int sequences - pad out with '.' characters to end 
            
                    if (j + firstRes - 1 <= alignPtr->getSeqLength(i))
                    {
                        //val = alignPtr.getResidue(i, j + firstRes - 1);
                        val = (*alignment)[i][j + firstRes - 1]; // NOTE june29
                    }
                    else
                    {
                        val = -3;
                    }
                    if((val == -3) || (val == 253))
                    {
                        residue = '.';
                    }
                    else if((val < 0) || (val > userParameters->getMaxAA()))
                    {
                        residue = '.';
                    }
                    else 
                    {
                        residue = userParameters->getAminoAcidCode(val);
                    }
                    (*gcgOutFile) << residue;
                    if(j % 10 == 0)
                    { 
                        (*gcgOutFile) << " ";
                    }
                }
            }
        }
        (*gcgOutFile) << "\n\n";
        gcgOutFile->close();   
    }
    catch(const bad_alloc& e)
    {
        gcgOutFile->close();
        cerr << "A bad_alloc exception has occured in the gcgOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(VectorOutOfRange e)
    {
        gcgOutFile->close();
        cerr << "An exception has occured in the gcgOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(...)
    {
        gcgOutFile->close();
        cerr << "An exception has occured in the gcgOut function.\n";
        throw 1;
    }    
}

void AlignmentOutput::phylipOut(Alignment* alignPtr, outputRegion partToOutput)
{
    int firstRes = partToOutput._firstRes;
    int lastRes = partToOutput._lastRes;
    int firstSeq = partToOutput._firstSeq;
    int lastSeq = partToOutput._lastSeq;    
    try
    {
        char residue;
        int val;
        int i, ii, chunks, block;    
        int j, k, pos1, pos2;    
        int nameLen;
        bool warn;
        vector<string> _seqNames;
        const SeqArray* alignment = alignPtr->getSeqArray();
        rangeNum rnum;

        // Push on a blank string. This is because the index starts at 1 for the seqs etc..
        _seqNames.push_back("");

        nameLen = 0;
    
        for(i = firstSeq; i <= lastSeq; i++)  
        {
            _seqNames.push_back(alignPtr->getName(i));
            ii = _seqNames.size();
            if(nameLen < ii) 
            {
                nameLen = ii;
            }
        }
        if(nameLen > 10) 
        {
            warn = false;
            for(i = 0; i < (int)_seqNames.size() - 1; i++)  
            {
                for(j = i + 1; j < (int)_seqNames.size(); j++) 
                {
                    if (_seqNames[i].compare(_seqNames[j]) == 0) 
                    {
                        warn = true;
                    }
                }
            }
            if(warn)
            {
                utilityObject->warning("Truncating sequence names to 10 characters for PHYLIP output.\nNames in the PHYLIP format file are NOT unambiguous.");
            }
            else
            {
                utilityObject->warning("Truncating sequence names to 10 characters for PHYLIP output.");
            }
        }
  
        int length = lastRes - firstRes + 1; //- nige

        chunks = length / GCG_LINELENGTH; //- nige
        if(length % GCG_LINELENGTH != 0)  //- nige
        {
            ++chunks;
        }
  
        (*phylipOutFile) << setw(6) << alignPtr->getNumSeqs() << " "
                         << setw(6) << length; //-nige
        
        for(block = 1; block <= chunks; block++) 
        {
            pos1 = ((block - 1) * GCG_LINELENGTH) + 1;
            //pos2 = (lastRes < pos1 + GCG_LINELENGTH - 1)? lastRes : pos1 + GCG_LINELENGTH - 1; //- nige
            pos2 = (length < pos1 + GCG_LINELENGTH - 1)? length : pos1 + GCG_LINELENGTH - 1;
            for(ii = firstSeq; ii <= lastSeq; ii++)  
            {
                i = alignPtr->getOutputIndex(ii - 1);
                if(block == 1)  
                {
                    if(!userParameters->getSeqRange()) 
                    {
                        string name = alignPtr->getName(i);
                        (*phylipOutFile) << "\n" << setw(10) << left 
                                             << name.substr(0,10) << " ";
                    }
                    else
                    {
                        findRangeValues(alignPtr, &rnum, firstRes, lastRes, ii);
                        std::stringstream ss;
                        std::stringstream ss2;
                        string rangeStart;
                        string rangeEnd;
                        ss << rnum.start;
                        ss >> rangeStart;
                        ss2 << rnum.end;
                        ss2 >> rangeEnd;
                        string nameAndRange = nameonly(alignPtr->getName(i)) + "/";
                        nameAndRange += rangeStart + "-" + rangeEnd;
                        (*phylipOutFile) << "\n" << setw(alignPtr->getMaxNames() + 15)
                                         << left
                                         << nameAndRange;
                    }
                }
                else
                {
                    (*phylipOutFile) << "\n           ";
                }
                for(j = pos1, k = 1; j <= pos2; j++, k++) 
                {
                    if (j + firstRes - 1 <= alignPtr->getSeqLength(i))
                    {
                        //val = alignPtr.getResidue(i, j + firstRes - 1);
                        val = (*alignment)[i][j + firstRes - 1]; // NOTE june29
                    }
                    else
                    {
                        val = -3;
                    }
                    if((val == -3) || (val == 253))
                    {
                        break;
                    }
                    else if((val < 0) || (val > userParameters->getMaxAA()))
                    {
                        residue = '-';
                    }
                    else 
                    {
                        residue = userParameters->getAminoAcidCode(val);
                    }
                    (*phylipOutFile) << residue;
                    if(j % 10 == 0) 
                    {
                        (*phylipOutFile) << " ";
                    }
                }
            }
            (*phylipOutFile) << "\n";
        }  
    }
    catch(const bad_alloc& e)
    {
        phylipOutFile->close();
        cerr << "A bad_alloc exception has occured in the phylipOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(VectorOutOfRange e)
    {
        phylipOutFile->close();
        cerr << "An exception has occured in the phylipOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(...)
    {
        phylipOutFile->close();
        cerr << "An exception has occured in the phylipOut function.\n";
        throw 1;
    }      
}

void AlignmentOutput::nexusOut(Alignment* alignPtr, outputRegion partToOutput)
{
    int firstRes = partToOutput._firstRes;
    int lastRes = partToOutput._lastRes;
    int firstSeq = partToOutput._firstSeq;
    int lastSeq = partToOutput._lastSeq;    
    try
    {
        char residue;
        int val;
        int i, ii, chunks, block;    
        int j, k, pos1, pos2;    
        const SeqArray* alignment = alignPtr->getSeqArray();
        rangeNum rnum;

        int length = lastRes - firstRes + 1; //- nige

        chunks = length / GCG_LINELENGTH; //- nige
        if(length % GCG_LINELENGTH != 0) //- nige
        {
            ++chunks;
        }
  
        (*nexusOutFile) << "#NEXUS\n";
        (*nexusOutFile) << "BEGIN DATA;\n";
        (*nexusOutFile) << "dimensions ntax=" << alignPtr->getNumSeqs() 
                        << " nchar=" << length << ";\n"; //- nige
        (*nexusOutFile) << "format missing=?\n";

	// removed - bugzilla bug 204
	/*
        (*nexusOutFile) << "symbols=\"";
        
        for(i = 0; i <= userParameters->getMaxAA(); i++)
        {
            (*nexusOutFile) << userParameters->getAminoAcidCode(i);
        }
        (*nexusOutFile) << "\"\n";
	*/

        (*nexusOutFile) << "interleave datatype=";
        bool _dnaFlag = userParameters->getDNAFlag();
        string _type = _dnaFlag ? "DNA " : "PROTEIN ";
        (*nexusOutFile) << _type;
        (*nexusOutFile) << "gap= -;\n";
        (*nexusOutFile) << "\nmatrix";
  
        for(block = 1; block <= chunks; block++) 
        {
            pos1 = ((block - 1) * GCG_LINELENGTH) + 1;
            pos2 = (length < pos1 + GCG_LINELENGTH - 1)? length : pos1 + GCG_LINELENGTH - 1; //- nige
            for(ii = firstSeq; ii <= lastSeq; ii++)  
            {
                i = alignPtr->getOutputIndex(ii - 1);
                if (!userParameters->getSeqRange()) 
                {
                    (*nexusOutFile) << "\n" << setw(alignPtr->getMaxNames() + 1) << left
                                    << alignPtr->getName(i) << " ";
                }
                else 
                {
                    findRangeValues(alignPtr, &rnum, firstRes, lastRes, ii);
                    
                    std::stringstream ss;
                    std::stringstream ss2;
                    string rangeStart;
                    string rangeEnd;
                    ss << rnum.start;
                    ss >> rangeStart;
                    ss2 << rnum.end;
                    ss2 >> rangeEnd;
                    string nameAndRange = nameonly(alignPtr->getName(i)) + "/";
                    nameAndRange += rangeStart + "-" + rangeEnd;
                    (*nexusOutFile) << "\n" << setw(alignPtr->getMaxNames() + 15)
                                    << left
                                    << nameAndRange;                    
                }
                for(j = pos1, k = 1; j <= pos2; j++, k++) 
                {
                    if (j + firstRes - 1 <= alignPtr->getSeqLength(i))
                    {
                        //val = alignPtr.getResidue(i, j + firstRes - 1);
                        val = (*alignment)[i][j + firstRes - 1]; // NOTE june29
                    }
                    else
                    {
                        val = -3;
                    }
                    if((val == -3) || (val == 253))
                    {
                        break;
                    }
                    else if((val < 0) || (val > userParameters->getMaxAA()))
                    {
                        residue = '-';
                    }
                    else 
                    {
                        residue = userParameters->getAminoAcidCode(val);
                    }
                    (*nexusOutFile) << residue;
                }
            }
            (*nexusOutFile) << "\n";
        }
        (*nexusOutFile) << ";\nend;\n";   
    }
    catch(const bad_alloc& e)
    {
        nexusOutFile->close();
        cerr << "A bad_alloc exception has occured in the nexusOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(VectorOutOfRange e)
    {
        nexusOutFile->close();
        cerr << "An exception has occured in the nexusOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(...)
    {
        nexusOutFile->close();
        cerr << "An exception has occured in the nexusOut function.\n";
        throw 1;
    }      
}

/*
 * gdeOut: print out the alignment in gde file format.
 */
void AlignmentOutput::gdeOut(Alignment* alignPtr, outputRegion partToOutput)
{
    int firstRes = partToOutput._firstRes;
    int lastRes = partToOutput._lastRes;
    int firstSeq = partToOutput._firstSeq;
    int lastSeq = partToOutput._lastSeq;    

    int length = lastRes - firstRes + 1; //- nige

    try
    {
        char residue;
        int val;
        vector<char> _ssMask1, _ssMask2, seq;
        int i, ii;
        int j, slen;    
        int lineLength;
        const SeqArray* alignment = alignPtr->getSeqArray();
        rangeNum rnum;

        seq.assign(alignPtr->getMaxAlnLength() + 1, '0');
        
        // decide the line length for this alignment - maximum is LINELENGTH 
        lineLength = PAGEWIDTH - alignPtr->getMaxNames();
        lineLength = lineLength - lineLength % 10; // round to a multiple of 10
        if (lineLength > LINELENGTH || lineLength <= 0) 
        {
            lineLength = LINELENGTH;
        }
  
        // If we are using the secondary structure info from profile 1, set it up
        if (userParameters->getStructPenalties1() == SECST && 
            userParameters->getUseSS1() == true) 
        {
            int _lengthSeq1 = alignPtr->getSeqLength(1);
            vector<char>* _secStructMask1 = alignPtr->getSecStructMask1();
            _ssMask1.assign(_lengthSeq1 + 10, 0);
            
            for (i = 0;i < _lengthSeq1;i++)
            {
                _ssMask1[i] = _secStructMask1->at(i);
            }
            printSecStructMask(_lengthSeq1, _secStructMask1, &_ssMask1);
        }
        
        // If we are using the secondary structure info from profile 2, set it up
        if (userParameters->getStructPenalties2() == SECST && 
            userParameters->getUseSS2() == true) 
        {
            int indexProfile2FirstSeq = alignPtr->getProfile1NumSeqs() + 1;
            int lengthSeqProfile2 = alignPtr->getSeqLength(indexProfile2FirstSeq);
            _ssMask2.assign(lengthSeqProfile2 + 10, 0);
            vector<char>* _secStructMask2 = alignPtr->getSecStructMask2();
            
            for (i=0;i < lengthSeqProfile2;i++)
            {
                _ssMask2[i] = _secStructMask2->at(i);
            }
            printSecStructMask(lengthSeqProfile2, _secStructMask2, &_ssMask2);  
        }

        bool _dnaFlag = userParameters->getDNAFlag();
        string _prefix = _dnaFlag ? "#" : "%";
        
        for(ii = firstSeq; ii <= lastSeq; ii++) 
        {
            i = alignPtr->getOutputIndex(ii - 1);
            
            (*gdeOutFile) << _prefix;
            if(!userParameters->getSeqRange()) 
            {
                (*gdeOutFile) << alignPtr->getName(i) << "\n";
            }
            else 
            {
                findRangeValues(alignPtr, &rnum, firstRes, lastRes, ii);
                (*gdeOutFile) << nameonly(alignPtr->getName(i)) << "/"
                              << rnum.start << "-" << rnum.end << "\n";
            }
            slen = 0;

            //for(j = firstRes; j < firstRes + lastRes; j++) //- nige
            for(j = firstRes; j <= lastRes; j++) //- nige
            {
                if (j <= alignPtr->getSeqLength(i))
                {
                    //val = alignPtr.getResidue(i, j);
                    val = (*alignment)[i][j]; // NOTE june29
                }
                else
                { 
                    val = -3;
                }
                if((val == -3) || (val == 253))
                {
                    break;
                }
                else if((val < 0) || (val > userParameters->getMaxAA()))
                {
                    residue = '-';
                }
                else 
                {
                    residue = userParameters->getAminoAcidCode(val);
                }
                if (userParameters->getLowercase())
                {
                    seq[j - firstRes] = (char)tolower((int)residue);
                }
                else
                {
                    seq[j - firstRes] = residue;
                }
                slen++;
            }
            for(j = 1; j <= slen; j++) 
            {
                (*gdeOutFile) << seq[j-1];
                if((j % lineLength == 0) || (j == slen)) 
                {
                    (*gdeOutFile) << "\n";
                }
            }
        }
  
        if (userParameters->getOutputStructPenalties() == 0 || 
            userParameters->getOutputStructPenalties() == 2) 
        {
            if (userParameters->getStructPenalties1() == SECST && 
                userParameters->getUseSS1() == true) 
            {
                (*gdeOutFile) << "\"SS_" << setw(alignPtr->getMaxNames()) << left
                              << alignPtr->getSecStructName1() << "\n";
                              
                //for(i = firstRes; i < firstRes + lastRes; i++) //- nige
                for(i = firstRes; i <= lastRes; i++) //- nige
                {
                    val = _ssMask1[i - 1];
                    if (val == userParameters->getGapPos1() || 
                        val == userParameters->getGapPos2())
                    {
                        seq[i - firstRes] = '-';
                    }
                    else
                    {
                        seq[i - firstRes] = val;
                    }
                }

                for(i = 1; i <= lastRes; i++) 
                {
                    (*gdeOutFile) << seq[i-1];
                    if((i % lineLength == 0) || (i == lastRes)) 
                    {
                        (*gdeOutFile) << "\n";
                    }
                }
            }
    
            if (userParameters->getStructPenalties2() == SECST && 
                userParameters->getUseSS2() == true) 
            {
                (*gdeOutFile) << "\"SS_" << setw(alignPtr->getMaxNames()) << left
                              << alignPtr->getSecStructName2() << "\n";
                                              
                //for(i = firstRes; i < firstRes + lastRes; i++) //- nige
                for(i = firstRes; i <= lastRes; i++) //- nige
                {
                    val = _ssMask2[i - 1];
                    if (val == userParameters->getGapPos1() || 
                        val == userParameters->getGapPos2())
                    {
                        seq[i - firstRes] = '-';
                    }
                    else
                    {
                        seq[i - firstRes] = val;
                    }
                }

                //for(i = 1; i <= lastRes; i++) //- nige
                for(i = 1; i <= length; i++) //- nige
                {
                    (*gdeOutFile) << seq[i - 1];
                    //if((i % lineLength == 0) || (i == lastRes)) //- nige
                    if((i % lineLength == 0) || (i == length)) //- nige
                    {
                        (*gdeOutFile) << "\n";
                    }
                }
            }
        }
        if (userParameters->getOutputStructPenalties() == 1 || 
            userParameters->getOutputStructPenalties() == 2) 
        {
            if (userParameters->getStructPenalties1() != NONE && 
                userParameters->getUseSS1() == true) 
            {
                (*gdeOutFile) << "\"GM_" << setw(alignPtr->getMaxNames()) << left
                              << alignPtr->getSecStructName1() << "\n";
                
                vector<char>* _gapPenaltyMask1 = alignPtr->getGapPenaltyMask1();          
                //for(i = firstRes; i < firstRes + lastRes; i++) //- nige
                for(i = firstRes; i <= lastRes; i++) //- nige
                {
                    val = _gapPenaltyMask1->at(i - 1);
                    if (val == userParameters->getGapPos1() || 
                        val == userParameters->getGapPos2())
                    {
                        seq[i - firstRes] = '-';
                    }
                    else
                    {
                        seq[i - firstRes] = val;
                    }
                }

                //for(i = 1; i <= lastRes; i++) //- nige
                for(i = 1; i <= length; i++) //- nige
                {
                    (*gdeOutFile) << seq[i - 1];
                    //if((i % lineLength == 0) || (i == lastRes)) //- nige
                    if((i % lineLength == 0) || (i == length)) //- nige
                    {
                        (*gdeOutFile) << "\n";
                    }
                }
            }
            if (userParameters->getStructPenalties2() != NONE && 
                userParameters->getUseSS2() == true) 
            {
                (*gdeOutFile) << "\"GM_" << setw(alignPtr->getMaxNames()) << left
                              << alignPtr->getSecStructName2() << "\n";
                
                vector<char>* _gapPenaltyMask2 = alignPtr->getGapPenaltyMask2();    
                //for(i = firstRes; i < firstRes + lastRes; i++) //- nige
                for(i = firstRes; i < length; i++) //- nige
                {
                    val = _gapPenaltyMask2->at(i-1);
                    if (val == userParameters->getGapPos1() || 
                        val == userParameters->getGapPos2())
                    {
                        seq[i - firstRes] = '-';
                    }
                    else
                    {
                        seq[i - firstRes] = val;
                    }
                }

                //for(i = 1; i <= lastRes; i++) //- nige
                for(i = 1; i <= length; i++) //- nige
                {
                    (*gdeOutFile) << seq[i - 1];
                    //if((i % lineLength == 0) || (i == lastRes)) //- nige
                    if((i % lineLength == 0) || (i == length)) //- nige
                    {
                        (*gdeOutFile) << "\n";
                    }
                }
            }
        }
        gdeOutFile->close();   
    }
    catch(const bad_alloc& e)
    {
        gdeOutFile->close();
        cerr << "A bad_alloc exception has occured in the gdeOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(VectorOutOfRange e)
    {
        gdeOutFile->close();
        cerr << "An exception has occured in the gdeOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(...)
    {
        gdeOutFile->close();
        cerr << "An exception has occured in the gdeOut function.\n";
        throw 1;
    }     
}

void AlignmentOutput::nbrfOut(Alignment* alignPtr, outputRegion partToOutput)
{
    char residue;
    int val;
    int i, ii;
    int j, slen;    
    int lineLength;
    rangeNum rnum;
    int firstRes = partToOutput._firstRes;
    int lastRes = partToOutput._lastRes;
    int firstSeq = partToOutput._firstSeq;
    int lastSeq = partToOutput._lastSeq;
            
    try
    {
        int _maxLength = alignPtr->getMaxAlnLength();
        vector<char> sequence;
        sequence.assign(_maxLength + 1, '0');
        const SeqArray* alignment = alignPtr->getSeqArray();
        // decide the line length for this alignment - maximum is LINELENGTH 
        lineLength = PAGEWIDTH - alignPtr->getMaxNames();
        lineLength = lineLength - lineLength % 10; // round to a multiple of 10
        if (lineLength > LINELENGTH || lineLength <= 0) // Mark, 18-7-07
        {
            lineLength = LINELENGTH;
        }
        
        // Get name prefix. DL if DNA, P1 if protein
        string namePrefix;
        bool _dnaFlag = userParameters->getDNAFlag();
        namePrefix = _dnaFlag ? ">DL;" : ">P1;";
        
        //int length = lastRes - firstRes + 1; //- nige

        for(ii = firstSeq; ii <= lastSeq; ii++) 
        {
            i = alignPtr->getOutputIndex(ii - 1);
            
            (*nbrfOutFile) << namePrefix;
            
            if (!userParameters->getSeqRange()) 
            {
                (*nbrfOutFile) << alignPtr->getName(i) << "\n"
                               << alignPtr->getTitle(i) << "\n";
            }
            else 
            {
                findRangeValues(alignPtr, &rnum, firstRes, lastRes, ii);
                (*nbrfOutFile) << nameonly(alignPtr->getName(i)) << "/" << rnum.start
                               << "-"<< rnum.end << "\n"
                               << alignPtr->getTitle(i) << "\n";
            }
            slen = 0;
            //for(j = firstRes; j < firstRes + lastRes; j++) //- nige
            for(j = firstRes; j <= lastRes; j++) //- nige
            {
                // NOTE I changed this here!!!!!
                if (j <= alignPtr->getSeqLength(i))
                {
                    //val = alignPtr.getResidue(i, j);
                    val = (*alignment)[i][j]; // NOTE june29
                }
                else
                { 
                    val = -3;
                }
                if((val == -3) || (val == 253))
                {
                    break;
                }
                else if((val < 0) || (val > userParameters->getMaxAA()))
                {
                    residue = '-';
                }
                else 
                {
                    residue = userParameters->getAminoAcidCode(val);
                }
                sequence[j - firstRes] = residue;
                slen++;
            }
            for(j = 1; j <= slen; j++) 
            {
                (*nbrfOutFile) << sequence[j - 1];
                if((j % lineLength == 0) || (j == slen))
                { 
                    (*nbrfOutFile) << "\n";
                }
            }
            (*nbrfOutFile) << "*\n";
        }    
        nbrfOutFile->close();
    }
    catch(const bad_alloc& e)
    {
        nbrfOutFile->close();
        cerr << "A bad_alloc exception has occured in the nbrfOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(VectorOutOfRange e)
    {
        nbrfOutFile->close();
        cerr << "An exception has occured in the nbrfOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(...)
    {
        nbrfOutFile->close();
        cerr << "An exception has occured in the nbrfOut function.\n";
        throw 1;
    }    
}

/* 
 * The function clustalOut is used to ouput the Alignment into a clustal format
 * file.
 */
void AlignmentOutput::clustalOut(Alignment* alignPtr, outputRegion partToOutput, ClustalWOutput *output)
{
	stringstream line;
    int firstRes = partToOutput._firstRes;
    int lastRes = partToOutput._lastRes;
    int firstSeq = partToOutput._firstSeq;
    int lastSeq = partToOutput._lastSeq;
        
    try
    {
        vector<char> seq1;    
        vector<int> seqNo;
        vector<int> printSeqNo;    
        vector<char> ss_mask1, ss_mask2;
        const SeqArray* alignment = alignPtr->getSeqArray();
        char  temp[MAXLINE];
        char c;
        int val;
        int ii, lv1, catident1[NUMRES], catident2[NUMRES], ident, chunks;
        int i, j, k, l;
        int pos, ptr;
        int lineLength;

        rangeNum rnum;
        string tmpStr = "";
        string sequenceLine = "";


        int _numSequences = alignPtr->getNumSeqs();
        // PMcG revert to Clustalw1.83 output format for name width
        const int NAMESWIDTH=std::max(std::min(MAXNAMESTODISPLAY,alignPtr->getMaxNames()) , MINNAMESTODISPLAY); 

        seqNo.assign(_numSequences + 1, 0);
        printSeqNo.assign(_numSequences + 1, 0);
    
        // check that lastSeq <= _numSequences
        if(lastSeq > _numSequences)
        {
            lastSeq = _numSequences;
        }
    
        for (i = firstSeq; i <= lastSeq; i++)
        {
            printSeqNo[i] = seqNo[i] = 0;
            for(j = 1; j < firstRes; j++) 
            {
                //val = alignPtr.getResidue(i, j);
                val = (*alignment)[i][j]; // NOTE june29
                if((val >= 0) || (val <= userParameters->getMaxAA()))
                { 
                    seqNo[i]++;
                }
            }
        }

        seq1.assign(alignPtr->getMaxAlnLength() + 1, 0);
        // Check if we have secondary structure in file 1 and if we want to output it.
        if (userParameters->getStructPenalties1() == SECST && 
            userParameters->getUseSS1() == true) 
        {
            int lengthSeq1 = alignPtr->getSeqLength(1);
            ss_mask1.assign(lengthSeq1 + 10, 0);
            vector<char>* _secStructMask1 = alignPtr->getSecStructMask1();
        
            for (i = 0; i < lengthSeq1; i++)
            {
                ss_mask1[i] = _secStructMask1->at(i);
            }
            printSecStructMask(lengthSeq1, _secStructMask1, &ss_mask1);
          
        }
        
        // Check if we have secondary structure in file 2 and if we want to output it.
        if (userParameters->getStructPenalties2() == SECST && 
            userParameters->getUseSS2() == true &&
            userParameters->getProfile2Empty() == false)
        {
            // AW added test getProfile2Empty was meant to fix bug
            // 141, but only prevents segfault, output still faulty.
            // Has to be fixed somewhere else
            int indexProfile2FirstSeq = alignPtr->getProfile1NumSeqs() + 1;
            int lengthSeqProfile2 = alignPtr->getSeqLength(indexProfile2FirstSeq);
            ss_mask2.assign(lengthSeqProfile2 + 10, 0);
            vector<char>* _secStructMask2 = alignPtr->getSecStructMask2();
            for (i = 0; i < lengthSeqProfile2; i++)
            {
                ss_mask2[i] = _secStructMask2->at(i);
            }
            printSecStructMask(lengthSeqProfile2, _secStructMask2, &ss_mask2);
        }
    
        //(*clustalOutFile) << "CLUSTAL "<< userParameters->getRevisionLevel()
        //                  << " multiple sequence alignment\n\n";
        line << "CLUSTAL "<< userParameters->getRevisionLevel()
			 << " ";
        output->msa.push_back(line.str());
        output->msa.push_back("");
        // decide the line length for this alignment - maximum is LINELENGTH

        //PMcG  9-2-2008 make line output same as 1.83
        //lineLength = PAGEWIDTH - alignPtr->getMaxNames();
        lineLength = PAGEWIDTH - NAMESWIDTH;
        lineLength = lineLength - lineLength % 10; // round to a multiple of 10
        if (lineLength > LINELENGTH || lineLength <= 0) // Mark 18-7-07
        { 
            lineLength = LINELENGTH;
        }
    
        int length = lastRes - firstRes + 1; //- nige

        chunks = length / lineLength; //- nige
        if(length % lineLength != 0) //- nige
        {
          ++chunks;
        }
        //printf("firstRes=%d,lastRes=%d,length=%d,chunks=%d\n",
        //       firstRes, lastRes, length, chunks);

        // This will loop through each of the blocks.
        line.str("");
        for(lv1 = 1; lv1 <= chunks; ++lv1) 
        {
            // pos is begining of chunk, ptr is the end of the chunk to be displayed.
            pos = ((lv1 - 1) * lineLength) + 1; 
            ptr = (length < pos + lineLength - 1) ? length : pos + lineLength - 1; //- nige
   
            //(*clustalOutFile) << "\n";
            output->msa.push_back(line.str());
            int _outStructPenalty = userParameters->getOutputStructPenalties();
            
            string secStructName1 = alignPtr->getSecStructName1(); // Mark 18-7-07
            
            if((int)secStructName1.size() > MAXNAMESTODISPLAY)
            {
                //secStructName1 = secStructName1.substr(0, MAXNAMESTODISPLAY);
                secStructName1 = secStructName1.substr(0, NAMESWIDTH);
            }            
            if (_outStructPenalty == 0 || _outStructPenalty == 2) 
            {
                if (userParameters->getStructPenalties1() == SECST && 
                    userParameters->getUseSS1() == true) 
                {
                    for(i = pos; i <= ptr; ++i) 
                    {
                        // Check if we can access mask position first
                        if((int)ss_mask1.size() > i + firstRes - 2)  
                        {
                            val = ss_mask1[i + firstRes - 2]; 
                            if (val == userParameters->getGapPos1() || 
                                val == userParameters->getGapPos2())
                            {
                                temp[i - pos] = '-';
                            }
                            else
                            {
                                temp[i - pos] = val;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                    temp[i - pos] = EOS;  
                                      
                    if(userParameters->getSeqRange())
                    {
                        //(*clustalOutFile) << "!SS_" << setw(MAXNAMESTODISPLAY +
                        //(*clustalOutFile) << "!SS_" << setw(NAMESWIDTH +
                        //                                    clusSecStructOffset)
                        //                  << left << secStructName1
                        //                  << "  " << temp << "\n";
                        line.str("");
                        line << "!SS_" << setw(NAMESWIDTH + clusSecStructOffset)
                        	 << left << secStructName1 << "  " << temp;
                        output->msa.push_back(line.str());
                    }
                    else
                    {
                        //(*clustalOutFile) << "!SS_" << setw(MAXNAMESTODISPLAY)
                        //(*clustalOutFile) << "!SS_" << setw(NAMESWIDTH)
                        //                  << left << secStructName1 << "  "
                        //                  << temp << "\n";
                        line.str("");
                        line << "!SS_" << setw(NAMESWIDTH)
                             << left << secStructName1 << "  "
                             << temp;
                        output->msa.push_back(line.str());
                    }
                }
            }
            if (_outStructPenalty == 1 || _outStructPenalty == 2) 
            {
                if (userParameters->getStructPenalties1() != NONE && 
                    userParameters->getUseSS1() == true) 
                {
                    vector<char>* _gapPenaltyMask1 = alignPtr->getGapPenaltyMask1();       
                    for(i = pos; i <= ptr; ++i) 
                    {
                        if((int)_gapPenaltyMask1->size() > i + firstRes - 2) 
                        {
                            val = _gapPenaltyMask1->at(i + firstRes - 2); 
                            if (val == userParameters->getGapPos1() || 
                                val == userParameters->getGapPos2())
                            {
                                temp[i - pos] = '-';
                            }
                            else
                            {
                                temp[i - pos] = val;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                    temp[i - pos] = EOS;
                    
                    //(*clustalOutFile) << "!GM_" << setw(MAXNAMESTODISPLAY) << left
                    //(*clustalOutFile) << "!GM_" << setw(NAMESWIDTH) << left
                    //                  << secStructName1 << "  "
                    //                  << temp << "\n";
                    line.str("");
					line << "!GM_" << setw(NAMESWIDTH)
						 << left << secStructName1 << "  "
						 << temp;
					output->msa.push_back(line.str());
                }
            }
            
            
            string secStructName2 = alignPtr->getSecStructName2(); // Mark 18-7-07
            //if((int)secStructName2.size() > MAXNAMESTODISPLAY)
            if((int)secStructName2.size() > NAMESWIDTH)
            {
                //secStructName2 = secStructName2.substr(0, MAXNAMESTODISPLAY);
                secStructName2 = secStructName2.substr(0, NAMESWIDTH);
            }
                                
            if (_outStructPenalty == 0 || _outStructPenalty == 2) 
            {
                if (userParameters->getStructPenalties2() == SECST && 
                    userParameters->getUseSS2() == true) 
                {
                    for(i = pos; i <= ptr; ++i) 
                    {
                        if((int)ss_mask2.size() > i + firstRes - 2)
                        {
                            val = ss_mask2[i + firstRes - 2]; 
                            if (val == userParameters->getGapPos1() || 
                                val == userParameters->getGapPos2())
                            {
                                temp[i - pos] = '-';
                            }
                            else
                            {
                                temp[i - pos] = val;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                    temp[i - pos] = EOS;
                    
                                        
                    if (userParameters->getSeqRange())
                    {
                        //(*clustalOutFile) << "!SS_" << setw(MAXNAMESTODISPLAY +
                        //(*clustalOutFile) << "!SS_" << setw(NAMESWIDTH +
                        //                                    clusSecStructOffset)
                        //                  << left << secStructName2;
                        //(*clustalOutFile) << "  " << temp << "\n";
                        line.str("");
						line << "!SS_" << setw(NAMESWIDTH + clusSecStructOffset)
							 << left << secStructName2 << "  "
							 << temp;
						output->msa.push_back(line.str());
                    }
                    else
                    {
                        //(*clustalOutFile) << "!SS_" << setw(MAXNAMESTODISPLAY)
                        //(*clustalOutFile) << "!SS_" << setw(NAMESWIDTH)
                        //                  << left << secStructName2 << "  "
                        //                  << temp << "\n";
                        line.str("");
						line << "!SS_" << setw(NAMESWIDTH)
							 << left << secStructName2 << "  "
							 << temp;
						output->msa.push_back(line.str());
                    }
                }
            }
            
            if (_outStructPenalty == 1 || _outStructPenalty == 2) 
            {
                if (userParameters->getStructPenalties2() != NONE && 
                    userParameters->getUseSS2() == true) 
                {
                    vector<char>* _gapPenaltyMask2 = alignPtr->getGapPenaltyMask2();
                    for(i = pos; i <= ptr; ++i) 
                    {
                        if((int)_gapPenaltyMask2->size() > i + firstRes - 2)
                        {
                            val = _gapPenaltyMask2->at(i + firstRes - 2); 
                            if (val == userParameters->getGapPos1() || 
                                val == userParameters->getGapPos2())
                            {
                                temp[i - pos] = '-';
                            }
                            else
                            {
                                temp[i - pos] = val;
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                    temp[i - pos] = EOS;
                    
                    //(*clustalOutFile) << "!GM_" << setw(MAXNAMESTODISPLAY) << left
                    //(*clustalOutFile) << "!GM_" << setw(NAMESWIDTH) << left
                    //                  << secStructName2 << "  "
                    //                  << temp << "\n";
                    line.str("");
					line << "!GM_" << setw(NAMESWIDTH)
						 << left << secStructName2 << "  "
						 << temp;
					output->msa.push_back(line.str());
                }
            }
        
            for(ii = firstSeq; ii <= lastSeq; ++ii) 
            {
                i = alignPtr->getOutputIndex(ii - 1); 
                printSeqNo[i] = 0;
                for(j = pos; j <= ptr; ++j) 
                {
                    if (j + firstRes - 1 <= alignPtr->getSeqLength(i))
                    {
                        //val = alignPtr.getResidue(i, j + firstRes - 1);
                        val = (*alignment)[i][j + firstRes - 1]; // NOTE june29
                    }
                    else
                    { 
                        val = -3;
                    }
                    if((val == -3) || (val == 253)) 
                    {
                        break;
                    }
                    else if((val < 0) || (val > userParameters->getMaxAA()))
                    {
                        seq1[j] = '-';
                    }
                    else 
                    {
                        seq1[j] = userParameters->getAminoAcidCode(val);
                        seqNo[i]++;
                        printSeqNo[i] = 1;
                    } 
                }
                for(; j <= ptr; ++j) 
                {
                    seq1[j]='-';
                }
            
                sequenceLine = "";

                for(int index = pos; index < ptr + 1; index++)
                {
                    sequenceLine += seq1[index];
                }
                
                string seqName = alignPtr->getName(i);
                //if((int)seqName.size() > MAXNAMESTODISPLAY)
                if((int)seqName.size() > NAMESWIDTH)
                {
                    //seqName = seqName.substr(0, MAXNAMESTODISPLAY);
                    seqName = seqName.substr(0, NAMESWIDTH);
                }
                
                if (!userParameters->getSeqRange()) 
                {
                    // NOTE I made a change here from + 5, to + 6.
                    //(*clustalOutFile) << setw(alignPtr->getMaxNames() + 6) << left
                    //                  << alignPtr->getName(i);

                    //(*clustalOutFile) << setw(MAXNAMESTODISPLAY + 6) << left
                    //                  << seqName; 
                    // PMcG changed this to revert to behaviour of clustalw1.83
                    // 
                    //(*clustalOutFile) << setw(std::max(std::min(MAXNAMESTODISPLAY,alignPtr->getMaxNames()) , MINNAMESTODISPLAY) + 6)
                    //(*clustalOutFile) << setw(NAMESWIDTH + 6) << left << seqName;
                    line.str("");
                    line << setw(NAMESWIDTH + 6) << left << seqName;


                }
                else // Put the sequence range information in! 
                {
                    findRangeValues(alignPtr, &rnum, firstRes, lastRes, ii);
                    string nameOnly = nameonly(seqName);
                    std::stringstream ss;
                    std::stringstream ss2;
                    string rangeStart;
                    string rangeEnd;
                    ss << rnum.start;
                    ss >> rangeStart;
                    ss2 << rnum.end;
                    ss2 >> rangeEnd;
                    tmpStr = nameOnly + "/" + rangeStart + "-" + rangeEnd;
                    //(*clustalOutFile) << setw(MAXNAMESTODISPLAY + clusSequenceOffset) 
                    //(*clustalOutFile) << setw(NAMESWIDTH + clusSequenceOffset)
                    //                  << left << tmpStr;
                    line.str("");
                    line << setw(NAMESWIDTH + clusSequenceOffset) << left << tmpStr;
                }
            
                //(*clustalOutFile) << sequenceLine;
                line << sequenceLine;

                if (userParameters->getClSeqNumbers() && printSeqNo[i])
                {
                    //(*clustalOutFile) << " " << seqNo[i];
                	line << " " << seqNo[i];
                }
            
                //(*clustalOutFile) << "\n";
                output->msa.push_back(line.str());
                line.str("");
            }
        
            // Now print out the conservation information!
            for(i = pos; i <= ptr; ++i) 
            {
                seq1[i] = ' ';
                ident = 0;

                for(j = 1; j <= (int)strongGroup.size(); j++)
                {
                    catident1[j - 1] = 0;
                }
                for(j = 1; j <= (int)weakGroup.size(); j++) 
                {
                    catident2[j - 1] = 0;
                }
                
                for(j = firstSeq; j <= lastSeq; ++j) 
                {
                    if((i + firstRes - 1 <= alignPtr->getSeqLength(j)) &&
                       (i + firstRes - 1 <= alignPtr->getSeqLength(firstSeq)))
                    {// NOTE june29
                        if(((*alignment)[firstSeq][i + firstRes - 1] >= 0) && 
                            ((*alignment)[firstSeq][i + firstRes - 1] <=
                            userParameters->getMaxAA())) 
                        {
                            // Count how many are identical to the first sequence
                            if((*alignment)[firstSeq][i + firstRes - 1] == 
                                (*alignment)[j][i + firstRes - 1])
                            {
                                ++ident;
                            }
                            // Count how many are in the same category.
                            for(k = 1; k <= (int)strongGroup.size(); k++) 
                            {
                                for(l = 0; (c = strongGroup[k - 1][l]); l++) 
                                {
                                    int resCode = (*alignment)[j][i + firstRes - 1];
                                    if(resCode <= userParameters->getMaxAA() + 1)
                                    {
                                        if (userParameters->getAminoAcidCode(resCode) == c)
                                        {
                                            catident1[k - 1]++;
                                            break;
                                        }
                                    }
                                }
                            }
                            for(k = 1; k <= (int)weakGroup.size(); k++) 
                            {
                                for(l = 0; (c = weakGroup[k - 1][l]); l++) 
                                {//NOTE june29
                                    int resCode = (*alignment)[j][i + firstRes - 1];
                                    if(resCode <= userParameters->getMaxAA() + 1)
                                    {
                                        if (userParameters->getAminoAcidCode(resCode) == c)
                                        {
                                            catident2[k - 1]++;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }    
                }
            
                // Now do the conservation part for each block. 
                if(ident == lastSeq - firstSeq + 1)
                {
                    seq1[i] = '*'; // All residues the same!
                }
                else if (!userParameters->getDNAFlag()) 
                {
                    for(k = 1; k <= (int)strongGroup.size(); k++) 
                    {
                        if (catident1[k - 1] == lastSeq - firstSeq + 1) 
                        {
                            seq1[i] = ':'; // All residues member of the same category
                            break;
                        }
                    }
                    if(seq1[i] == ' ')
                    {
                        for(k = 1; k <= (int)weakGroup.size(); k++) 
                        {
                            if (catident2[k - 1] == lastSeq - firstSeq + 1) 
                            {
                                seq1[i] = '.'; // All residues member of the same category
                                break;
                            }
                        }
                    }
                }
            }

            sequenceLine = "";
            for(int index = pos; index < ptr + 1; index++)
            {
                sequenceLine += seq1[index];
            }

            //for(k = 0; k < MAXNAMESTODISPLAY + 6; k++)
            for(k = 0; k < NAMESWIDTH + 6; k++)
            { 
                //(*clustalOutFile) << " ";
            	line << " ";
            }
            if(userParameters->getSeqRange()) 
            {

                //(*clustalOutFile) << "          ";
            	line << "          ";
            }

            //(*clustalOutFile) << sequenceLine << "\n";
            line << sequenceLine;
            output->msa.push_back(line.str());
            line.str("");
        }
        //clustalOutFile->close();
    }
    catch(const bad_alloc& e)
    {
        //clustalOutFile->close();
        cerr << "A bad_alloc exception has occured in the clustalOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(VectorOutOfRange e)
    {
        //clustalOutFile->close();
        cerr << "An exception has occured in the clustalOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(exception e)
    {
        //clustalOutFile->close();
        cerr << "An exception has occured in the clustalOut function.\n"
             << e.what() << "\n";
        throw 1;
    }
    catch(...)
    {
        //clustalOutFile->close();
        cerr << "An exception has occured in the clustalOut function.\n";
        throw 1;
    }    
    return;
}

/*
 * The function printSecStructMask takes in the 2 mask vectors. It makes changes to 
 * mask and puts the result in structMask. This is the one that will be printed out.
 * NOTE I have had to use (*structMask)[i] syntax. This is not good. But there is no
 * way to do random access setting in vectors otherwise.
 */
void AlignmentOutput::printSecStructMask(int prfLength, vector<char>* mask, 
                                         vector<char>* structMask)
{
    int i, j;
//    calculate the gap penalty mask from the secondary structures

    i = 0;
    try
    {
        while (i < prfLength) 
        {
            if (tolower(mask->at(i)) == 'a' || mask->at(i) == '$') 
            {
                for (j = 0; j < userParameters->getHelixEndMinus(); j++) 
                {
                    if (i + j >= prfLength || (tolower(mask->at(i+j)) != 'a'
                                    && mask->at(i+j) != '$')) 
                    {    
                        break;
                    }
                    (*structMask)[i+j] = 'a';
                }
                i += j;
                while (tolower(mask->at(i)) == 'a' || mask->at(i) == '$') 
                {
                    if (i >= prfLength) 
                        break;
                    if (mask->at(i) == '$') 
                    {
                        (*structMask)[i] = 'A';
                        i++;
                        break;
                    }
                    else 
                        (*structMask)[i] = (*mask)[i];
                    i++;
                }
                for (j = 0; j < userParameters->getHelixEndMinus(); j++) 
                {
                    if ((i-j-1>=0) && (tolower(mask->at(i-j-1)) == 'a'
                                        || mask->at(i-j-1) == '$'))
                    {
                        (*structMask)[i - j - 1] = 'a';
                    }
                }
            }
            else if (tolower(mask->at(i)) == 'b' || mask->at(i) == '%') 
            {
                for (j = 0; j < userParameters->getStrandEndMinus(); j++) 
                {
                    if (i + j >= prfLength || (tolower(mask->at(i+j)) != 'b'
                                        && mask->at(i+j) != '%')) 
                    {
                        break;
                    }
                    (*structMask)[i+j] = 'b';
                }
                i += j;
                while (tolower(mask->at(i)) == 'b' || mask->at(i) == '%') 
                {
                    if (i >= prfLength) 
                        break;
                    if (mask->at(i) == '%') 
                    {
                        (*structMask)[i] = 'B';
                        i++;
                        break;
                    }
                    else
                    { 
                        (*structMask)[i] = mask->at(i);
                    }
                    i++;
                }
                for (j = 0; j < userParameters->getStrandEndMinus(); j++) 
                {
                    if ((i - j - 1 >= 0) && (tolower(mask->at(i - j - 1)) == 'b' || 
                        mask->at(i-j-1) == '%'))
                    {
                        (*structMask)[i-j-1] = 'b';
                    }
                }
            }
            else
            { 
                i++;
            }
        }
    }
    catch(const exception& e)
    {
        // Catch all the exceptions
        cerr << "There has been an exception in the function printSecStructMask\n"
             << e.what() << "\n";
        throw 1;
    }
}

/*
 * The function findRangeValues is used to find the range to be printed out.
 *
 *
 */
void AlignmentOutput::findRangeValues(Alignment* alignPtr, rangeNum *rnum, int firstRes, 
                                int len, int firstSeq)
{
    assert(rnum);
    assert(firstRes > 0);
    assert(len > 0);
    assert(firstSeq > 0);
    
    try
    {
        int val;
        int i, ii;
        int j, slen;    

        char tmpName[FILENAMELEN + 15];
        int iStart = 0;
        int iEnd = 0; // to print sequence start-end with names
        int found = 0;
        int ngaps = 0;
        int tmpStart = 0; 
        int tmpEnd = 0;
        int ntermgaps = 0;
        int pregaps = 0;
        int tmpk = 0;
        int isRange = 0;
        int formula = 0;
        const SeqArray* alignment = alignPtr->getSeqArray();
        
        tmpName[0] = '\0';
        slen = 0;

        ii = firstSeq ;
        i = alignPtr->getOutputIndex(ii - 1); // NOTE Same as elsewhere! 
        string name = alignPtr->getName(i);
        
        if( (sscanf(name.c_str(),"%[^/]/%d-%d",tmpName, &tmpStart, &tmpEnd) == 3)) 
        {
            isRange = 1;
        }
        
        for(tmpk = 1; tmpk < firstRes; tmpk++)
        { 
            // do this irrespective of above sscanf
            //val = alignPtr.getResidue(i, tmpk);
            val = (*alignment)[i][tmpk]; // NOTE june29
            if ((val < 0) || (val > userParameters->getMaxAA()))
            { 
                //it is gap
                pregaps++;
            }
        }
        
        for(j = firstRes; j < firstRes + len; j++) 
        {
            if(j > alignPtr->getSeqLength(i))
            {
                val = -3; // Cant get it so break out.
            }
            else
            {
                //val = alignPtr.getResidue(i, j);
                val = (*alignment)[i][j]; // NOTE june29
            }
            if((val == -3) || (val == 253))
            {
                break;
            }
            else if((val < 0) || (val > userParameters->getMaxAA())) 
            {
                ngaps++;
            }
            else 
            {
                found = j;
            }
            if ( found && (iStart == 0) ) 
            {
                iStart = found;
                ntermgaps = ngaps;
            }
            slen++;
        }
        if(userParameters->getSeqRange()) 
        {
            cout << "Name : " << alignPtr->getName(i) << " "
                 << "\n  firstRes = "<< firstRes << " "
                 << "   len = " << len << " "
                 << "\n  iStart = " << iStart << " "
                 << "\n  tmpStart = " << tmpStart << " "
                 << "\n  ngaps = " << ngaps << " "
                 << "\n  pregaps = " << pregaps << " ";
            if (!isRange)
            {
                formula = iStart - pregaps;
            }
            else
            {
                formula = iStart - pregaps +  ( tmpStart == 1 ? 0: tmpStart-1) ;
            }

            cout << "\n\nsuggestion  iStart - pregaps + tmpStart - ntermgaps = "
                 << iStart << " - " << pregaps << " + " << tmpStart << " - " << ntermgaps
                 << " formula " << formula << " ";
        }
        else 
        {
            cerr << "\n no range found .... strange,  iStart = " << iStart;
            formula = 1;
        }
        if (pregaps == firstRes - 1) // all gaps -  now the conditions........ 
        {
            formula = tmpStart ; // keep the previous start... 
        }
        formula = (formula <= 0) ? 1: formula;
        if (pregaps == 0 && tmpStart == 0) 
        {
            formula = firstRes;
        }
        iEnd = formula + len - ngaps -1;
        rnum->start = formula;
        rnum->end = iEnd;
        cout << "\n check... " << alignPtr->getName(i) << " " << rnum->start 
             << " - " << rnum->end;
        cout << " Done checking.........";
    }
    catch(VectorOutOfRange e)
    {
        cerr << "An exception has occured in the findRangeValues function.\n"
             << e.what() << "\n";
        throw 1;
    }    
    catch(...)
    {
        cerr << "An exception has occured in findRangeValues function\n";
        throw 1;
    }
}

string AlignmentOutput::nameonly(string s)
{
    string tmp;
    try
    {
        int i = 0;

        while (i < (int)s.size()) 
        {
            if(s.at(i) != '/')
            {
                tmp += s.at(i);
                i++;
            }
            else
            {
                break;
            }
        }
    }
    catch(const exception& e)
    {
        cerr << "An exception has occured in the function nameonly\n"
             << e.what();
        throw 1;
    }
    return tmp;
}

int AlignmentOutput::SeqGCGCheckSum(vector<char>* seq, int length)
{
    int i;
    long check;
    int seqResIndex = 1;
    //int seqLength = seq->size();
    
    for (i = 0, check = 0; i < length; i++, seqResIndex++)
    {
        char _toUpperCase = (*seq)[seqResIndex];
        check += ((i % 57) + 1) * toupper(_toUpperCase);
    }

    return (check % 10000);
}

void AlignmentOutput::showAlign()
{
    //FILE *file;
    int  numLines;
    char temp[MAXLINE + 1];
    temp[0] = '\0';
    string fileName;
    string answer;
    
    if(userParameters->getOutputClustal()) 
    {
        fileName = clustalOutName;
    }
    else if(userParameters->getOutputNbrf()) 
    {
        fileName = nbrfOutName;
    }
    else if(userParameters->getOutputGCG()) 
    {
        fileName = gcgOutName;
    }
    else if(userParameters->getOutputPhylip()) 
    {
        fileName = phylipOutName;
    }
    else if(userParameters->getOutputGde()) 
    {
        fileName = gdeOutName;
    }
    else if(userParameters->getOutputNexus()) 
    {
        fileName = nexusOutName;
    }
    else if(userParameters->getOutputFasta()) 
    {
        fileName = fastaOutName;
    }
    else
    {
        return; // No file output type!
    }

    ifstream _fileIn;
    _fileIn.open(fileName.c_str(), ios::in);        
    _fileIn.seekg(0, std::ios::beg); // start at the beginning
    
    cout << "\n\n";
    numLines = 0;

    while(_fileIn.getline(temp, MAXLINE + 1, '\n')) 
    {
       //fputs(temp,stdout);
       cout << temp << "\n";
       ++numLines;
       if(numLines >= PAGE_LEN) 
       {
           cout << "\n";
           utilityObject->getStr(string("Press [RETURN] to continue or  X  to stop"), answer);
           if(toupper(answer[0]) == 'X') 
           {
               _fileIn.close();
               return;
           }
           else
           {
               numLines = 0;
           }
       }
    }
    _fileIn.close();
    cout << "\n";
    utilityObject->getStr(string("Press [RETURN] to continue"), answer);
}

}

