/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/*
 * The class AlignmentOutput is used to output the Alignment in all the different
 * formats that have been selected. It will output all the different file types if
 * these have been selected from the menu or the commandline.
 * To use this class we must call openAlignmentOutput first. Then we call the function
 * createAlignmentOutput with an Alignment to be output and the first and last sequence
 * to be output as well. 
 */
#ifndef ALIGNMENTOUTPUT_H
#define ALIGNMENTOUTPUT_H

#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <sstream>
#include <exception>
#include <cassert>
#include "Alignment.h"

#include "../RClustalW.h"


namespace clustalw
{

typedef struct rangeNum 
{
    int start;
    int end;
} rangeNum;

typedef struct outputRegion 
{
    int _firstSeq;
    int _lastSeq;
    int _firstRes;
    int _lastRes;
} outputRegion;

class AlignmentOutput
{
    public:
        /* Functions */
        AlignmentOutput();
        bool openAlignmentOutput(string path);
        bool QTOpenFilesForOutput(AlignmentFileNames fileNames);
        void createAlignmentOutput(Alignment* alignPtr, int firstSeq, int lastSeq, ClustalWOutput *output);
        void printSecStructMask(int prfLength, vector<char>* mask, vector<char>* structMask);
        /* Attributes */

    private:
        /* Functions */
        void fastaOut(Alignment* alignPtr, outputRegion partToOutput, ClustalWOutput *output);
        void clustalOut(Alignment* alignPtr, outputRegion partToOutput, ClustalWOutput *output);
        void gcgOut(Alignment* alignPtr, outputRegion partToOutput);
        void nexusOut(Alignment* alignPtr, outputRegion partToOutput);
        void phylipOut(Alignment* alignPtr, outputRegion partToOutput);
        void nbrfOut(Alignment* alignPtr, outputRegion partToOutput);
        void gdeOut(Alignment* alignPtr, outputRegion partToOutput);
        string nameonly(string s);
        
        void findRangeValues(Alignment* alignPtr, rangeNum *rnum, int firstRes, int lastRes, 
                             int firstSeq);
        bool openExplicitFile(auto_ptr<ofstream>& outFile, string fileName);
        string openOutputFile(auto_ptr<ofstream>& outFile, string prompt, string path, 
                              string fileExtension);
        int SeqGCGCheckSum(vector<char>* sequence, int length);
        void showAlign();
        /* Attributes */
   
        auto_ptr<ofstream> clustalOutFile;
        auto_ptr<ofstream> gcgOutFile;
        auto_ptr<ofstream> nbrfOutFile;
        auto_ptr<ofstream> phylipOutFile;
        auto_ptr<ofstream> gdeOutFile;
        auto_ptr<ofstream> nexusOutFile;
        auto_ptr<ofstream> fastaOutFile;
        
        string clustalOutName;
        string gcgOutName;
        string phylipOutName;
        string nbrfOutName;
        string gdeOutName;
        string nexusOutName;
        string fastaOutName;
        vector<string> strongGroup; 
        vector<string> weakGroup;
        int clusSecStructOffset;
        int clusSequenceOffset;        
};

}
#endif

