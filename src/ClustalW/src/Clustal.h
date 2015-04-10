/**
 * Author: Mark Larkin
 *
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * The class Clustal is the main class in the program. It is used by the interactive
 * menu, command line parser and clustal x to perform the algorithmic part of the
 * program. 
 */
#ifndef CLUSTAL_H
#define CLUSTAL_H

#include <string>
#include "general/clustalw.h"
#include "general/utils.h"
#include "general/userparams.h"
#include "fileInput/FileReader.h"
#include "alignment/Alignment.h"
#include "alignment/AlignmentOutput.h"

using namespace std;

namespace clustalw
{

class Clustal
{
    public:
        /* Functions */
        Clustal();
        void align(string* phylipName, ClustalWOutput *output, bool createOutput = true);
        void sequencesAlignToProfile(string* phylipName, ClustalWOutput *output);
        void profileAlign(string* p1TreeName, string* p2TreeName, ClustalWOutput *output);
        void doGuideTreeOnly(string* phylipName);
        void doAlignUseOldTree(string* phylipName, ClustalWOutput *output);
        void getHelp(string helpPointer, bool printTitle = false);
        void getHelp(char helpPointer, bool printTitle = false);
        void getFullHelp();
        int sequenceInput(bool append, string *offendingSeq);
        int profile1Input(string profile1Name = "");
        int profile2Input(string profile2Name = "");
        int commandLineReadSeq(int firstSeq, ClustalWInput *input);
        void outputNow(ClustalWOutput *output);
        void phylogeneticTree(string* phylip_name, string* clustal_name, string* dist_name,
                              string* nexus_name, string pimName);
        void bootstrapTree(string* phylip_name, string* clustal_name, string* nexus_name);
        Alignment* getAlignmentPtr(){return &alignmentObj;} 
        void QTcalcLowScoreSegments(LowScoreSegParams* params);
        void QTcalcWeightsForLowScoreSeg(LowScoreSegParams* params);
        void QTremoveShortSegments(LowScoreSegParams* params);
        void QTSetFileNamesForOutput(AlignmentFileNames fileNames);
        bool QTRealignSelectedRange(AlignmentFileNames fileNames, int beginPos, int endPos,
                                    bool realignEndGapPen, ClustalWOutput *output);
        void test();
        /* Attributes */

    private:
        /* Functions */
        void initInterface();
        void calcGapPenaltyMask(int prfLength, vector<char>* mask, vector<char>* gapMask);
        bool useExistingGuideTree(int type, string* phylipName, const string& path);
        void promptForNewGuideTreeName(int type, string* treeName, const string& path);
        //bool removeFirstIterate(Alignment* alnPtr, DistMatrix* distMat);
        /* Attributes */
        enum{Sequences, Profile1, Profile2};
        string sequencesMsg, profile1Msg, profile2Msg;
        string newProfile1TreePrompt, newProfile2TreePrompt;
        
        Alignment alignmentObj;
        string helpFileName;
        int newSeq;
        bool checkTree;
        AlignmentFileNames QTFileNames;
};
}
#endif
