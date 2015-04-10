/**
 * Author: Mark Larkin
 *
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.
 */
/**
 * Changes:
 *
 * 10-02-07,Nigel Brown(EMBL): changed ifstream to InFileStream to handle
 * cross-platform end-of-lines (for dendrograms, not for help file). Remerged
 * this change 13-04-07.
 *
 * Mark 10-5-2007: Bug fix # 42. call getWeightsForQtLowScore function in
 * QTcalcWeightsForLowScoreSeg instead of getWeightsFromDistMat.
 *
 * Mark 22-5-07, changed the distmatrix to be the size of alignObject.numSeqs
 * Mark Change 20-6-07, added call to calculateMaxLengths()
 * Mark, 3-7-07, Changed getHelp.
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <cstdio>
#include "Clustal.h"
#include "pairwise/FullPairwiseAlign.h"
#include "pairwise/FastPairwiseAlign.h"
#include "multipleAlign/MSA.h"
#include "multipleAlign/LowScoreSegProfile.h"
#include "multipleAlign/Iteration.h"
#include "tree/TreeInterface.h"
#include "general/debuglogObject.h"
#include "general/statsObject.h"
#include "alignment/ObjectiveScore.h"
#include "general/ClustalWResources.h"
#include "Help.h"
#include <ctime>

using namespace std;

namespace clustalw
{

Clustal::Clustal()
{
#ifdef WINDOWS
    helpFileName = string("clustalw.hlp");
#else
    helpFileName = string("clustalw_help");
#endif

    checkTree = true;
    newSeq = 0;
    sequencesMsg = "\nUse the existing GUIDE TREE file,  ";
    profile1Msg = "\nUse the existing GUIDE TREE file for Profile 1,  " ;
    profile2Msg = "\nUse the existing GUIDE TREE file for Profile 2,  ";

    newProfile1TreePrompt = "\nEnter name for new GUIDE TREE file for profile 1 [";
    newProfile2TreePrompt = "\nEnter name for new GUIDE TREE file for profile 2 [";

    initInterface();
}

/** ***********************************************************************
 * The function align is used to do a full multiple sequence alignment.   *
 **************************************************************************/
// FIXME: merge doAlignUseOldTree in here
void Clustal::align(string* phylipName, ClustalWOutput *output, bool createOutput)
{
    //time_t start, end;
    //double dif;
    //start = time (NULL);
    //ObjectiveScore score;
    //double _score = score.getSagaScore(&alignmentObj);
    //cout << "SAGA score " << _score << "\n";

    string path;
    int count;
    AlignmentOutput alignOutput;

    if(userParameters->getEmpty() && userParameters->getMenuFlag())
    {
        utilityObject->error("No sequences in memory. Load sequences first.");
        return;
    }

    userParameters->setStructPenalties1(NONE);
    userParameters->setStructPenalties2(NONE);

    alignmentObj.clearSecStruct1();
    alignmentObj.clearSecStruct2();

    utilityObject->getPath(userParameters->getSeqName(), &path);

    if(createOutput && (userParameters->getMenuFlag() || !userParameters->getInteractive()))
    {
        if(!alignOutput.openAlignmentOutput(path))
        {
            return;
        }
    }
    else if(createOutput)
    {
        // We are using clustalQT.
        // Open all the files.
        if(!alignOutput.QTOpenFilesForOutput(QTFileNames))
        {
            return; // could not open the files.
        }
    }

    if (userParameters->getSaveParameters())
    {
        userParameters->createParameterOutput();
    }

    if(userParameters->getResetAlignmentsNew() || userParameters->getResetAlignmentsAll())
    {
        alignmentObj.resetAlign();
    }
    if(userParameters->getDisplayInfo())
    {
        cout << "Start of Pairwise alignments\n";
        cout << "Aligning...\n";
    }
    if(userParameters->getDNAFlag())
    {
        userParameters->setDNAParams();
    }
    else
    {
        userParameters->setProtParams();
    }

    if (statsObject->isEnabled()) {
        statsObject->logInputSeqStats(&alignmentObj);
    }


    /// STEP 1: PAIRWISE ALIGNMENT TO GENERATE DISTANCE MATRIX

    SymMatrix distMat;
    int _numSeqs = alignmentObj.getNumSeqs();
    distMat.ResizeRect(_numSeqs + 1);

    PairwiseAlignBase* pairwiseDist;
    if (userParameters->getQuickPairAlign())
    {
        pairwiseDist = new FastPairwiseAlign();
    }
    else
    {
        pairwiseDist = new FullPairwiseAlign();
    }
    // Generate distance matrix!
    pairwiseDist->pairwiseAlign(&alignmentObj, &distMat, 0, _numSeqs, 0, _numSeqs);
    delete pairwiseDist;

    bool success = false;
    auto_ptr<AlignmentSteps> progSteps;
    vector<int> seqWeight(_numSeqs + 1);

    #if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("Calling getWeightsAndStepsFromDistMat\n");
        }
    #endif

    TreeInterface calcSteps;
    progSteps = calcSteps.getWeightsAndStepsFromDistMat(&seqWeight, &distMat, &alignmentObj,
                                                        1, _numSeqs, phylipName, &success);
    //cout << "weights and steps calculated!\n";
    //end = time (NULL);
    //dif = difftime(end, start);
    //cout << "It took " << dif << " seconds so Far\n";

    if(!success)
    {
        #if DEBUGFULL
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("Unsuccessful!!!\n");
            }
        #endif
        return;
    }
    #if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("doing multiSeqAlign\n");
        }
    #endif
    MSA* msaObj = new MSA();
    count = msaObj->multiSeqAlign(&alignmentObj, &distMat, &seqWeight, progSteps.get(), 0);
    delete msaObj;
    //cout << "alignment finished!\n";

    //end = time (NULL);
    //dif = difftime(end, start);
    //cout << "It took " << dif << " seconds so Far\n";

    if (count <= 0)
    {
        return;
    }

    if (userParameters->getMenuFlag())
    {
        cout << "\n\n\n";
    }

    /// Do iteration to improve alignment!!!
    if(userParameters->getDoRemoveFirstIteration() == ALIGNMENT)
    {
        //userParameters->setNumIterations(_numSeqs * 2);
        Iteration iterateObj;
        iterateObj.removeFirstIterate(&alignmentObj);
        alignmentObj.calculateMaxLengths(); // Mark Change 20-6-07
        if(userParameters->getDisplayInfo())
            cout << "Finished iteration\n";
    }

    if (statsObject->isEnabled()) {
        statsObject->logAlignedSeqStats(&alignmentObj);
    }



    /// STEP 4: OUTPUT THE ALIGNMENT
    //if(createOutput)
    {
        alignOutput.createAlignmentOutput(&alignmentObj, 1, _numSeqs, output);
    }
    (*phylipName) = "";

    //end = time (NULL);
    //dif = difftime(end, start);
    //cout << "It took " << dif << " seconds\n";
    return;
}

/** ****************************************************************************
 * The function sequencesAlignToProfile is used to align a set of sequences to *
 * a profile                                                                   *
 *******************************************************************************/
void Clustal::sequencesAlignToProfile(string* phylipName, ClustalWOutput *output)
{

	cout << "sequencesAlignToProfile called";

    string path;
    string treeName;
    bool useTree;
    int i, j, count;
    float dscore;
    bool saveSS2;
    AlignmentOutput alignOutput;

    if(userParameters->getProfile1Empty() && userParameters->getMenuFlag())
    {
        utilityObject->error("No profile in memory. Input 1st profile first.\n");
        return;
    }

    if(userParameters->getProfile2Empty() && userParameters->getMenuFlag())
    {
        utilityObject->error("No sequences in memory. Input sequences first.\n");
        return;
    }

    utilityObject->getPath(userParameters->getProfile2Name(), &path);

    if(userParameters->getMenuFlag() || !userParameters->getInteractive())
    {
        if(!alignOutput.openAlignmentOutput(path))
        {
            return;
        }
    }
    else
    {
        // We are using clustalQT.
        // Open all the files.
        if(!alignOutput.QTOpenFilesForOutput(QTFileNames))
        {
            return; // could not open the files.
        }
    }

    newSeq = alignmentObj.getProfile1NumSeqs() + 1;

    // check for secondary structure information for list of sequences

    saveSS2 = userParameters->getUseSS2();

    if (userParameters->getStructPenalties2() != NONE && userParameters->getUseSS2() == true
        && (alignmentObj.getNumSeqs() - alignmentObj.getProfile1NumSeqs() > 1))
    {
        if (userParameters->getStructPenalties2() == SECST)
        {
            utilityObject->warning("\n\nWARNING: ignoring secondary structure for a list of sequences\n\n");
        }
        else if (userParameters->getStructPenalties2() == GMASK)
        {
            utilityObject->warning("\n\nWARNING: ignoring gap penalty mask for a list of sequences\n\n");
        }
        userParameters->setUseSS2(false);
    }

    DistMatrix distMat;
    int _numSeqs = alignmentObj.getNumSeqs();
    distMat.ResizeRect(_numSeqs + 1);

    //
    // Convert to similarities!!!!!!!!
    // This is calcualting the similarities of the sequences in the profile part
    //
    for (i = 1; i <= newSeq; i++)
    {
        for (j = i + 1; j <= newSeq; j++)
        {
           dscore = alignmentObj.countid(i,j);
           distMat(i, j) = ((double)100.0 - (double)dscore)/(double)100.0;
           distMat(j, i) = distMat(i, j);
        }
    }
    //distMat.printArray();

    InFileStream _treeFile;  //nige

    useTree = false;
    //
    // Note put this into a separate function!!!!!
   //
    if (_numSeqs >= 2)
    {
        useTree = useExistingGuideTree(Sequences, phylipName, path);
    }

    if (userParameters->getSaveParameters())
    {
        userParameters->createParameterOutput();
    }

    if(userParameters->getResetAlignmentsNew() || userParameters->getResetAlignmentsAll())
    {
        alignmentObj.resetProfile2();
    }
    else
    {
        alignmentObj.fixGaps();
    }

    int _length = 0;
    if (userParameters->getStructPenalties1() == SECST)
    {
        _length = alignmentObj.getSeqLength(1);
        calcGapPenaltyMask(_length, alignmentObj.getSecStructMask1(),
                           alignmentObj.getGapPenaltyMask1());
    }
    if (userParameters->getStructPenalties2() == SECST)
    {
        _length = alignmentObj.getSeqLength(alignmentObj.getProfile1NumSeqs() + 1);
        calcGapPenaltyMask(_length, alignmentObj.getSecStructMask2(),
                           alignmentObj.getGapPenaltyMask2());
    }

    // PROGRESSIVE ALIGNMENT ALGORITHM //
    vector<int> seqWeights(_numSeqs + 1);

    bool success = false;

    if (useTree == false) // create the new tree file, if necessary
    {
        if(userParameters->getDisplayInfo())
        {
            cout << "Start of Pairwise alignments\n";
            cout << "Aligning...\n";
        }

        if(userParameters->getDNAFlag())
        {
            userParameters->setDNAParams();
        }
        else
        {
            userParameters->setProtParams();
        }

        // STEP 1: CALCULATE DISTANCE MATRIX USING PAIRWISE ALIGNMENT //
        PairwiseAlignBase* pairwiseDist;
        if (userParameters->getQuickPairAlign())
        {
            pairwiseDist = new FastPairwiseAlign();
        }
        else
        {
            pairwiseDist = new FullPairwiseAlign();
        }

        // Generate distance matrix!
        pairwiseDist->pairwiseAlign(&alignmentObj, &distMat, 0, _numSeqs, newSeq-2, _numSeqs);
        delete pairwiseDist;

        if(userParameters->getDisplayInfo())
            cout << "\n\n";

        TreeInterface calcSeqWeights;
        calcSeqWeights.getWeightsFromDistMat(&seqWeights, &distMat, &alignmentObj, 1,
                                             _numSeqs, phylipName, &success);
    }
    else
    {
        TreeInterface calcSeqWeights;
        calcSeqWeights.getWeightsFromGuideTree(&alignmentObj, &distMat, phylipName,
        &seqWeights, 1, _numSeqs, &success);
    }

    if(!success)
    {
        return;
    }
    // If it is true, call the function here to get the seqweights etc from it.

    // if users request only the guide tree, return
    if (userParameters->getNewTreeFile())
    {
        return;
    }


    MSA* msaObj = new MSA();
    count = msaObj->seqsAlignToProfile(&alignmentObj, &distMat, &seqWeights, newSeq - 2,
                                       *phylipName);
    delete msaObj;

    userParameters->setUseSS2(saveSS2);

    if (count <= 0)
    {
        return;
    }

    if (userParameters->getMenuFlag())
    {
        cout << "\n\n\n";
    }

    /// STEP 4: OUTPUT THE ALIGNMENT  //
    alignOutput.createAlignmentOutput(&alignmentObj, 1, _numSeqs, output);

    (*phylipName) = "";
}

/** ****************************************************************************
 * The function profileAlign is used to align two profiles                     *
 *******************************************************************************/
void Clustal::profileAlign(string* p1TreeName, string* p2TreeName, ClustalWOutput *output)
{

	cout << "profileAlign called";

    string path;
    //string treeName;
    bool useTree1, useTree2;
    int count, i, j, dscore;
    int _profile1NumSeqs = alignmentObj.getProfile1NumSeqs();
    AlignmentOutput alignOutput;

    if(userParameters->getProfile1Empty() || userParameters->getProfile2Empty())
    {
        utilityObject->error("No sequences in memory. Load sequences first.\n\n");
        return;
    }

    utilityObject->getPath(userParameters->getProfile1Name(), &path);

    if(userParameters->getMenuFlag() || !userParameters->getInteractive())
    {
        if(!alignOutput.openAlignmentOutput(path))
        {
            return;
        }
    }
    else
    {
        // We are using clustalQT.
        // Open all the files.
        if(!alignOutput.QTOpenFilesForOutput(QTFileNames))
        {
            return; // could not open the files.
        }
    }

    if(userParameters->getResetAlignmentsNew() || userParameters->getResetAlignmentsAll())
    {
        alignmentObj.resetProfile1();
        alignmentObj.resetProfile2();
    }
    else
    {
        alignmentObj.fixGaps();
    }

    // Check if there exists a tree for profile1
    useTree1 = false;

    if (_profile1NumSeqs >= 2)
    {
        useTree1 = useExistingGuideTree(Profile1, p1TreeName, path);
    }

    // Check if there exists a tree for profile2
    useTree2 = false;

    utilityObject->getPath(userParameters->getProfile2Name(), &path);

    if (alignmentObj.getNumSeqs() - _profile1NumSeqs >= 2)
    {
        useTree2 = useExistingGuideTree(Profile2, p2TreeName, path);
    }

    if (userParameters->getSaveParameters())
    {
        userParameters->createParameterOutput();
    }
    int _length = 0;
    if (userParameters->getStructPenalties1() == SECST)
    {
        _length = alignmentObj.getSeqLength(1);
        calcGapPenaltyMask(_length, alignmentObj.getSecStructMask1(),
                           alignmentObj.getGapPenaltyMask1());
    }

    if (userParameters->getStructPenalties2() == SECST)
    {
        _length = alignmentObj.getSeqLength(_profile1NumSeqs + 1);
        calcGapPenaltyMask(_length, alignmentObj.getSecStructMask2(),
                           alignmentObj.getGapPenaltyMask2());
    }

    // Declare the distance matrix
    DistMatrix distMat;
    int _numSeqs = alignmentObj.getNumSeqs();
    distMat.ResizeRect(_numSeqs + 1);

    if (useTree1 == false)
    {
        if (_profile1NumSeqs >= 2)
        {
            for (i = 1; i <= _profile1NumSeqs; i++)
            {
                for (j = i + 1; j <= _profile1NumSeqs; j++)
                {
                    dscore = static_cast<int>(alignmentObj.countid(i,j));
                    distMat(i, j) = (100.0 - dscore)/100.0;
                    distMat(j, i) = distMat(i, j);
                }
            }
            utilityObject->getPath(userParameters->getProfile1Name(), &path);

            // We need to get the name of the file first because the message is different!
            if(userParameters->getMenuFlag())
            {
                promptForNewGuideTreeName(Profile1, p1TreeName, path);
            }
            else
            {
                string treeName;
                treeName = path + "dnd";
                p1TreeName = new string(treeName);
            }
        }

        if (useTree2 == false)
        {
            if(_numSeqs - _profile1NumSeqs >= 2)
            {
                for (i = 1 + _profile1NumSeqs; i <= _numSeqs; i++)
                {
                    for (j = i + 1; j <= _numSeqs; j++)
                    {
                        dscore = static_cast<int>(alignmentObj.countid(i,j));
                        distMat(i, j) = (100.0 - dscore) / 100.0;
                        distMat(j, i) = distMat(i, j);
                    }
                }

                utilityObject->getPath(userParameters->getProfile2Name(), &path);

                if(userParameters->getMenuFlag())
                {
                    promptForNewGuideTreeName(Profile2, p2TreeName, path);
                }
                else
                {
                    string treeName;
                    treeName = path + "dnd";
                    p2TreeName = new string(treeName);
                }

            }
        }
    }

    bool success = false;
    vector<int> prof1Weight, prof2Weight;
    prof1Weight.resize(_profile1NumSeqs);
    prof2Weight.resize(_numSeqs);
    TreeInterface tree;
    tree.getWeightsForProfileAlign(&alignmentObj, &distMat, p1TreeName, &prof1Weight,
                                   p2TreeName, &prof2Weight, _numSeqs, _profile1NumSeqs,
                                   useTree1, useTree2, &success);
    if(!success)
    {
        return;
    }

    // do an initial alignment to get the pairwise identities between the two
    // profiles - used to set parameters for the final alignment
    MSA* msaObj = new MSA();

    alignmentObj.resetProfile1();
    alignmentObj.resetProfile2();

    count = msaObj->doProfileAlign(&alignmentObj, &distMat, &prof1Weight, &prof2Weight);
    delete msaObj;

    if (count == 0)
    {
        return;
    }
    if(userParameters->getMenuFlag())
    {
        cout << "\n\n\n";
    }

    alignOutput.createAlignmentOutput(&alignmentObj, 1, _numSeqs, output);

    (*p1TreeName) = "";
    (*p2TreeName) = "";
}

void Clustal::doGuideTreeOnly(string* phylipName)
{
    string path;
    if(userParameters->getEmpty())
    {
        utilityObject->error("No sequences in memory. Load sequences first.\n");
        return;
    }

    userParameters->setStructPenalties1(NONE);
    userParameters->setStructPenalties2(NONE);

    alignmentObj.clearSecStruct1();
    alignmentObj.clearSecStruct2();

    if(userParameters->getResetAlignmentsNew() || userParameters->getResetAlignmentsAll())
    {
        alignmentObj.resetAlign();
    }

    utilityObject->getPath(userParameters->getSeqName(), &path);
    int _numSeqs = alignmentObj.getNumSeqs();

    if (_numSeqs < 2)
    {
        utilityObject->error("Less than 2 sequences in memory. Phylogenetic tree cannot be built.\n");
        return;
    }

    if (userParameters->getSaveParameters())
    {
        userParameters->createParameterOutput();
    }

    if(userParameters->getDisplayInfo())
    {
        cout << "Start of Pairwise alignments\n";
        cout << "Aligning...\n";
    }

    if(userParameters->getDNAFlag())
    {
        userParameters->setDNAParams();
    }
    else
    {
        userParameters->setProtParams();
    }

    ///STEP 1: PAIRWISE ALIGNMENT TO GENERATE DISTANCE MATRIX //

    DistMatrix distMat;
    distMat.ResizeRect(_numSeqs + 1);

    PairwiseAlignBase* pairwiseDist;
    if (userParameters->getQuickPairAlign())
    {
        pairwiseDist = new FastPairwiseAlign();
    }
    else
    {
        pairwiseDist = new FullPairwiseAlign();
    }
    // Generate distance matrix!
    pairwiseDist->pairwiseAlign(&alignmentObj, &distMat, 0, _numSeqs, 0, _numSeqs);
    delete pairwiseDist;

    /* AW DEBUG
       fprintf(stdout, "\nDEBUG: distance matrix following...\n");
       distMat.printArray();
    */

    bool success = false;
    TreeInterface tree;
    tree.generateTreeFromDistMat(&distMat, &alignmentObj, 1, _numSeqs, phylipName, &success);

    if(userParameters->getResetAlignmentsNew() || userParameters->getResetAlignmentsAll())
    {
        alignmentObj.resetAlign();
    }

    (*phylipName) = "";
}


// FIXME this is to  90% identical with align(), please merge
void Clustal::doAlignUseOldTree(string* phylipName, ClustalWOutput *output)
{
	cout << "doAlignUseOldTree called";

    string path;
    int count;
    AlignmentOutput alignOutput;

    if(userParameters->getEmpty())
    {
        utilityObject->error("No sequences in memory. Load sequences first.\n");
        return;
    }
    
    userParameters->setStructPenalties1(NONE);
    userParameters->setStructPenalties2(NONE);

    alignmentObj.clearSecStruct1();
    alignmentObj.clearSecStruct2();

    utilityObject->getPath(userParameters->getSeqName(), &path);

    if(userParameters->getMenuFlag() || !userParameters->getInteractive())
    {
        if(!alignOutput.openAlignmentOutput(path))
        {
            return;
        }
    }
    else
    {
        // We are using clustalQT.
        // Open all the files.
        if(!alignOutput.QTOpenFilesForOutput(QTFileNames))
        {
            return; // could not open the files.
        }
    }

    if(userParameters->getResetAlignmentsNew() || userParameters->getResetAlignmentsAll())
    {
        alignmentObj.resetAlign();
    }

    int _numSeqs = alignmentObj.getNumSeqs();
    DistMatrix distMat;
    distMat.ResizeRect(_numSeqs + 1);
    utilityObject->getPath(userParameters->getSeqName(), &path);

    if (_numSeqs >= 2)
    {

        if(userParameters->getMenuFlag())
        {
            phylipName = new string(path);
            *phylipName = *phylipName + "dnd";

            string message, answer;
            message = "\nEnter a name for the guide tree file [" + *phylipName + "]";
            utilityObject->getStr(message, answer);

            if(!answer.empty())
            {
                phylipName = new string(answer);
            }
        }

        if(userParameters->getMenuFlag() || !userParameters->getInteractive())
        {
            InFileStream _treeFile;  //nige
            _treeFile.open(phylipName->c_str());
            if(!_treeFile.is_open())
            {
                utilityObject->error("Cannot open tree file [%s]\n", phylipName->c_str());
                return;
            }
            _treeFile.close();
        }
    }
    else
    {
        if(userParameters->getDisplayInfo())
        {
            cout << "Start of Pairwise alignments\n";
            cout << "Aligning...\n";
        }
        if(userParameters->getDNAFlag())
        {
            userParameters->setDNAParams();
        }
        else
        {
            userParameters->setProtParams();
        }

        PairwiseAlignBase* pairwiseDist;
        if (userParameters->getQuickPairAlign())
        {
            pairwiseDist = new FastPairwiseAlign();
        }
        else
        {
            pairwiseDist = new FullPairwiseAlign();
        }
        // Generate distance matrix!
        pairwiseDist->pairwiseAlign(&alignmentObj, &distMat, 0, _numSeqs, 0, _numSeqs);
        delete pairwiseDist;
    }

    if (userParameters->getSaveParameters())
    {
        userParameters->createParameterOutput();
    }
    auto_ptr<AlignmentSteps> progSteps;
    vector<int> seqWeights(_numSeqs + 1);
    bool success = false;
    TreeInterface tree;
    progSteps = tree.getWeightsAndStepsFromTree(&alignmentObj, &distMat, phylipName,
                                                &seqWeights, 1, _numSeqs, &success);

    if(!success)
    {
        return;
    }

    MSA* msaObj = new MSA();
    count = msaObj->multiSeqAlign(&alignmentObj, &distMat, &seqWeights, progSteps.get(), 0);
    delete msaObj;

    if (count <= 0)
    {
        return;
    }

    if (userParameters->getMenuFlag())
    {
        cout << "\n\n\n";
    }

    // same as in align()
    if(userParameters->getDoRemoveFirstIteration() == ALIGNMENT)
    {
        //userParameters->setNumIterations(_numSeqs * 2);
        Iteration iterateObj;
        iterateObj.removeFirstIterate(&alignmentObj);
        alignmentObj.calculateMaxLengths(); // Mark Change 20-6-07
        if(userParameters->getDisplayInfo())
            cout << "Finished iteration\n";
    }

    alignOutput.createAlignmentOutput(&alignmentObj, 1, _numSeqs, output);

    phylipName = new string("");
}



void Clustal::getFullHelp()
{
    vector<string> markers;
    Help myhelp;
    bool showtitle = true;

    markers = myhelp.ListSectionMarkers();
    for (unsigned int i=0; i<markers.size(); i++) {
        string m =  markers[i];
        getHelp(m, showtitle);
    }
}

void Clustal::getHelp(char helpPointer, bool printTitle)
{
    string s = "";
    s += helpPointer;
    Clustal::getHelp(s, printTitle);
}



/*
 * Andreas Wilm (UCD): edited to support new help system (separate
 * file now which is compiled in)
 *
 * Author?
 * The clustal help file is called clustalw_help. Should be in the same
 * directory. I have changed it to a C++ style implementation. I am taking
 * out support for VMS until later. We will see if we are still supporting that platform.
 */
void Clustal::getHelp(string helpPointer, bool printTitle)
{
    Help myhelp;
    string helpString;

    helpString =  myhelp.GetSection(helpPointer);

    if (printTitle) {
        helpString = "\n\n>> HELP " +
            helpPointer +
            " <<             " +
            myhelp.GetSectionTitle(helpPointer) +
            "\n\n" +
            helpString;
    }


    if(! userParameters->getMenuFlag())
    {
        cout << helpString;
    }
    else
    {
        string::size_type lastPos = 0;
        string::size_type pos = helpString.find_first_of("\n", lastPos);
        int nlines = 0;

        while (pos != string::npos)
        {
            cout << helpString.substr(lastPos, pos - lastPos);
            nlines++;

            if(nlines >= PAGE_LEN)
            {
                char tempChar;
                cout << "\nPress [RETURN] to continue or  X  to stop ";
                cin.get(tempChar);
                if(toupper(tempChar) == 'X')
                {
                    return;
                }
                else
                {
                    nlines = 0;
                }
            }

            lastPos = pos; //helpString.find_first_not_of("\n", pos);
            pos = helpString.find_first_of("\n", lastPos+1);
            //cerr << "DEBUG: pos=" << pos << " lastPos=" << lastPos << "/" << helpString.length() << "\n";
        }
    }
}



/**
 * The wrap functions will be used with interface classes to perform the tasks
 * that were previously done there.
 */
int Clustal::sequenceInput(bool append, string *offendingSeq)
{
    int code;
    // If we are not appending, we need to clear the Alignment object.
    if(!append)
    {
        alignmentObj.clearAlignment();
    }

    FileReader readSeqFile;
    code = readSeqFile.seqInput(&alignmentObj, append, offendingSeq);
    return code;
}

/**
 * profile1Input is a wrapper function for the profileInput. This is because the
 * other classes dont have access to FileReader
 */
int Clustal::profile1Input(string profile1Name)
{
    int code;
    // I need to clear out the Alignment object.
    alignmentObj.clearAlignment();
    userParameters->setProfileNum(1);

    userParameters->setSeqName(profile1Name);
    userParameters->setProfile1Name(profile1Name);

    FileReader readProfileFile;
    code = readProfileFile.profileInput(&alignmentObj);

    // If we are using the commandline check if there are seqs!
    if(!userParameters->getInteractive())
    {
        // AW: FIXME code should be handled higher up the stack and check all codes
        // also, shouldnt we use  utilityObject->error()?
        if(code != OK)
        {
            if (code==NOSEQUENCESINFILE)
                cerr << "ERROR: There are no sequences in profile2 file." << std::endl;
            else if (code==ALLNAMESNOTDIFFERENT)
                cerr << "ERROR: Not all sequence names are different" << std::endl;
            else
                cerr << "ERROR: Unhandled error code (" << code << ") returned from profileInput.\n";
            // AW: should we really exit here? What if called from clustalx?
            throw 2;
        }
    }

    return code;
}

/**
 * profile2Input is a wrapper function for the profileInput. This is because the
 * other classes dont have access to FileReader
 */
int Clustal::profile2Input(string profile2Name)
{
    int code;

    if(userParameters->getProfileNum() == 2)
    {
        // Remove the sequences from the previous one.
        int numSeqsProfile1 = alignmentObj.getProfile1NumSeqs();
        alignmentObj.resizeSeqArray(numSeqsProfile1 + 1);
        // Clear anything from profile2 in alignment.
        alignmentObj.clearSecStruct2();
    }
    userParameters->setProfileNum(2);

    userParameters->setSeqName(profile2Name);
    userParameters->setProfile2Name(profile2Name);


    FileReader readProfileFile;
    cout << "before profileInput\n";
    code = readProfileFile.profileInput(&alignmentObj);
    cout << "after profileInput\n";

    if(!userParameters->getInteractive())
    {
        // AW: FIXME code should be handled higher up the stack and check all codes
        // also, shouldnt we use  utilityObject->error()?
        if(code != OK) {
            if (code==NOSEQUENCESINFILE)
                cerr << "ERROR: There are no sequences in profile2 file." << std::endl;
            else if (code==ALLNAMESNOTDIFFERENT)
                cerr << "ERROR: Not all sequence names are different" << std::endl;
            else
                cerr << "ERROR: Unhandled error code (" << code << ") returned from profileInput.\n";
            // AW: should we really exit here? What if called from clustalx?
            // DD: fixed
            if(!userParameters->getGui())
                throw 2;
        }
    }
    return code;
}


/**
 * The function commandline_readseq is called by the command line to
 * read in the sequences. It also prints out the names. This was
 * previously done by the command line parser.
 */
int Clustal::commandLineReadSeq(int firstSeq, ClustalWInput *input)
{
    // Clear the alignment, although obviously there shouldnt be anything in it.
    alignmentObj.clearAlignment();
    userParameters->setProfileNum(0);
    int code = 0;
    string offendingSeq;
    FileReader readInputFile;
    string line = userParameters->getSeqName();

    if (strcmp(line.c_str(), "internalRsequence") == 0) {
    	code = readInputFile.readCharacterSeqs(&alignmentObj, firstSeq, &offendingSeq, input);
    } else {
    	code = readInputFile.readSeqs(&alignmentObj, firstSeq, &offendingSeq);
    }

    if(code != OK)
    {
        if(code == CANNOTOPENFILE)
        {
            utilityObject->error("Cannot open input file. No alignment!\n");
        }
        else if(code == NOSEQUENCESINFILE)
        {
            utilityObject->error("No sequences in file. No alignment!\n");
        }
        else if(code == ALLNAMESNOTDIFFERENT)
        {
            utilityObject->error("Multiple sequences found with same name (found %s at least twice)!", offendingSeq.c_str());
        }
        else if(code == EMPTYSEQUENCE)
        {
            utilityObject->error("Empty sequences found: %s\n", offendingSeq.c_str());
        }
        else if(code == SEQUENCETOOBIG)
        {
            utilityObject->error("Sequence(s) too big: %s\n", offendingSeq.c_str());
        }
        else if(code == BADFORMAT)
        {
            utilityObject->error("Sequences are badly formatted!\n");
        }
        else
        {
            utilityObject->error("\nThere was a problem reading in the file. No alignment!\n");
        }
        throw -1;
    }

    alignmentObj.printSequencesAddedInfo();

    userParameters->setEmpty(false);
    return code;
}

/*
 * The function outputNow is used to output the alignment. It can be called by the
 * menu or command line parser.
 */
void Clustal::outputNow(ClustalWOutput *output)
{
	cout << "outputNow called";
    if(alignmentObj.getNumSeqs() > 0)
    {
        string path = "";
        if(!userParameters->getMenuFlag())
        {
            string _seqName = userParameters->getSeqName();
            utilityObject->getPath(_seqName, &path);
        }
        AlignmentOutput alignmentOutput;
        alignmentOutput.openAlignmentOutput(path);

        alignmentOutput.createAlignmentOutput(&alignmentObj, 1, alignmentObj.getNumSeqs(), output);
    }
    else
    {
        utilityObject->error("No sequences have been loaded\n");
    }
    return;
}

void Clustal::phylogeneticTree(string* phylipName, string* clustalName, string* distName,
                               string* nexusName, string pimName)
{
    TreeNames treeNames;
    treeNames.clustalName = *clustalName;
    treeNames.distName = *distName;
    treeNames.nexusName = *nexusName;
    treeNames.phylipName = *phylipName;
    treeNames.pimName = pimName;
    TreeInterface tree;
    tree.treeFromAlignment(&treeNames, &alignmentObj);
}

void Clustal::bootstrapTree(string* phylipName, string* clustalName, string* nexusName)
{
    TreeNames treeNames;
    treeNames.clustalName = *clustalName;
    treeNames.nexusName = *nexusName;
    treeNames.phylipName = *phylipName;
    TreeInterface tree;
    tree.bootstrapTree(&treeNames, &alignmentObj);
}

void Clustal::initInterface()
{
    userParameters->setEmpty(true);
    userParameters->setProfile1Empty(true);
    userParameters->setProfile2Empty(true);
}

/**
 * The function calcGapPenaltyMask is used to calculate the gapMask from the secondary
 * structure information.
 */
void Clustal::calcGapPenaltyMask(int prfLength, vector<char>* mask, vector<char>* gapMask)
{
    int i,j;

    vector<char> structMask;
    structMask.resize(prfLength + 1);

    int _helixEndPlus = userParameters->getHelixEndPlus();
    int _helixEndMinus = userParameters->getHelixEndMinus();
    int _strandEndPlus = userParameters->getStrandEndPlus();
    int _strandEndMinus = userParameters->getStrandEndMinus();

    i = 0;
    while (i < prfLength)
    {
        if (tolower((*mask)[i]) == 'a' || (*mask)[i] == '$')
        {
            for (j = -_helixEndPlus; j < 0; j++)
            {
                if(i + j < (int)structMask.size())
                {
                    if ((i + j >= 0) && (tolower(structMask[i + j]) != 'a') &&
                        (tolower(structMask[i + j]) != 'b'))
                    {
                        structMask[i + j] = 'a';
                    }
                }
            }
            for (j = 0; j < _helixEndMinus; j++)
            {
                if (i + j >= prfLength || (tolower((*mask)[i + j]) != 'a'
                                    && (*mask)[i + j] != '$'))
                {
                    break;
                }
                structMask[i + j] = 'a';
            }
            i += j;
            while (tolower((*mask)[i]) == 'a' || (*mask)[i] == '$')
            {
                if (i >= prfLength)
                {
                    break;
                }
                if ((*mask)[i] == '$')
                {
                    structMask[i] = 'A';
                    i++;
                    break;
                }
                else
                {
                    structMask[i] = (*mask)[i];
                }
                i++;
            }
            for (j = 0; j < _helixEndMinus; j++)
            {
                if ((i - j - 1 >= 0) && (tolower((*mask)[i - j - 1]) == 'a'
                                         || (*mask)[i - j - 1] == '$'))
                {
                    structMask[i - j - 1] = 'a';
                }
            }
            for (j = 0; j < _helixEndPlus; j++)
            {
                if (i + j >= prfLength)
                {
                    break;
                }
                structMask[i+j] = 'a';
            }
        }
        else if (tolower((*mask)[i]) == 'b' || (*mask)[i] == '%')
        {
            for (j = -_strandEndPlus; j < 0; j++)
            {
                if ((i + j >= 0) && (tolower(structMask[i + j]) != 'a')
                                 && (tolower(structMask[i + j]) != 'b'))
                {
                    structMask[i + j] = 'b';
                }
            }
            for (j = 0; j < _strandEndPlus; j++)
            {
                if (i + j >= prfLength || (tolower((*mask)[i + j]) != 'b'
                                          && (*mask)[i + j] != '%'))
                {
                    break;
                }
                structMask[i+j] = 'b';
            }
            i += j;
            while (tolower((*mask)[i]) == 'b' || (*mask)[i] == '%')
            {
                if (i >= prfLength)
                {
                    break;
                }
                if ((*mask)[i] == '%')
                {
                    structMask[i] = 'B';
                    i++;
                    break;
                }
                else
                {
                    structMask[i] = (*mask)[i];
                }
                i++;
            }
            for (j = 0; j < _strandEndMinus; j++)
            {
                if ((i-j-1>=0) && (tolower((*mask)[i-j-1]) == 'b' || (*mask)[i-j-1] == '%'))
                {
                    structMask[i - j - 1] = 'b';
                }
            }
            for (j = 0; j < _strandEndPlus; j++)
            {
                if (i + j >= prfLength)
                {
                    break;
                }
                structMask[i+j] = 'b';
            }
        }
        else
        {
            i++;
        }
    }

    for(i = 0; i < prfLength;i++)
    {
        switch (structMask[i])
        {
            case 'A':
                (*gapMask)[i] = userParameters->getHelixPenalty() + '0';
                break;
            case 'a':
                (*gapMask)[i] = userParameters->getHelixEndPenalty() + '0';
                break;
            case 'B':
                (*gapMask)[i] = userParameters->getStrandPenalty() +'0';
                break;
            case 'b':
                (*gapMask)[i] = userParameters->getStrandEndPenalty() + '0';
                break;
            default:
                (*gapMask)[i] = userParameters->getLoopPenalty() + '0';
                break;
        }
    }
}


/**
 * The function QTcalcLowScoreSegments is used to calculate the residues in the sequences
 * that score badly.
 * @param firstSeq first seq in the alignment or profile
 * @param nSeqs the number of sequences
 * @param nCols the length of the longest seq
 * @param seqWeight a vector to hold the sequence weights
 * @param lowScoreRes
 * @param seqWeightCalculated
 */
void Clustal::QTcalcLowScoreSegments(LowScoreSegParams* params)
{
    int i, j;
    float sum, prevSum;
    float gscale;
    vector<int> weight;
    int sweight;
    vector<int> gaps;
    int matrix[NUMRES][NUMRES];
    vector<float> fsum;
    vector<float> bsum;
    vector<float> pscore;
    int _maxAA = userParameters->getMaxAA();

    // STEP 1: Calculate the sequence weights
    QTcalcWeightsForLowScoreSeg(params);

    subMatrix->getQTMatrixForLowScoreSeg(matrix);

    const SeqArray* seqArray = alignmentObj.getSeqArray();

    gaps.resize(params->nCols + 1);

    for (j = 1; j <= params->nCols; j++)
    {
        gaps[j - 1] = 0;
        for(i = params->firstSeq + 1; i < params->firstSeq + params->nSeqs; i++)
        {
            if (j < alignmentObj.getSeqLength(i))
            {
                if (((*seqArray)[i][j] < 0) || ((*seqArray)[i][j] > _maxAA))
                {
                    gaps[j-1]++;
                }
            }
        }
    }

    // STEP 2: Calculate the profile
    LowScoreSegProfile lowScoreProfile(params->nCols, params->firstSeq,
                                       params->firstSeq + params->nSeqs);

    weight.resize(params->firstSeq + params->nSeqs + 1);

    for(i = params->firstSeq;i < params->firstSeq + params->nSeqs;i++)
    {
        weight[i] = (*(params->seqWeight))[i - params->firstSeq];
    }

    lowScoreProfile.calcLowScoreSegProfile(seqArray, matrix, &weight);

    const SeqArray* profile = lowScoreProfile.getProfilePtr();

    sweight = 0;
    for(i = params->firstSeq; i < params->firstSeq + params->nSeqs; i++)
    {
        sweight += weight[i];
    }

    //Now, use the profile scores to mark segments of each sequence which score badly.

    fsum.resize(params->nCols + 2);
    bsum.resize(params->nCols + 2);
    pscore.resize(params->nCols + 2);

    for(i = params->firstSeq + 1; i < params->firstSeq + params->nSeqs + 1; i++)
    {
    // In a forward phase, sum the profile scores. Mark negative sums as exceptions.
    //If the sum is positive, then it gets reset to 0.
        sum = 0.0;
        for(j = 1; j <= alignmentObj.getSeqLength(i); j++)
        {
            gscale = (float)(params->nSeqs - gaps[j - 1]) / (float)params->nSeqs;
            if((*seqArray)[i][j] < 0 || (*seqArray)[i][j] >= _maxAA)
            {
                pscore[j - 1] = 0.0;
                sum = 0.0;
            }
            else
            {
                pscore[j-1]=((*profile)[j][(*seqArray)[i][j]]- weight[i - 1] *
                             matrix[(*seqArray)[i][j]][(*seqArray)[i][j]]) * gscale / sweight;
            }
            sum += pscore[j - 1];
            if(sum > 0.0)
            {
                sum = 0.0;
            }
            fsum[j - 1] = sum;
        }
// trim off any positive scoring residues from the end of the segments
        prevSum = 0;
        for(j = alignmentObj.getSeqLength(i) - 1; j >= 0; j--)
        {
            if(prevSum >= 0.0 && fsum[j] < 0.0 && pscore[j] >= 0.0)
            {
                fsum[j] = 0.0;
            }
            prevSum = fsum[j];
        }

// Now, in a backward phase, do the same summing process.
        sum = 0.0;
        for(j = alignmentObj.getSeqLength(i); j >= 1; j--)
        {
            if((*seqArray)[i][j] < 0 || (*seqArray)[i][j] >= _maxAA)
            {
                sum = 0;
            }
            else
            {
                sum += pscore[j - 1];
            }
            if(sum > 0.0)
            {
                sum = 0.0;
            }
            bsum[j - 1] = sum;
        }
// trim off any positive scoring residues from the start of the segments
        prevSum = 0;
        for(j = 0; j < alignmentObj.getSeqLength(i); j++)
        {
            if(prevSum >= 0.0 && bsum[j] < 0.0 && pscore[j] >= 0.0)
            {
                bsum[j] = 0.0;
            }
            prevSum = bsum[j];
        }
//Mark residues as exceptions if they score negative in the forward AND backward directions.
        for(j = 1; j <= alignmentObj.getSeqLength(i); j++)
        {
            if(fsum[j - 1] < 0.0 && bsum[j - 1] < 0.0)
            {
                if((*seqArray)[i][j] >= 0 && (*seqArray)[i][j] < _maxAA)
                {
                    (*(params->lowScoreRes))[i - params->firstSeq - 1][j - 1] = -1;
                }
            }
        }
    }

// Finally, apply the length cutoff to the segments - removing segments shorter
//than the cutoff

    QTremoveShortSegments(params);

}

void Clustal::QTremoveShortSegments(LowScoreSegParams* params)
{
    int i,j,k,start;
    //panel_data data;

    //GetPanelExtra(p,&data);
    if(params->nSeqs <= 0)
        return;

// Reset all the exceptions - a value of 1 indicates an exception that
// will be displayed. A value of -1 is used to remember exceptions that
// are temporarily hidden in the display

    for(i = 0; i < params->nSeqs; i++)
    {
        for(j = 0; j < params->nCols; j++)
        {
            if((*(params->lowScoreRes))[i][j] == -1)
            {
                (*(params->lowScoreRes))[i][j] = 1;
            }
        }
    }

    for(i = 0; i < params->nSeqs; i++)
    {
        start = -1;
        for(j = 0; j <= params->nCols; j++)
        {
            if(start == -1)
            {
                if((*(params->lowScoreRes))[i][j] == 1)
                    start = j;
            }
            else
            {
                if(j == params->nCols || (*(params->lowScoreRes))[i][j] == 0)
                {
                    if(j - start < userParameters->getQTminLenLowScoreSegment())
                    {
                        for(k = start; k < j; k++)
                        {
                            (*(params->lowScoreRes))[i][k] = -1;
                        }
                    }
                    start = -1;
                }
            }
        }
    }
}

/**
 * Change: Mark 22-5-07, changed the distmatrix to be the size of alignObject.numSeqs
 */
void Clustal::QTcalcWeightsForLowScoreSeg(LowScoreSegParams* params)
{
    int i, j;
    vector<int> weight;
    float dscore;
    DistMatrix distMat(alignmentObj.getNumSeqs() + 1); // Mark: changed size
    // Aw potential trouble here: what if we don't have write
    // permission to current directory?
#ifdef UNIX
    char treeName[FILENAMELEN]=".score.ph";
#else
    char treeName[FILENAMELEN]="tmp.ph";
#endif

    if(params->nSeqs <= 0)
    {
        return;
    }

// if sequence weights have been calculated before - don't bother
//doing it again (it takes too long). data.seqweight is set to NULL when
// new sequences are loaded. /
    if(params->seqWeightCalculated == true)
    {
        return;
    }

    utilityObject->info("Calculating sequence weights...");

    // count pairwise percent identities to make a phylogenetic tree //
    if(params->nSeqs >= 2)
    {       //i=firstSeq + 1; i <= firstSeq + nSeqs
        for (i = params->firstSeq + 1; i <= params->firstSeq + params->nSeqs; i++)
        {   //j=i + 1; i <= firstSeq + nSeqs
            for (j = i + 1; j <= params->firstSeq + params->nSeqs; j++) // Mark 22-5-07
            {
                dscore = alignmentObj.countid(i, j); // Mark 22-5-07
                distMat(i, j) = (100.0 - dscore) / 100.0;
                distMat(j, i) = distMat(i, j);
            }
        }
        string name = string(treeName);
        bool success = false;
        weight.resize(params->firstSeq + params->nSeqs + 1);
        TreeInterface tree;
        tree.getWeightsForQtLowScore(&weight, &distMat, &alignmentObj, params->firstSeq + 1,
                                   params->nSeqs, &name, &success); // Mark change 10-5-07
        if(!success)
        {
            return;
        }

        for(i = params->firstSeq;i < params->firstSeq + params->nSeqs;i++)
        {
            (*(params->seqWeight))[i - params->firstSeq] = weight[i];
        }

        utilityObject->info("Done.");
    }
}

void Clustal::QTSetFileNamesForOutput(AlignmentFileNames fileNames)
{
    QTFileNames = fileNames;
}

bool Clustal::QTRealignSelectedRange(AlignmentFileNames fileNames, int beginPos, int endPos, bool realignEndGapPen, ClustalWOutput *output)
{
	cout << "QTRealignSelectedRange called";
    bool alignEndGapPen = userParameters->getEndGapPenalties();

    Alignment saveOldAlign = alignmentObj; // Take a copy of it. Note provided copy
                                           // constructor is ok.
    bool ok;
    ok = alignmentObj.removeAllOutsideRange(beginPos, endPos);

    if(!ok)
    {
        alignmentObj = saveOldAlign;
        return false;
    }
    // Temporarily set the alignment output to be input
    int saveOutOrder = userParameters->getOutputOrder();
    userParameters->setOutputOrder(INPUT);

    //set end gap penalties to be realignEndGapPenalties
    userParameters->setEndGapPenalties(realignEndGapPen);
    //do the alignment
    if(alignmentObj.getNumSeqs() <= 0)
    {
        alignmentObj = saveOldAlign;
        return false;
    }
    QTSetFileNamesForOutput(fileNames);
    string phylipName = fileNames.treeFile;
    align(&phylipName, output, false);

    userParameters->setOutputOrder(saveOutOrder);

    // reset the end gap penalties
    userParameters->setEndGapPenalties(alignEndGapPen);

    // remove postions that only contain gaps
    int nSeqs = alignmentObj.getNumSeqs();
    alignmentObj.removeAllGapOnlyColumns(1, nSeqs, 0);

    // save it to a temporary area.
    SeqArray realignedArea = *(alignmentObj.getSeqArray());

    // Paste it back into the original alignment.
    alignmentObj = saveOldAlign;
    bool result;
    result = alignmentObj.updateRealignedRange(realignedArea, beginPos, endPos);
    if(!result)
    {
        utilityObject->error("something went wrong while updating the realigned range\n");
    }

    // output the alignments
    AlignmentOutput alignOutput;
    if(!alignOutput.QTOpenFilesForOutput(QTFileNames))
    {
        return false; // could not open the files.
    }

    alignOutput.createAlignmentOutput(&alignmentObj, 1, nSeqs, output);

    return true;
}

void Clustal::test()
{
    cout << "RUNNING TEST\n";
    ClustalWOutput *output = new ClustalWOutput();
    AlignmentOutput alignOutput;
    string path;
    utilityObject->getPath(userParameters->getSeqName(), &path);

    if(!alignOutput.openAlignmentOutput(path))
    {
        cerr << "could not open the file\n";
        return;
    }

    vector<int> selected;
    int nSeqs = alignmentObj.getNumSeqs();
    selected.resize(nSeqs + 1, 0);
    selected[9] = 1;
    selected[10] = 1;
    //selected[1] = 1;
    alignmentObj.removeGapOnlyColsFromSelectedSeqs(&selected);
    alignOutput.createAlignmentOutput(&alignmentObj, 1, nSeqs, output);
}

/**
 * This function is used to ask the user if they would like to use an existing guide tree.
 * Note: It expects that phylipName actually points to a string that has been allocated.
 */
bool Clustal::useExistingGuideTree(int type, string* phylipName, const string& path)
{
    bool useTree = false;
    string treeName;
    InFileStream _treeFile;  //nige
    bool paramUseTree;

    string* ptrToMsg;
    if(type == Sequences)
    {
        ptrToMsg = &sequencesMsg;
        paramUseTree = userParameters->getUseTreeFile();
    }
    else if(type == Profile1)
    {
        ptrToMsg = &profile1Msg;
        paramUseTree = userParameters->getUseTree1File();;
    }
    else if(type == Profile2)
    {
        ptrToMsg = &profile2Msg;
        paramUseTree = userParameters->getUseTree2File();;
    }
    else
    {
        ptrToMsg = &sequencesMsg;
        paramUseTree = userParameters->getUseTreeFile();
    }

    if (checkTree && userParameters->getMenuFlag())
    {
        treeName = path + "dnd";

        _treeFile.open(treeName.c_str());
        _treeFile.seekg(0, std::ios::beg);

        if(_treeFile.is_open())
        {
            string message = *ptrToMsg + treeName + "  (y/n) ? [y]";
            string answer;
            utilityObject->getStr(message, answer);

            if(answer[0] != 'n' && answer[0] != 'N')
            {
                if(!phylipName)
                {
                    phylipName = new string(treeName);
                }
                else
                {
                    *phylipName = treeName;
                }
                useTree = true;
            }
            _treeFile.close();
        }
    }
    else if (!userParameters->getMenuFlag() && paramUseTree)
    {
        useTree = true;
    }

    return useTree;
}

void Clustal::promptForNewGuideTreeName(int type, string* treeName, const string& path)
{
    string* ptrToMsg;
    if(type == Profile1)
    {
        ptrToMsg = &newProfile1TreePrompt;
    }
    else if(type == Profile2)
    {
        ptrToMsg = &newProfile2TreePrompt;
    }
    else
    {
        ptrToMsg = &newProfile1TreePrompt;
    }

    if(!treeName)
    {
        treeName = new string("");
    }

    while(treeName->empty())
    {
        string message = *ptrToMsg + path + "dnd]";
        string answer;
        utilityObject->getStr(message, answer);
        if(answer.empty())
        {
            answer = path + "dnd";
            *treeName = answer;
        }
        else
        {
            *treeName = answer;
        }
    }
}

}
