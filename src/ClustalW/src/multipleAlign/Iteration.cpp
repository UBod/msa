/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.
 *
 * Changes:
 * Mark 30-5-2007: Changed iterationOnTreeNode function as it was adding in extra gaps
 * at the end of an alignment.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "Iteration.h"
#include "../alignment/ObjectiveScore.h"
#include "../general/utils.h"
#include "../general/userparams.h"
#include "../tree/TreeInterface.h"
#include "../clustalw_version.h"
#include "MSA.h"

namespace clustalw
{

bool Iteration::iterationOnTreeNode(int numSeqsProf1, int numSeqsProf2, int& prfLength1,
                                    int& prfLength2, SeqArray* seqArray)
{
    Alignment alignmentToIterate;
    int numSeqsInProfiles = numSeqsProf1 + numSeqsProf2;

    if(numSeqsInProfiles <= 2)
    {
        return false;
    }
        
    SeqArray profileSeqs; 
    profileSeqs.resize(numSeqsInProfiles + 1);

    // Copy the SeqArray!
    for (int j = 0; ((j < numSeqsProf1 + numSeqsProf2) && 
		     (j < (int)seqArray->size())); j++)
    {
        profileSeqs[j + 1].clear();
        profileSeqs[j + 1].resize(prfLength1 + 1);
        for (int i = 0; i < prfLength1 && i < (int)(*seqArray)[j].size(); i++)
        {
            profileSeqs[j + 1][i + 1] = (*seqArray)[j][i];
        }
    }
        
    alignmentToIterate.addSequences(&profileSeqs);
    //userParameters->setNumIterations(numSeqsInProfiles * 2);
        
    bool changed = false;
    changed = removeFirstIterate(&alignmentToIterate);
        
    if(changed)
    {          
        SeqArray* iteratedSeqs = alignmentToIterate.getSeqArrayForRealloc();
        string aaCodes = userParameters->getAminoAcidCodes();
        
        int newPrf1Length = 0, newPrf2Length = 0;    
        
        for (int j = 0; j < numSeqsProf1 + numSeqsProf2; j++)
        {            
            if(j < numSeqsProf1)
            { 
                if(alignmentToIterate.getSeqLength(j + 1) > newPrf1Length)
                {
                    newPrf1Length = alignmentToIterate.getSeqLength(j + 1);
                }
            }
            else if(j < numSeqsProf1 + numSeqsProf2)
            {
                if(alignmentToIterate.getSeqLength(j + 1) > newPrf2Length)
                {
                    newPrf2Length = alignmentToIterate.getSeqLength(j + 1);
                }
            }              
        }
        
        prfLength1 = newPrf1Length; // mark 30-5-2007
        prfLength2 = newPrf2Length; // mark 30-5-2007

        for (int j = 0; j < numSeqsProf1 + numSeqsProf2; j++)
        {        
            // I need to recalculate the prfLength1 and prfLength2
            (*seqArray)[j].clear();
            (*seqArray)[j].assign((*iteratedSeqs)[j + 1].begin() + 1, 
                                  (*iteratedSeqs)[j + 1].end()); 
            (*seqArray)[j].resize(prfLength1 + extraEndElemNum, 31);
            (*seqArray)[j][prfLength1] = ENDALN;
        }
    }   

    return true;
}

void Iteration::printSeqArray(SeqArray* arrayToPrint)
{
    cout << "HERE IS THE SEQARRAY\n";
    // I need to use iterators for everything here.
    SeqArray::iterator mainBeginIt = arrayToPrint->begin();
    SeqArray::iterator mainEndIt = arrayToPrint->end();
    vector<int>::iterator begin, end;
    string aaCodes = userParameters->getAminoAcidCodes();
    
    for(; mainBeginIt != mainEndIt; mainBeginIt++)
    {
        if(mainBeginIt->size() > 0)
        {
            begin = mainBeginIt->begin() + 1;
            end = mainBeginIt->end();
            for(; begin != end; begin++)
            {
                if(*begin < (int)aaCodes.size())
                {
                    cout << aaCodes[*begin];
                }
                else
                {
                    cout << "-";
                }
            }
            cout << "\n";
        }
    }
    cout << "\n\n";
}

/**
 * The function removeFirstIterate is used to perform the remove first iteration
 * strategy on the finished alignment. It optimises the score give in alignScore.
 *   For iter = 1 to numIterations
 *      For seq i = 1 to numSeqs
 *         remove seq i
 *         if either of the profiles has all gaps, remove this column.
 *         realign using profileAlign
 *         if its better, keep it. If its not better, dont keep it. 
 * @param alnPtr The alignment object.
 * @return true if it has been successful, false if it has not been successful.
 */
bool Iteration::removeFirstIterate(Alignment* alnPtr)
{   
    if(!alnPtr)
    {
        return false;
    }
    
    string p1TreeName;
    p1TreeName = "";
    string p2TreeName;
    int nSeqs = alnPtr->getNumSeqs();
    
    if(nSeqs <= 2)
    {
        return false;
    }
    DistMatrix distMat;    
    distMat.ResizeRect(nSeqs + 1);
    
    ObjectiveScore scoreObj;
    int iterate = userParameters->getDoRemoveFirstIteration();
    userParameters->setDoRemoveFirstIteration(NONE);
           
    double firstScore = scoreObj.getScore(alnPtr);
    //cout << "firstScore = " << firstScore << "\n";
    double score = 0;
    double bestScore = firstScore;
    double dscore;
    int count;
    bool scoreImproved = false;
    bool scoreImprovedAnyIteration = false;
    int prof1NumSeqs = 1;
    
    // This will be used for removing gaps!!!
    vector<int> profile1;
    vector<int> profile2;
    profile1.resize(nSeqs + 1, 0);
    profile1[1] = 1;
    profile2.resize(nSeqs + 1, 1);
    profile2[0] = 0;
    profile2[1] = 0;
    vector<int> prof1Weight, prof2Weight;
    int iterations = userParameters->getNumIterations();
    //cout << "Max num iterations = " << iterations << "\n";
    Alignment bestAlignSoFar;
    TreeInterface tree;
    // One iteration consists of removing each of the sequences, reseting all the gap
    // only columns. If the score is better, the new alignment is kept.
    for(int n = 1; n <= iterations; n++)
    {
        scoreImproved = false;
        cout << "ITERATION " << n << " OF " << iterations << "\n";
        for(int i = 1; i <= nSeqs; i++)
        {
            vector<Sequence> seqVector;
            Alignment iterateAlign = *alnPtr;

            iterateAlign.setProfile1NumSeqs(1);

            // We remove the sequence i from the profile, and paste into the first position
            // This is to make it easy to do the profile alignment.
            vector<int> selected;
            selected.resize(nSeqs + 1, 0);
            selected[i] = 1;
            seqVector = iterateAlign.cutSelectedSequencesFromAlignment(&selected);
            iterateAlign.pasteSequencesIntoPosition(&seqVector, 0);

            // Remove any gap only columns
            iterateAlign.removeGapOnlyColsFromSelectedSeqs(&profile1);
            iterateAlign.removeGapOnlyColsFromSelectedSeqs(&profile2);

            // Calculate a simple distance matrix.
            if(nSeqs - 1 >= 2) 
            {
                for (int i = 1; i <= nSeqs; i++) 
                {
                    for (int j = i + 1; j <= nSeqs; j++) 
                    {
                        dscore = iterateAlign.countid(i, j);
                        distMat(i, j) = (100.0 - dscore)/100.0;
                    }
                }

                /* temporary tree file
                 *  
                 * can't use the safer mkstemp function here, because
                 * we just pass down the filename :(
                 */
                char buffer[L_tmpnam];
                tmpnam (buffer);
                p2TreeName = buffer + string(".dnd");
                // should test here if file is writable
            }
            bool success = false;
            prof1Weight.clear();
            prof1Weight.resize(prof1NumSeqs);
            prof2Weight.clear();
            prof2Weight.resize(nSeqs);

            tree.getWeightsForProfileAlign(&iterateAlign, &distMat, &p1TreeName, &prof1Weight,
                                           &p2TreeName, &prof2Weight, nSeqs, prof1NumSeqs,
                                           false, false, &success);
            remove(p2TreeName.c_str());
            if(!success)
            {
                /* returning false only means alignment hasn't
                 * changed, but here getWeightsForProfileAlign failed,
                 * most likely because p2TreeName couldn't be read. an
                 * error will be printed to console.  clustalw should
                 * then exit, FIXME: clustalx users have to sit
                 * through all error messages until someone
                 * implements a way to return an exit code and react
                 * appropriately
                 */
                // does anyone know how to use
                // (userParameters->getMenuFlag() ||
                // !userParameters->getInteractive() instead?
                char buf[1024];
                utilityObject->myname(buf);
                if (strcasecmp(buf, "clustalw")==0) {
                    throw EXIT_FAILURE;
                } else {
                    // the next two lines were here before the exit
                    // was added. keeping it for clustalx although it
                    // doesnt seem to make any sens
                    userParameters->setDoRemoveFirstIteration(iterate);
                    return false;
                }
            }
                        
            MSA* msaObj = new MSA();

            iterateAlign.resetProfile1();
            iterateAlign.resetProfile2();
            // Do the profile alignment.
            count = msaObj->doProfileAlign(&iterateAlign, &distMat, 
                                            &prof1Weight, &prof2Weight);   
            delete msaObj;
            // Check if its better
            score = scoreObj.getScore(&iterateAlign);
            iterateAlign.setProfile1NumSeqs(0);
            
            if(score < bestScore) // Might be a problem with this.
            {
                //cout << "**********************************************\n";
                //cout << "***** Better score found using iteration *****\n";
                //cout << "**********************************************\n";
                bestScore = score;
                bestAlignSoFar = iterateAlign;
                scoreImproved = true;
                scoreImprovedAnyIteration = true;
            }
            distMat.clearArray();
            distMat.ResizeRect(nSeqs + 1);   
        }
        if(scoreImproved == false)
        {
            cout << "Score was not improved in last iteration. Exiting...\n";
            break;
        }
    }
    
    //
    // NOTE if we have improved it, then we need to update the sequences in alnPtr
    // 1) get the unique id of seq i
    // 2) get the sequence from new object using id
    // 3) update the sequence in alnPtr
    //
    if(scoreImprovedAnyIteration) // If we need to update the alnPtr object.
    {
        cout << "Iteration improved Align score: " << bestScore << "\n";
        int seqId;
        const vector<int>* improvedSeq;
        for(int i = 1; i <= nSeqs; i++) // For each seq in alnPtr
        {
            seqId = alnPtr->getUniqueId(i);
            improvedSeq = bestAlignSoFar.getSequenceFromUniqueId(seqId);
            alnPtr->updateSequence(i, improvedSeq);
        }
    }
    cout << "FINAL score: " << bestScore << "\n";
    userParameters->setDoRemoveFirstIteration(iterate);
    return true; // It was successful.
}

}

