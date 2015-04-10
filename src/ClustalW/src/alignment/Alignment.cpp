/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * Changes: 
 *
 * 16-02-07,Nigel Brown(EMBL): Added friend NameIterator to allow a caller to
 * process the name vector.
 *
 * 05-03-07,Nigel Brown(EMBL): modified searchForString() to skip over gaps in
 * sequence while matching.
 *
 * 26-03-07,Nigel Brown(EMBL): suppressed error message box for name clashes;
 * added testUniqueNames() predicate, which compares new sequence names with
 * those in the alignment vector BEFORE appending them.
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <exception>
#include <cmath>
#include <sstream>
#include "Alignment.h"

using namespace std;

namespace clustalw
{

Alignment::Alignment()
 : gapPos1(userParameters->getGapPos1()),
   gapPos2(userParameters->getGapPos2())
{
    maxNames = 0;
    numSeqs = 0;
    maxAlignmentLength = 0;
    lengthLongestSequence = 0;
    profile1NumSeqs = 0;
}

// Andreas Wilm (UCD): shouldn't resetProfile1 and resetProfile2 be merged?
/** remove gaps from older alignments (code = gap_pos1) */
void Alignment::resetProfile1()
{
    register int sl;                /* which have  code = gap_pos2  */
    int i,j;
    int _gapPos1 = userParameters->getGapPos1();
    int _gapPos2 = userParameters->getGapPos2();
    bool _resetAlignNew = userParameters->getResetAlignmentsNew();
    bool _resetAlignAll = userParameters->getResetAlignmentsAll();
    
    if (userParameters->getStructPenalties1() != NONE) 
    {
        sl = 0;
        for (j = 0; j < (int)gapPenaltyMask1.size(); ++j) 
        {
            if (gapPenaltyMask1[j] == _gapPos1 && (_resetAlignNew ||
                           _resetAlignAll)) 
            {
                continue;
            }
            if (gapPenaltyMask1[j] == _gapPos2 && (_resetAlignAll)) 
            {
                continue;
            }
            gapPenaltyMask1[sl] = gapPenaltyMask1[j];
            ++sl;
        }
    }
  
    if (userParameters->getStructPenalties1() == SECST) 
    {
        sl = 0;
        for (j = 0; j < (int)secStructMask1.size(); ++j) 
        {
            if (secStructMask1[j] == _gapPos1 && (_resetAlignNew ||
                          _resetAlignAll)) 
            {
                continue;
            }
            if (secStructMask1[j] == _gapPos2 && (_resetAlignAll)) 
            {
                continue;
            }
            secStructMask1[sl] = secStructMask1[j];
            ++sl;
        }
    }
  
    for(i = 1; i <= profile1NumSeqs; ++i) 
    {
        sl = 0;
        for(j = 1; j <= getSeqLength(i); ++j) 
        {
            if(seqArray[i][j] == _gapPos1 && (_resetAlignNew || _resetAlignAll)) 
            {
                continue;
            }
            if(seqArray[i][j] == _gapPos2 && (_resetAlignAll)) 
            {
                continue;
            }
            ++sl;
            seqArray[i][sl] = seqArray[i][j];
        }

        // Andreas Wilm (UCD): added 2008-03-07:
        // Remove excess bit at end of sequence
        int numExtraElements = seqArray[i].size() - 1 - sl;
        for(int k = 0; k < numExtraElements; k++)
        {
            seqArray[i].pop_back();
        }
    }
}

// Andreas Wilm (UCD): shouldn't resetProfile1 and resetProfile2 be merged?
void Alignment::resetProfile2()
{
    register int sl;                /* which have  code = gap_pos2  */
    int i, j;
    int _gapPos1 = userParameters->getGapPos1();
    int _gapPos2 = userParameters->getGapPos2();
    bool _resetAlignNew = userParameters->getResetAlignmentsNew();
    bool _resetAlignAll = userParameters->getResetAlignmentsAll();
    int _profile1NumSeqs = profile1NumSeqs;

    
    if (userParameters->getStructPenalties2() != NONE) 
    {
        sl = 0;
        for (j = 0; j < (int)gapPenaltyMask2.size(); ++j) 
        {
            if (gapPenaltyMask2[j] == _gapPos1 && (_resetAlignNew || _resetAlignAll)) 
            {
                continue;
            }
            if (gapPenaltyMask2[j] == _gapPos2 && (_resetAlignAll)) 
            {
                continue;
            }
            gapPenaltyMask2[sl] = gapPenaltyMask2[j];
            ++sl;
        }
    }
  
    if (userParameters->getStructPenalties2() == SECST) 
    {
        sl = 0;
        for (j = 0; j < (int)secStructMask2.size(); ++j) 
        {
            if (secStructMask2[j] == _gapPos1 && (_resetAlignNew ||
                         _resetAlignAll)) 
            {
                continue;
            }
            if (secStructMask2[j] == _gapPos2 && (_resetAlignAll)) 
            {
                continue;
            }
            secStructMask2[sl] = secStructMask2[j];
            ++sl;
        }
    }
  
    for(i = _profile1NumSeqs + 1; i <= numSeqs; ++i) 
    {
        sl = 0;
        for(j = 1; j <= getSeqLength(i); ++j) 
        {
            if(seqArray[i][j] == _gapPos1 && (_resetAlignNew || _resetAlignAll)) 
            {
                continue;
            }
            if(seqArray[i][j] == _gapPos2 && (_resetAlignAll)) 
            {
                continue;
            }
            ++sl;
            seqArray[i][sl] = seqArray[i][j];
        }


        // Andreas Wilm (UCD) added 2008-03-07:
        // Remove excess bit at end of sequence
        int numExtraElements = seqArray[i].size() - 1 - sl;
        for(int k = 0; k < numExtraElements; k++)
        {
            seqArray[i].pop_back();
        }
    }

}

void Alignment::resetAllSeqWeights()
{
    seqWeight.clear();
    seqWeight.resize(numSeqs + 1, 100);
}
 
/*
 * The function addOutputIndex must check that the number of outputindex
 * Not sure what to do with this one.
 */
bool Alignment::addOutputIndex(vector<int>* outputIndexToAdd)
{
    // First need to clear the old outputIndex
    outputIndex.clear();
    // Check if the size is correct
    if((int)outputIndexToAdd->size() == numSeqs)
    {
        // Add the output index
        outputIndex = *outputIndexToAdd;
        return true;
    }
    else
    {
        clearAlignment();
        return false; // Could not add them
    }    
    
}

/*
 * The function appendOutputIndex is used when we are appending the outputIndex
 * to the current one. We do this when we are adding a profile2.
 */
bool Alignment::appendOutputIndex(vector<int>* outputIndexToAppend)
{
    // Check if the size is correct
    if((int)(outputIndexToAppend->size() + outputIndex.size()) == numSeqs)
    {
        // Append the outputIndex
        vector<int>::iterator outIndexIter = outputIndexToAppend->begin();
        while(outIndexIter != outputIndexToAppend->end())
        {
            outputIndex.push_back(*outIndexIter);
            outIndexIter++;
        }
        if((int)outputIndex.size() == numSeqs)
        {
            return true;
        }
        else
        {
            clearAlignment();
            cerr << "There is a problem with adding the sequences\n";
            return false;
        }
    }
    else
    {
        clearAlignment();
        return false;
    }
}

/*
 * The function addSequences takes a vector of sequences and adds it to the Alignment.
 * This is used if there is nothing in memory that we need.
 *
 */
void Alignment::addSequences(vector<Sequence>* seqVector)
{
    clearAlignment();
    
    numSeqs = seqVector->size();
    vector<int> emptyVec;
    
    // NOTE push dummy sequences on first!
    /********************************************************
     * It was decided to stay with the seqs and seqArray    *
     * index begining at 1. This is because it is difficult * 
     * to change in the code, so we push dummy seqs on      *
     ********************************************************/
    seqArray.push_back(emptyVec); // EMPTY sequence
    names.push_back(string(""));
    titles.push_back(string(""));
    sequenceIds.push_back(0);
    
    addSequencesToVector(seqVector);
    
    calculateMaxLengths();
    seqWeight.resize(numSeqs + 1, 100);
}

/*
 * The function appendSequences is used when we have a profile2.
 *
 */
void Alignment::appendSequences(vector<Sequence>* seqVector)
{
    numSeqs += seqVector->size(); // Add to the number we already have!
    addSequencesToVector(seqVector);
    seqWeight.clear();
    seqWeight.resize(numSeqs + 1, 100);
    calculateMaxLengths();
}

void Alignment::pasteSequencesIntoPosition(vector<Sequence>* seqVector, int pos, 
                                           bool explicitPasteToProfile2)
{
    SeqArray::iterator seqArrayIterator;
    vector<string>::iterator namesIterator;
    vector<string>::iterator titlesIterator;
    vector<unsigned long>::iterator sequenceIdsIterator;
    
    int profNum = userParameters->getProfileNum();
    int numSeqsToInsert = seqVector->size();
    if(numSeqsToInsert == 0 || pos < 0)
    {
        return;
    }
    if(pos == numSeqs)
    {
        seqArrayIterator = seqArray.end();
        namesIterator = names.end();
        titlesIterator = titles.end();
        sequenceIdsIterator = sequenceIds.end();
    }
    else
    {
        seqArrayIterator = seqArray.begin() + pos + 1;
        namesIterator = names.begin() + pos + 1;
        titlesIterator = titles.begin() + pos + 1;
        sequenceIdsIterator = sequenceIds.begin() + pos + 1;            
    }
    
    int prof1NumSeqs = profile1NumSeqs;

    for(int i = numSeqsToInsert - 1; i >= 0; i--)
    {
        seqArray.insert(seqArrayIterator, *(*seqVector)[i].getSequence());
        names.insert(namesIterator, (*seqVector)[i].getName());
        titles.insert(titlesIterator, (*seqVector)[i].getTitle());
        sequenceIds.insert(sequenceIdsIterator, (*seqVector)[i].getIdentifier());
        
        numSeqs++;
        if(profNum != 0 && !explicitPasteToProfile2 && pos <= prof1NumSeqs)
        {
            prof1NumSeqs++;
        }
    }
    
    if(profNum != 0 && pos <= prof1NumSeqs)
    {
        profile1NumSeqs = prof1NumSeqs;
    }
    
    resetAllSeqWeights();
    setDefaultOutputIndex();    
}

void Alignment::debugPrintAllNames()
{
    vector<string>::iterator nameIter = names.begin();
    while(nameIter != names.end())
    {
        cout << *nameIter << endl;
        nameIter++;
    }
}

void Alignment::NameIterator::begin(Alignment *a) {
    //cout << "begin()\n";
    if (a) {
        alignment = a;
        i = alignment->names.begin();
    }
}

const string Alignment::NameIterator::next() {
    //cout << "next()\n";
    if (!alignment)
        return string();
    if (i == alignment->names.end())
        return string();
    return *i++;
}

bool Alignment::NameIterator::end() {
    //cout << "end()\n";
    if (!alignment)
        return true;
    if (i == alignment->names.end())
        return true;
    return false;
}

// 26-03-03,nige: test before appending seqVector
bool Alignment::testUniqueNames(vector<Sequence>* seqVector, string *offendingSeq)
{
    vector<string>::iterator   oldName;
    vector<Sequence>::iterator newName;
    bool unique = true;
    
    //iterate over new candidate names
    for (newName = seqVector->begin(); unique && newName != seqVector->end(); newName++) {
        //iterate over old stored names
        for (oldName = names.begin() + 1; unique && oldName != names.end(); oldName++) {
            if (*oldName == newName->getName()) {
                offendingSeq->assign(*oldName);
                unique = false;
            }
        }
    }
    return unique;
}


/*
 * The function addSequencesToVector is used to add sequences to our seqVector
 * It is used by both addSequences and appendSequences. This is to reduce code
 * duplication.
 */
void Alignment::addSequencesToVector(vector<Sequence>* seqVector)
{
    std::vector<Sequence>::iterator itSeq; 
    
    for(itSeq = seqVector->begin(); itSeq != seqVector->end(); ++itSeq)
    {
        seqArray.push_back(*(*itSeq).getSequence());
        names.push_back((*itSeq).getName());
        titles.push_back((*itSeq).getTitle());
        sequenceIds.push_back((*itSeq).getIdentifier());
    }
    
    if(!(((int)seqArray.size() == numSeqs + 1) && ((int)names.size() == numSeqs + 1)
          && ((int)titles.size() == numSeqs + 1) && ((int)sequenceIds.size() == numSeqs + 1)))
    {
        cerr << "There has been an error adding the sequences to Alignment.\n"
             << "Must terminate the program. EaddSequencesrror occured in addSequences.\n";
        throw 1;
    }
}

void Alignment::clearAlignment()
{
    // Erase all the elements from the vector!
    clearSeqArray();
    names.clear();
    titles.clear();
    outputIndex.clear();
    sequenceIds.clear();
    clearSecStruct1();
    clearSecStruct2();
    seqWeight.clear();
    maxNames = 0;
    numSeqs = 0;
    maxAlignmentLength = 0;
    lengthLongestSequence = 0; 
    userParameters->setProfileNum(0);
    userParameters->setProfile1Empty(true);
    userParameters->setProfile2Empty(true);   
}

void Alignment::clearSeqArray()
{
    for(int i = 0; i < (int)seqArray.size(); i++)
    {
        seqArray[i].clear();
    }
    seqArray.clear();
}

void Alignment::clearSecStruct1()
{
    gapPenaltyMask1.clear();
    secStructMask1.clear();
    secStructName1 = "";
}

void Alignment::clearSecStruct2()
{
    gapPenaltyMask2.clear();
    secStructMask2.clear();
    secStructName2 = "";
}

/*
 * The function alignScore is used to score the alignment!
 * This is used by other classes to see what score the alignment has gotten.
 */
int Alignment::alignScore(void)
{
    int maxRes;
    int seq1, seq2, res1, res2;
    int ngaps;
    int i, len1, len2;
    int score;
    int _maxAA = userParameters->getMaxAA();
    float _gapOpen = userParameters->getGapOpen();
    int matrix[NUMRES][NUMRES];
    
    //
    // calculate an overall score for the alignment by summing the
    // scores for each pairwise alignment
    //
    maxRes = subMatrix->getAlnScoreMatrix(matrix);
    if (maxRes == 0)
    {
        utilityObject->error("Matrix for alignment scoring not found\n");
        return 0;
    }

    score = 0;
    for (seq1 = 1; seq1 <= numSeqs; seq1++) 
    {
        for (seq2 = 1; seq2 < seq1; seq2++)
        {
            len1 = seqArray[seq1].size() - 1; 
            len2 = seqArray[seq2].size() - 1; 
            for (i = 1; i < len1 && i < len2; i++)
            {
                res1 = seqArray[seq1][i];
                res2 = seqArray[seq2][i];

                if ((res1 >= 0) && (res1 <= _maxAA) && (res2 >= 0) && (res2 <= _maxAA))
                {
                    score += matrix[res1][res2];
                }
            }

            ngaps = countGaps(seq1, seq2, len1);

            score = static_cast<int>(score - (100 * _gapOpen * ngaps)); // Mark change 17-5-07

        }
    }

    score /= 100;

    utilityObject->info("Alignment Score %d\n", score);    
    return score;
}

int Alignment::countGaps(int seq1, int seq2, int len)
{
    int i, g;
    int q, r;//,  *Q,  *R;
 
    vector<int> Q, R;
    Q.resize(len + 2, 0);
    R.resize(len + 2, 0);
    
    int _maxAA = userParameters->getMaxAA();
    
    try
    {
        Q[0] = R[0] = g = 0;

        for (i = 1; i < len; i++)
        {
            if (seqArray[seq1][i] > _maxAA)
            {
                q = 1;
            }
            else
            {
                q = 0;
            }

            if (seqArray[seq2][i] > _maxAA)
            {
                r = 1;
            }
            else
            {
                r = 0;
            }
        
            // NOTE I havent a clue what this does!!!!
            if (((Q[i - 1] <= R[i - 1]) && (q != 0) && (1-r != 0)) || 
                ((Q[i - 1] >= R[i - 1]) && (1-q != 0) && (r != 0)))
            {
                g += 1;
            }

            if (q != 0)
            {
                Q[i] = Q[i - 1] + 1;
            }
            else
            {
                Q[i] = 0;
            }

            if (r != 0)
            {
                R[i] = R[i - 1] + 1;
            }
            else
            {
                R[i] = 0;
            }
        }
    }
    catch(const exception &ex)
    {
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue. Function = countGaps\n";
        throw 1;
    }

    return (g);
}

void Alignment::resetAlign()
{
    /* remove gaps from older alignments (code = gap_pos1) EXCEPT for
       gaps that were INPUT with the seqs.  which have code =
       gap_pos2
    */
    
    register int sl;
    int i, j;
    int _gapPos1 = userParameters->getGapPos1();
    int _gapPos2 = userParameters->getGapPos2();
    bool _resetAlignNew = userParameters->getResetAlignmentsNew();
    bool _resetAlignAll = userParameters->getResetAlignmentsAll();

    
    for(i = 1; i <= numSeqs;++i) 
    {
        sl = 0;
        for(j = 1; j <= getSeqLength(i); ++j) 
        {
            if(seqArray[i][j] == _gapPos1 && ( _resetAlignNew || _resetAlignAll)) 
            {
                continue;
            }
            if(seqArray[i][j] == _gapPos2 && (_resetAlignAll)) 
            {
                continue;
            }
            ++sl;
            seqArray[i][sl] = seqArray[i][j];
        }

        
        // Andreas Wilm (UCD) added 2008-03-07:
        // Remove excess bit at end of sequence
        int numExtraElements = seqArray[i].size() - 1 - sl;
        for(int k = 0; k < numExtraElements; k++)
        {
            seqArray[i].pop_back();
        }
    }
}

void Alignment::fixGaps()
{
    int i,j;
  
    if (userParameters->getStructPenalties1() != NONE) 
    {
        for (j = 0; j < getSeqLength(1); ++j) 
        {
            if (gapPenaltyMask1[j] == userParameters->getGapPos1())
            {
                gapPenaltyMask1[j] = userParameters->getGapPos2();
            }
        }
    }
  
    if (userParameters->getStructPenalties1() == SECST) 
    {
        for (j = 0; j < getSeqLength(1); ++j) 
        {
            if (secStructMask1[j] == userParameters->getGapPos1())
            {
                secStructMask1[j] = userParameters->getGapPos2();
            }
        }
    }
  
    for(i = 1; i <= numSeqs; ++i) 
    {
        for(j = 1; j <= getSeqLength(i); ++j) 
        {
            if(seqArray[i][j] == userParameters->getGapPos1())
            {
                seqArray[i][j] = userParameters->getGapPos2();
            }
        }
    }

}

float Alignment::countid(int s1, int s2)
{
    int c1, c2;
    int i;
    int count, total;
    float score;
    int shorterSeqLength = (getSeqLength(s1) < getSeqLength(s2)) ? getSeqLength(s1) : getSeqLength(s2);

    count = total = 0;
    for (i = 1; i <= shorterSeqLength; i++) // NOTE june29
    {
        c1 = seqArray[s1][i];
        c2 = seqArray[s2][i];
        if ((c1 >= 0) && (c1 < userParameters->getMaxAA()))
        {
            total++;
            if (c1 == c2)
            {
                count++;
            }
        }

    }

    if (total == 0)
    {
        score = 0;
    }
    else
    {
        score = 100.0 *(float)count / (float)total;
    }

    return (score);
}

void Alignment::debugPrintSequences()
{
    cout << std::endl;
    for(int i = 0; i < (int)seqArray.size(); i++)
    {
        for(int j = 0; j < (int)seqArray[i].size(); j++)
        {
            if(seqArray[i][j] > 9)
                cout << " " << seqArray[i][j];
            else
                cout << "  " << seqArray[i][j];
        }
        cout << std::endl;
    } 
}

/*
 * Note the max_aln_length is now called maxAlignmentLength, and it will be stored
 * and calculated in this class. Mostly it is used for allocating arrays. But not always.
 */        
void Alignment::calculateMaxLengths()
{
    maxAlignmentLength = 0;
    lengthLongestSequence = 0;
    if(seqArray.size() > 0)
    {
        SeqArray::iterator seqArrayIter = seqArray.begin();
        
        while(seqArrayIter != seqArray.end())
        {
            // NOTE I needed to change this to >= for a bug I had!!!!!!!
            if((int)(*seqArrayIter).size() - 1 >= lengthLongestSequence)
            {
                lengthLongestSequence = (*seqArrayIter).size();
            }
            ++seqArrayIter;
        }

        if(lengthLongestSequence > 0)
        {
            maxAlignmentLength = (lengthLongestSequence * 2) - 2;
            lengthLongestSequence -= 1; // MADE A CHANGE HERE AS WELL!! 
        }
        else
        {
            lengthLongestSequence = 0;
            maxAlignmentLength = 0;
        }

    }
    else
    {
        maxAlignmentLength = 0;
    }
    maxNames = 0;
    if(names.size() > 0)
    {
        vector<string>::iterator nameVecIter = names.begin();
        
        while(nameVecIter != names.end())
        {
            if((int)(*nameVecIter).size() > maxNames)
            {
                maxNames = (*nameVecIter).size();
            }
            ++nameVecIter;
        }
        if(maxNames < 10)
        {
            maxNames = 10;
        }
    }
    else
    {
        maxNames = 0;
    }
}

/*
 * This function checks to see if all names are different. It returns true if they
 * are all different, and false if there are 2 the same.
 * currently quadratic - would be nice to speed up
 */
bool Alignment::checkAllNamesDifferent(string *offendingSeq)
{
    bool different = true;
    // NOTE I added the + 1 here because, if we had a sequence with a name as blank
    // this would be the same as the first one!
    vector<string>::iterator namesIter1 = names.begin() + 1;
    vector<string>::iterator namesIter2;
    int counter1 = 1;
    int counter2 = 2;

    
    while(namesIter1 != names.end())
    {
        namesIter2 = namesIter1 + 1;
        while(namesIter2 != names.end())
        {
            if((*namesIter1).compare(*namesIter2) == 0) // If we have 2 strings the same.
            {
                different = false;     
                /* 23-03-2007,nige: let someone up the stack deal with this - GUI is too deeply entangled.
                 * utilityObject->error("Multiple sequences found with same name '%s' (first %d chars are significant)\n", namesIter1->c_str(), MAXNAMES);
                 */

                offendingSeq->assign((*namesIter1));
                clearAlignment();
                return different; // Not all different!
            }
            namesIter2++;
            counter2++;
        }
        namesIter1++;
        counter1++;
        counter2 = counter1 + 1;
    }
    return different;
}






void Alignment::addSecStructMask1(vector<char>* secStructMaskToAdd)
{
    secStructMask1 = *secStructMaskToAdd;
}

void Alignment::addSecStructMask2(vector<char>* secStructMaskToAdd)
{
    secStructMask2 = *secStructMaskToAdd;
}

void Alignment::addGapPenaltyMask1(vector<char>* gapPenaltyMaskToAdd)
{
    gapPenaltyMask1 = *gapPenaltyMaskToAdd;
}

void Alignment::addGapPenaltyMask2(vector<char>* gapPenaltyMaskToAdd)
{
    gapPenaltyMask2 = *gapPenaltyMaskToAdd;
}

void Alignment::addSecStructName1(string nameToAdd)
{
    secStructName1 = nameToAdd;
}

void Alignment::addSecStructName2(string nameToAdd)
{
    secStructName2 = nameToAdd;
}

void Alignment::addSeqWeight(vector<int>* _seqWeight)
{
    if(seqWeight.size() == _seqWeight->size())
    {
        int size = seqWeight.size();
        
        for(int i = 0; i < size; i++)
        {
            seqWeight[i] = (*_seqWeight)[i];
        }
    }
}

void Alignment::printSequencesAddedInfo()
{
    if(userParameters->getDisplayInfo())
    {
        int startSeq = userParameters->getProfile2Empty() ? 1:
                       profile1NumSeqs + 1;
        string dnaFlag = userParameters->getDNAFlag() ? "bp" : "aa";


        for(int i = startSeq; i <= numSeqs; i++) 
        {
            cout << "Sequence " << i << ": " 
                 << std::left << setw(maxNames) << names.at(i) 
                 << std::right << setw(6) << getSequenceLength(i) 
                 << " " << dnaFlag << std::endl;
        }
    }
}

void Alignment::debugPrintOutAlignInfo()
{
    for(int i = 1; i <= numSeqs; i++) 
    {
        cout << "seq-no=" << i << ": name=" 
             << std::left << setw(maxNames) << names.at(i)
             << " length=" 
             << std::right << setw(6) << getSequenceLength(i) 
             << std::endl;
    }
}

int Alignment::getSequenceLength(int index)
{
    return seqArray.at(index).size() - 1;
}

int Alignment::getLengthLongestSequence()
{
    int _lengthLongestSequence = 0;

    for(int i = 1; i <= numSeqs; i++)
    {
        if(getSeqLength(i) > _lengthLongestSequence)
        {
            _lengthLongestSequence = getSeqLength(i);
        }
    }    
    return _lengthLongestSequence;
}

/**
 * This function is used to find the length of the longest sequence in the range.
 */
int Alignment::getLengthLongestSequence(int firstSeq, int lastSeq)
{
    int _lengthLongestSequence = 0;
    
    if(firstSeq >= 1 && lastSeq <= numSeqs)
    {
        for(int i = firstSeq; i <= lastSeq; i++)
        {
            if(getSeqLength(i) > _lengthLongestSequence)
            {
                _lengthLongestSequence = getSeqLength(i);
            }
        }
    }    
    return _lengthLongestSequence; // Will return 0 if cant check seqs
}

int Alignment::getMaxNames()
{
    return maxNames;
}

string Alignment::getSecStructName1()
{
    return secStructName1;
}

string Alignment::getSecStructName2()
{
    return secStructName2;
}

vector<char>* Alignment::getGapPenaltyMask1()
{
    return &gapPenaltyMask1;
}

vector<char>* Alignment::getGapPenaltyMask2()
{
    return &gapPenaltyMask2;
}

vector<char>* Alignment::getSecStructMask1()
{
    return &secStructMask1;
}

vector<char>* Alignment::getSecStructMask2()
{
    return &secStructMask2;
}


int Alignment::getSecStructMask1Element(int index)
{
    if(index > 0 && index < (int)secStructMask1.size())
    {
        return secStructMask1[index];
    }
    else
    {
        throw VectorOutOfRange(string("secStructMask1"), index, secStructMask1.size() - 1);
    }    
}

int Alignment::getSecStructMask2Element(int index)
{
    if(index > 0 && index < (int)secStructMask2.size())
    {
        return secStructMask2[index];
    }
    else
    {
        throw VectorOutOfRange(string("secStructMask2"), index, secStructMask2.size() - 1);
    } 
}

int Alignment::getGapPenaltyMask1Element(int index)
{
    if(index > 0 && index < (int)gapPenaltyMask1.size())
    {
        return gapPenaltyMask1[index];
    }
    else
    {
        throw VectorOutOfRange(string("gapPenaltyMask1"), index, gapPenaltyMask1.size() - 1);
    }   
}

int Alignment::getGapPenaltyMask2Element(int index)
{
    if(index > 0 && index < (int)gapPenaltyMask2.size())
    {
        return gapPenaltyMask2[index];
    }
    else
    {
        throw VectorOutOfRange(string("gapPenaltyMask2"), index, gapPenaltyMask2.size() - 1);
    }   
}

string Alignment::getName(int index)
{
    if(index > 0 && index < (int)names.size())
    {
        return names[index];
    }
    else
    {
        throw VectorOutOfRange(string("names"), index, names.size() - 1);
    } 
}

/**
 * Change:
 * Mark 25-1-2007: This function was added to allow access to the unique id.
 */
unsigned long Alignment::getUniqueId(int seq)
{
    if(seq > 0 && seq < (int)sequenceIds.size())
    {
        return sequenceIds[seq];
    }
    else
    {
        throw VectorOutOfRange(string("sequenceIds"), seq, sequenceIds.size() - 1);
    } 
}

string Alignment::getTitle(int index)
{
    if(index > 0 && index < (int)titles.size())
    {
        return titles[index];
    }
    else
    {
        throw VectorOutOfRange(string("titles"), index, titles.size() - 1);
    } 
}

int Alignment::getOutputIndex(int index)
{
    if(index >= 0 && index < (int)outputIndex.size())
    {
        return outputIndex[index];
    }
    else
    {
        throw VectorOutOfRange(string("outputIndex"), index, outputIndex.size() - 1);
    } 
}

int Alignment::getSeqWeight(int index) const
{
    if(index >= 0 && index < (int)seqWeight.size())
    {
        return seqWeight[index];
    }
    else
    {
        throw VectorOutOfRange(string("seqWeight"), index, seqWeight.size() - 1);
    } 
}


/**
 * This function is used by Qt. It is used to calculate the histogram values for the
 * widget in clustalQt. The function is in here because it needs access to the SubMatrix
 * class.
 * @param matNum The number of the matrix to be used in the comparison.
 * @return A pointer to a vector containing the histogram values.
 */
vector<int>* Alignment::QTcalcHistColumnHeights(int firstSeq, int nSeqs,
                                              Array2D<int>* exceptionalRes)
{
    int n, i, s, p, r, r1;
    int numColumns = getLengthLongestSequence();
    int scoreScale = userParameters->getQTScorePlotScale();//5; // From ClustalX.
    int scoreCutOff = userParameters->getQTResExceptionCutOff();
    bool includeGaps = false;
    //short  *mat_xref, *matptr;
    float median, mean;
    float t, q1, q3, ul;
    vector<float> seqdist, sorteddist;
    float diff;
    vector<int> seqvector;
    vector<int> freq;
    vector<vector<int> > profile;
    int matrix[NUMRES][NUMRES];
    int _maxAA = userParameters->getMaxAA();
    int _gapPos1 = userParameters->getGapPos1();
    histogramColumnHeights.resize(numColumns);
    //panel_data data1;
    subMatrix->getQTMatrixForHistogram(matrix);
    
    profile.resize(numColumns + 2, vector<int>(_maxAA + 2));
    freq.resize(_maxAA + 2);
 
    for(p = 0; p < numColumns; p++)
    {
        for(r = 0; r < _maxAA; r++)
        {
            freq[r] = 0;
        }        
        for(s = firstSeq; s < firstSeq + nSeqs; s++)
        {
            if(p < getSeqLength(s + 1) && seqArray[s + 1][p + 1] >= 0 &&
                 seqArray[s + 1][p + 1] < _maxAA)
            {
                freq[seqArray[s + 1][p + 1]]++;
            }
        }
        for(r = 0; r < _maxAA; r++)
        {
            profile[p][r] = 0;
            for(r1 = 0; r1 < _maxAA; r1++)
            {
                profile[p][r] += freq[r1] * matrix[r1][r];
            }
            profile[p][r] = static_cast<int>(profile[p][r] / 
                             static_cast<float>(nSeqs)); // Mark change 17-5-07
        }
    }

    seqvector.resize(_maxAA + 2);
    seqdist.resize(nSeqs + 1);
    sorteddist.resize(nSeqs + 1);

    for(p = 0; p < numColumns; p++)
    {
        for(s = firstSeq; s < firstSeq + nSeqs; s++)
        {
            if (p < getSeqLength(s + 1))
            {
                for (r = 0; r < _maxAA; r++)
                {
                    seqvector[r] = matrix[r][seqArray[s + 1][p + 1]];
                }
            }
            else
            {
                for (r = 0; r < _maxAA; r++)
                {
                    seqvector[r] = matrix[r][_gapPos1];
                }
            }
            seqdist[s - firstSeq] = 0.0;
            for(r = 0; r < _maxAA; r++)
            {
                diff = profile[p][r] - seqvector[r];
                diff /= 1000.0;
                seqdist[s - firstSeq] += diff * diff;
            }
            seqdist[s - firstSeq] = sqrt((double)seqdist[s - firstSeq]);
        }

        // calculate mean,median and rms of seq distances
        mean = median = 0.0;
        if(includeGaps)
        {
            for(s = 0; s < nSeqs; s++)
            {
                mean += seqdist[s];
            }
            mean /= nSeqs;
            n = nSeqs;
            for(s = 0; s < nSeqs; s++)
            {
                sorteddist[s] = seqdist[s];
            }
        }
        else
        {
            n = 0;
            for(s = firstSeq; s < firstSeq + nSeqs; s++)
            {    
                if(p < getSeqLength(s + 1) && seqArray[s + 1][p + 1] >= 0 &&
                   seqArray[s + 1][p + 1] < _maxAA)
                {
                    mean += seqdist[s - firstSeq];
                    n++;
                }
            }
            if(n > 0) 
            {
                mean /= n;
            }
            for(s = firstSeq, i = 0; s < firstSeq + nSeqs; s++)
            {
                if(p < getSeqLength(s + 1) && seqArray[s + 1][p + 1] >= 0 &&
                   seqArray[s + 1][p + 1] < _maxAA)
                {
                    sorteddist[i++] = seqdist[s - firstSeq];
                }
            }
        }
        sortScores(&sorteddist, 0, n - 1);


        if(n == 0)
        {
            median = 0;
        }
        else if(n % 2 == 0)
        {
            median = (sorteddist[n / 2 - 1] + sorteddist[n / 2]) / 2.0;
        }
        else
        {
            median = sorteddist[n / 2];
        }
        
        if(scoreScale <= 5)
        {
            histogramColumnHeights[p] = static_cast<int>(exp((double)(-mean *
                                        (6 - scoreScale) / 4.0)) * 100.0 * n / nSeqs);
        }
        else
        {    
            histogramColumnHeights[p] = static_cast<int>(exp((double)(-mean / 
                                        (4.0 * (scoreScale - 4)))) * 100.0 * n / nSeqs);
        }
        
        if(n == 0)
        {
            ul = 0;
        }
        else
        {
            t = n/4.0 + 0.5;
            if(t - (int)t == 0.5)
            {
                q3 = (sorteddist[(int)t] + sorteddist[(int)t + 1]) / 2.0;
                q1 = (sorteddist[n-(int)t] + sorteddist[n - (int)t - 1]) / 2.0;
            }
            else if(t - (int)t > 0.5)
            {
                q3 = sorteddist[(int)t + 1];
                q1 = sorteddist[n - (int)t - 1];
            }
            else 
            {
                q3 = sorteddist[(int)t];
                q1 = sorteddist[n - (int)t];
            }
            if (n < 4)
            {
                ul = sorteddist[0];
            }
            else
            { 
                ul = q3 + (q3 - q1) * ((float)scoreCutOff / 2.0);
            }
        }
        
        if((exceptionalRes->getRowSize() >= nSeqs) &&
            exceptionalRes->getColSize() >= numColumns)
        {
            for(int s = firstSeq; s < firstSeq + nSeqs; s++)
            {
                if(seqdist[s - firstSeq] > ul && p < getSeqLength(s + 1) && 
                   seqArray[s + 1][p + 1] >= 0 && 
                   seqArray[s + 1][p + 1] < userParameters->getMaxAA())
                {
                    (*exceptionalRes)[s - firstSeq][p] = 1;
                }
                else
                {
                    (*exceptionalRes)[s - firstSeq][p] = 0;
                }
            }
        }
    }
    return &histogramColumnHeights;
}

void Alignment::sortScores(vector<float>* scores, int f, int l)
{
    int i,last;

    if(f >= l) 
        return;

    swap(scores, f, (f + l) / 2);
    last = f;
    for(i = f + 1; i <= l; i++)
    {
        if((*scores)[i] > (*scores)[f])
        {
            swap(scores, ++last, i);
        }
    }
    swap(scores, f, last);
    sortScores(scores, f, last - 1);
    sortScores(scores, last + 1, l);
}

void Alignment::swap(vector<float>* scores, int s1, int s2)
{
    float temp;

    temp = (*scores)[s1];
    (*scores)[s1] = (*scores)[s2];
    (*scores)[s2] = temp;
}

int Alignment::searchForString(bool* found, int seq, int beginRes, string search)
{
    int lengthSeq = getSeqLength(seq);
    if(beginRes > lengthSeq)
    {
        *found = false;
        return beginRes;
    }
    
    int res = beginRes;
    
    // First need to convert search into a vector of ints!
    vector<int> codedSearch;
    int size = search.size();
    codedSearch.resize(size);
    int code;
    for(int i = 0; i < size; i++)
    {
        code = userParameters->resIndex(userParameters->getAminoAcidCodes(), search[i]);
        codedSearch[i] = code;
    }

    int numSame = 0;
    int startPos = -1;
    int searchSize = codedSearch.size();
    // Now check for the string of ints!!!!
    for(; res <= lengthSeq; res++) //- nige: res starts at 1

    {
        if(seqArray[seq][res] == codedSearch[0])
        {
            startPos = res; //- nige
            for(int i = 0; i < searchSize && (i + res) <= lengthSeq; i++) //- nige: res starts at 1
            {
                if(seqArray[seq][res + i] == codedSearch[i])
                {
                    numSame++;
                }
                //nige: hack: encoded gap character: see also AlignmentScroll.cpp
                else if (seqArray[seq][res + i] == 31 || seqArray[seq][res + i] == 30)
                {
                    res++;
                    i--;
                }
                else
                {
                    break; // Not the same
                }
            }
            if(numSame == searchSize)
            {
                *found = true;
                return startPos;
            }
            else
            {
                numSame = 0;
            }
        }
    }
    *found = false;
    return startPos; // Not found!!!!!!!
}

void Alignment::removeGapsFromSelectedSeqs(vector<int>* selected)
{
    //getNumSeqs()
    int size = selected->size();
    int lengthOfSelectedSeq = 0;
    int gapPos1 = userParameters->getGapPos1();
    int gapPos2 = userParameters->getGapPos2();
    int s1;
    
    for(int i = 1; i <= getNumSeqs() && i < size; i++)
    {
        if((*selected)[i] == 1)
        {
            // remove gaps from this seq!
            lengthOfSelectedSeq = getSeqLength(i);
            s1 = 0;
            for(int j = 1; j <= lengthOfSelectedSeq; j++)
            {
                if(seqArray[i][j] == gapPos1 || seqArray[i][j] == gapPos2)
                {
                    continue;
                }
                ++s1;
                seqArray[i][s1] = seqArray[i][j];
            }
            // Then remove the excess bit at the end of the array
            int numExtraElements = lengthOfSelectedSeq - s1;
            
            if((int)seqArray[i].size() > numExtraElements)
            {
                for(int k = 0; k < numExtraElements; k++)
                {
                    seqArray[i].pop_back();
                }
            }
        }
    }

}

void Alignment::removeGapOnlyColsFromSelectedSeqs(vector<int>* selected)
{
    int numGaps = 0;
    int NoneSelected = -1;
    int numColumns = 0;
    int sizeSelected = selected->size();
    int firstSeqSelected = NoneSelected;
    int gapPos1 = userParameters->getGapPos1();
    int gapPos2 = userParameters->getGapPos2();
    int k; 
      
    for(int i = 1; i < sizeSelected; i++)
    {
        if((*selected)[i] == 1)
        {
            numColumns++;
            if(firstSeqSelected == -1)
            {
                firstSeqSelected = i;
            }
        }
    }
    
    if(firstSeqSelected == NoneSelected)
    {
        cout << "No Sequences have been selected\n";
        return;
    }
    
    for(int i = 1; i <= getSeqLength(firstSeqSelected);)
    {
        numGaps = 0;
        for(int j = firstSeqSelected; j < sizeSelected && (*selected)[j] == 1; j++)
        {
            if(getSeqLength(j) >= i)
            {
                if(seqArray[j][i] == gapPos1 || seqArray[j][i] == gapPos2)
                {
                    numGaps++;
                }
            }
        }
        if(numGaps == numColumns)
        {
            //cout << "                removing a gap column\n\n";
            for(int j = firstSeqSelected; j < sizeSelected && (*selected)[j] == 1; j++)
            {
                for(k = i + 1; k <= getSeqLength(j) + 1 && (int)seqArray[j].size() > k; k++)
                {
                    seqArray[j][k - 1] = seqArray[j][k];
                }
                seqArray[j].pop_back(); // Remove the last element!!
                
                if(getSeqLength(firstSeqSelected) <= 0)
                {
                    break;
                }                
            }
        }
        else
        {
            i++;
        }
    }
    
}

void Alignment::removeAllGapOnlyColumns(int fSeq, int lSeq, int profileNum)
{
    if(fSeq >= lSeq)
    {
        return;
    }
    int gapPos1 = userParameters->getGapPos1();
    int gapPos2 = userParameters->getGapPos2();
    
    int numGaps = 0;
    int numColumns = lSeq - fSeq + 1;
    int k;
    // We must cheack each column to see if it consists of only '-'
    for(int i = 1; i <= getSeqLength(fSeq);)
    {
        numGaps = 0;
        for(int j = fSeq; j <= lSeq; j++)
        {
            if(getSeqLength(j) >= i)
            {
                if(seqArray[j][i] == gapPos1 || seqArray[j][i] == gapPos2)
                {
                    numGaps++;
                }
            }
        }
        if(numGaps == numColumns)
        {
            for(int j = fSeq; j <= lSeq; j++)
            {
                for(k = i + 1; k <= getSeqLength(j) + 1; k++)
                {
                    seqArray[j][k - 1] = seqArray[j][k];
                }
                seqArray[j].pop_back(); // Remove the last element!!
                
                if(profileNum == 1)
                {
                    int lengthSecStruct = secStructMask1.size();
                    int lengthGapMask = gapPenaltyMask1.size();
                    
                    for(k = i; k <= getSeqLength(fSeq) && k < lengthSecStruct; k++)
                    {
                        secStructMask1[k - 1] = secStructMask1[k];
                    }
                    for(k = i; k <= getSeqLength(fSeq) && k < lengthGapMask; k++)
                    {
                        gapPenaltyMask1[k - 1] = gapPenaltyMask1[k];
                    }                    
                }
                
                if(profileNum == 2)
                {
                    int lengthSecStruct = secStructMask2.size();
                    int lengthGapMask = gapPenaltyMask2.size();
                    
                    for(k = i; k <= getSeqLength(fSeq) && k < lengthSecStruct; k++)
                    {
                        secStructMask2[k - 1] = secStructMask2[k];
                    }
                    for(k = i; k <= getSeqLength(fSeq) && k < lengthGapMask; k++)
                    {
                        gapPenaltyMask2[k - 1] = gapPenaltyMask2[k];
                    }                    
                }
                if(getSeqLength(fSeq) <= 0)
                {
                    break;
                }                
            }
        }
        else
        {
            i++;
        }
    }
}

vector<Sequence> Alignment::cutSelectedSequencesFromAlignment(vector<int>* selected)
{
    vector<Sequence> cutSequences;
    int sizeOfSelected = selected->size();
    SeqArray::iterator seqArrayIterator;
    vector<string>::iterator namesIterator;
    vector<string>::iterator titlesIterator;
    vector<int>::iterator seqWeightIterator;
    vector<unsigned long>::iterator sequenceIdsIterator;
     
    int newProfile1NumSeqs = profile1NumSeqs;
    int profNum = userParameters->getProfileNum();
    int prof1NumSeqs = profile1NumSeqs; 
    int numCutSoFar = 0;
    int intialNumSeqs = numSeqs;
      
    for(int i = 1; i < sizeOfSelected && i <= intialNumSeqs; i++)
    {
        if((*selected)[i] == 1)
        {
            // Cut the sequence from the alignment!
            seqArrayIterator = seqArray.begin() + i - numCutSoFar;
            namesIterator = names.begin() + i - numCutSoFar;
            titlesIterator = titles.begin() + i - numCutSoFar;
            seqWeightIterator = seqWeight.begin() + i - numCutSoFar;
            sequenceIdsIterator = sequenceIds.begin() + i - numCutSoFar;
            Sequence SeqToCut(&seqArray[i - numCutSoFar], *namesIterator, *titlesIterator, 
                              *sequenceIdsIterator);
            
            numCutSoFar++;
            
            seqArrayIterator->clear();
            seqArray.erase(seqArrayIterator);
            names.erase(namesIterator);
            titles.erase(titlesIterator);
            seqWeight.erase(seqWeightIterator);
            sequenceIds.erase(sequenceIdsIterator);
            
            if(numSeqs > 0)
            {
                numSeqs--;
            }
            else
            {
                numSeqs = 0;
                break;
            }
            if(profNum > 0 && i <= prof1NumSeqs)
            {
                if(newProfile1NumSeqs > 0)
                {
                    newProfile1NumSeqs--;
                }
                else
                {
                    newProfile1NumSeqs = 0;
                }
            }
            cutSequences.push_back(SeqToCut);
        }
    }
    profile1NumSeqs = newProfile1NumSeqs;
    setDefaultOutputIndex();
    resetAllSeqWeights();
    return cutSequences; 
}

void Alignment::setDefaultOutputIndex()
{
    outputIndex.clear();
    outputIndex.resize(numSeqs);
    for(int iseq = 1; iseq <= numSeqs; iseq++)
    {
        outputIndex[iseq - 1] = iseq;
    }
}

bool Alignment::removeAllOutsideRange(int beginPos, int endPos)
{
    bool ok;
    if(beginPos < 0 || endPos > getLengthLongestSequence())
    {
        return false; // cannot do it!!!!
    }
    
    // trim the seqArray
    ok = keepPortionOfSeqArray(beginPos, endPos);
    if(!ok)
    {
        cerr << "There was a problem removing a portion of the array\n";
        return false;
    }
    
    // recalculate the maxLengths
    calculateMaxLengths();
    
    // Clear the histogram columns
    histogramColumnHeights.clear();
    
    // reset the weights
    resetAllSeqWeights();
	return true;
}

bool Alignment::keepPortionOfSeqArray(int beginRangeIndex, int endRangeIndex)
{
    SeqArray sectionToRealign;
    vector<int> emptyVec;
    sectionToRealign.push_back(emptyVec); // EMPTY sequence 
       
    SeqArray::iterator posToAddTo = sectionToRealign.begin(); 
    // erase from all sequences the range specified here!!!!!
    if(beginRangeIndex < 0 || endRangeIndex < 0)
    {
        return false;
    }  
      
    SeqArray::iterator mainBeginIt = seqArray.begin() + 1;
    SeqArray::iterator mainEndIt = seqArray.end();
    
    vector<int>::iterator begin, end, beginRange, endRange, beginCopyRange, endCopyRange;
    
    for(; mainBeginIt != mainEndIt; mainBeginIt++)
    {
        vector<int> vecToAdd;
        begin = mainBeginIt->begin() + 1;
        end = mainBeginIt->end();
        beginRange = begin + beginRangeIndex;
        endRange = begin + endRangeIndex + 1;
        beginCopyRange = beginRange;
        endCopyRange = endRange;
        
        // We need to copy all of this into another vector.
        if(endCopyRange < end && beginCopyRange < end)
        {
            vecToAdd.push_back(0);
            for(; beginCopyRange != endCopyRange; beginCopyRange++)
            {
                vecToAdd.push_back(*beginCopyRange);
            }
            sectionToRealign.push_back(vecToAdd);
        }
        else
        {
            return false;        
        }
        
        if(endRange < end && beginRange < end)
        {
            mainBeginIt->erase(beginRange, endRange);
        }
        else
        {
            return false;       
        }
    }
    clearSeqArray();
    seqArray = sectionToRealign;    
    return true; 
}

void Alignment::debugPrintSeqArray(SeqArray* arrayToPrint)
{
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
}

void Alignment::debugPrintProfile1()
{
    cout << "************** PROFILE1 *********************\n";
    SeqArray::iterator mainBeginIt = seqArray.begin() + 1;
    SeqArray::iterator mainEndIt = mainBeginIt + profile1NumSeqs;
    vector<int>::iterator begin, end;
    string aaCodes = userParameters->getAminoAcidCodes();
            
    for(; mainBeginIt != mainEndIt; mainBeginIt++)
    {
        cout << "PROFILE1 SEQ: ";
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
}

void Alignment::debugPrintProfile2()
{
    cout << "************** PROFILE2 *********************\n";
    SeqArray::iterator mainBeginIt = seqArray.begin() + 1 + profile1NumSeqs;
    SeqArray::iterator mainEndIt = seqArray.end();
    vector<int>::iterator begin, end;
    string aaCodes = userParameters->getAminoAcidCodes();
            
    for(; mainBeginIt != mainEndIt; mainBeginIt++)
    {
        cout << "PROFILE2 SEQ: ";
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
}
                             
bool Alignment::updateRealignedRange(SeqArray realignedSeqs, int beginPos, int endPos)
{
    if(realignedSeqs.size() != seqArray.size())
    {
        return false;
    }
    if(beginPos < 0 || endPos < 0)
    {
        return false;
    }
    
    // erase from all sequences the range specified here!!!!!  
      
    SeqArray::iterator mainBeginIt = seqArray.begin() + 1;
    SeqArray::iterator mainEndIt = seqArray.end();

    SeqArray::iterator pasteBeginIt = realignedSeqs.begin() + 1;
    SeqArray::iterator pasteEndIt = realignedSeqs.end();
            
    vector<int>::iterator begin, end, beginRange, endRange;
    
    for(; mainBeginIt != mainEndIt && pasteBeginIt != pasteEndIt; mainBeginIt++)
    {
        vector<int> vecToAdd;
        begin = mainBeginIt->begin() + 1;
        end = mainBeginIt->end();
        beginRange = begin + beginPos;
        endRange = begin + endPos + 1;
        
        if(endRange < end && beginRange < end)
        {
            mainBeginIt->erase(beginRange, endRange);
            mainBeginIt->insert(beginRange, pasteBeginIt->begin() + 1,
                                pasteBeginIt->end());
        }
        else
        {
            return false;        
        }
        pasteBeginIt++;
    }
    return true;      
}

bool Alignment::reloadAlignment()
{
    if(getNumSeqs() <= 0)
    {
        return false;
    }
    if(userParameters->getOutputOrder() == INPUT)
    {
        return true;
    }
    if((int)outputIndex.size() != getNumSeqs())
    {
        outputIndex.clear();
        return false;
    }
    vector<int> emptyVec;
    string emptyString = "";
    SeqArray outputOrderSeqArray;
    outputOrderSeqArray.resize(getNumSeqs() + 1);
    outputOrderSeqArray[0] = emptyVec;
    vector<string> outputOrderNames;
    outputOrderNames.resize(getNumSeqs() + 1);
    outputOrderNames[0] = emptyString;
    vector<string> outputOrderTitles;
    outputOrderTitles.resize(getNumSeqs() + 1);
    outputOrderTitles[0] = emptyString;
    vector<unsigned long> outputOrderSequenceIds;
    outputOrderSequenceIds.resize(getNumSeqs() + 1);
    outputOrderSequenceIds[0] = 0;
    
    int size = seqArray.size();
    if((seqArray.size() != names.size()) || (seqArray.size() != titles.size()) ||
        sequenceIds.size() != names.size())
    {
        return false;
    }    
    
    int _outIndex;
    // Now for each seq,
    for(int i = 1; i < size; i++)
    { 
        if(i < (int)outputOrderSeqArray.size() && i - 1 < (int)outputIndex.size() &&
           outputIndex[i - 1] < size)
        {
            _outIndex = outputIndex[i - 1];
            outputOrderSeqArray[i] = seqArray[_outIndex];
            outputOrderNames[i] = names[_outIndex];
            outputOrderTitles[i] = titles[_outIndex];
            outputOrderSequenceIds[i] = sequenceIds[_outIndex];
        }
        else
        {
            return false;
        }
    }
    
    // Now we have a copy in the correct order.
    // Remove all the elements from the old ones and set them to be these arrays
    clearSeqArray();
    seqArray = outputOrderSeqArray;
    names.clear();
    names = outputOrderNames;
    titles.clear();
    titles = outputOrderTitles;
    sequenceIds.clear();
    sequenceIds = outputOrderSequenceIds; 

    return true;
}

const vector<int>* Alignment::getSequenceFromUniqueId(unsigned long id)
{
    for(int i = 0; i < (int)sequenceIds.size(); i++)
    {
        if(sequenceIds[i] == id)
        {
            return getSequence(i);
        }
    }
    
    // We have not found it, throw an exception!!!
    throw SequenceNotFoundException();
}

/**
 * Change:
 * Mark 26-1-2007: This function was added to allow access to the unique id.
 */
void Alignment::updateSequence(int index, const vector<int>* seq)
{
    if(index >= 1 && index < (int)seqArray.size())
    {
        seqArray[index] = *seq;
    }
}

/**
 * Change:
 * Mark 17-2-2007: This function was added to check if a res is a gap or not.
 */
bool Alignment::isGap(int seq, int col) const
{
    int res = seqArray[seq][col];
    if(res == gapPos1 || res == gapPos2)
    {
        return true;
    }
    else
    {
        return false;
    }
}
/**
 * This function will be used so that we can create an alignment object from a seqArray.
 * This will be used for the tree iteration.
 */
void Alignment::addSequences(SeqArray* seqVector)
{
    clearAlignment();
    numSeqs = seqVector->size() - 1;
    vector<int> emptyVec;

    seqArray.push_back(emptyVec); // EMPTY sequence
    names.push_back(string(""));
    titles.push_back(string(""));
    sequenceIds.push_back(0);
    cout << "\nThere are " << numSeqs << " in the alignment obj\n";
    for(int i = 1; i <= numSeqs; i++)
    {
        ostringstream name;
        seqArray.push_back((*seqVector)[i]);
        titles.push_back(string(""));
        sequenceIds.push_back(utilityObject->getUniqueSequenceIdentifier());
        name << "name" << numSeqs;
        names.push_back(name.str());        
    }
    
    calculateMaxLengths();
    seqWeight.resize(numSeqs + 1, 100);    
}

}

