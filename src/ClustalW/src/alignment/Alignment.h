/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * The Alignment class is used to store the alignment that is being constructed.
 * It also contains other information such as gap penalty masks etc.
 * An object of this type will be passed by reference to the FileReader. This FileReader
 * and the FileParsers will then set it up properly from the information given in the file.
 * I have decided to put everything into vectors, string etc. No more array*'s, gets rid 
 * of the memory allocation problem.
 *
 * CHANGE: 
 * Mark Jan 16th 2007. I have changed the pasteSequencesIntoPosition function to allow
 * explicit pastes into profile2.
 * Mark 25-1-2007. I have changed the class so that each of the sequences have a unique 
 * identifier. Several functions were changed to allow this.
 * 
 * 16-02-07,Nigel Brown(EMBL): Added friend NameIterator to allow a caller to
 * process the name vector.
 *
 * 23-03-07,Nigel Brown(EMBL): added testUniqueNames() predicate, which
 * compares new sequence names with those in the alignment vector BEFORE
 * appending them. 
 */

 // NOTE NOTE NOTE Very important! The list of sequences begins from 1 to numSeqs.
 // This is because of the fact that the code was written in Fortran where arrays begin at
 // 1. It has become difficult to change this. Ramu has tried before and had problems
 // so we decided to leave it this way.

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <vector>
#include <string>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include "Sequence.h"
#include "../substitutionMatrix/globalmatrix.h"
#include "../general/userparams.h"
#include "../general/VectorOutOfRange.h"
#include "../general/SequenceNotFoundException.h"


// FIXME because this object is used for aligned and unaligned
// sequences it would be nice to have a isAligned flag here (AW)

using namespace std;

namespace clustalw
{

typedef std::vector<vector <int> > SeqArray;

class Alignment
{
    public:
        /* Functions */
        Alignment();
        void addSequences(vector<Sequence>* seqVector);
        void addSequences(SeqArray* seqVector);
        void appendSequences(vector<Sequence>* seqVector);
        vector<Sequence> cutSelectedSequencesFromAlignment(vector<int>* selected);        
        void pasteSequencesIntoPosition(vector<Sequence>* seqVector, int pos, 
                                        bool explicitPasteToProfile2 = false);
                                        
        void resizeSeqArray(int size){seqArray.resize(size); numSeqs = size - 1;
                                      outputIndex.resize(size - 1); names.resize(size);
                                      titles.resize(size);};
        bool addOutputIndex(vector<int>* outputIndexToAdd);
        bool appendOutputIndex(vector<int>* outputIndexToAppend);
        void addSecStructMask1(vector<char>* secStructMaskToAdd);
        void addSecStructMask2(vector<char>* secStructMaskToAdd);
        void addSeqWeight(vector<int>* _seqWeight);
        void addGapPenaltyMask1(vector<char>* gapPenaltyMaskToAdd);
        void addGapPenaltyMask2(vector<char>* gapPenaltyMaskToAdd);
        vector<char>* getSecStructMask1();
        vector<char>* getSecStructMask2();
        const vector<int>* getOutputIndex();
        vector<char>* getGapPenaltyMask1();
        vector<char>* getGapPenaltyMask2();     
        void addSecStructName1(string nameToAdd);
        void addSecStructName2(string nameToAdd);
        int alignScore(void);
        int countGaps(int s1, int s2, int l);
        void resetAlign();
        void fixGaps();
        float countid(int s1, int s2);
        
        const vector<int>* getSequence(int index){return &seqArray[index];}; // For Pairwise!
        const vector<int>* getSequence(int index) const {return &seqArray[index];};
        const vector<int>* getSequenceFromUniqueId(unsigned long id); // For iteration        
        const SeqArray* getSeqArray() const {return &seqArray;}; // For multiple align!
        SeqArray* getSeqArrayForRealloc(){return &seqArray;};
        void updateSequence(int index, const vector<int>* seq);
        
        bool checkAllNamesDifferent(string *offendingSeq);
        bool testUniqueNames(vector<Sequence>* seqVector, string *offendingSeq);
        void clearAlignment();
        void clearSecStruct1();
        void clearSecStruct2();
        void printSequencesAddedInfo();
        
        string getSecStructName1();
        string getSecStructName2();        
        int getNumSeqs() const {return numSeqs;};
        int getMaxNames();
        int getMaxAlnLength(){return maxAlignmentLength;};
        void setMaxAlnLength(int len){maxAlignmentLength = len;};         
        int getLengthLongestSequence();
        int getLengthLongestSequence(int firstSeq, int lastSeq);
        int getSeqLength(int index) const {return seqArray[index].size() - 1;};
        int getSecStructMask1Element(int index);
        int getSecStructMask2Element(int index);
        int getGapPenaltyMask1Element(int index);
        int getGapPenaltyMask2Element(int index);
        int getOutputIndex(int index);
        int getSeqWeight(int index) const;
        const vector<int>* getSeqWeights() const{return &seqWeight;}
        string getName(int index);
        string getTitle(int index);
        vector<int>* QTcalcHistColumnHeights(int firstSeq, int nSeqs, 
                                           Array2D<int>* exceptionalRes); 
                                            // NOTE July 13, for Qt
        
        // NOTE the following functions are to be used when we are doing a profile
        // alignment. It resets the gaps from fixed.
        void resetProfile1();
        void resetProfile2();
        void resetAllSeqWeights();
        
        int searchForString(bool* found, int seq, int beginRes, string search);
        void removeGapsFromSelectedSeqs(vector<int>* selected);
        void removeGapOnlyColsFromSelectedSeqs(vector<int>* selected);
        void removeAllGapOnlyColumns(int fSeq, int lSeq, int profileNum);
        void setDefaultOutputIndex();
        bool removeAllOutsideRange(int beginPos, int endPos);
        bool updateRealignedRange(SeqArray realignedSeqs, int beginPos, int endPos);
        bool reloadAlignment();
        
        int getProfile1NumSeqs(){return profile1NumSeqs;};
        void setProfile1NumSeqs(int value){profile1NumSeqs = value;}
        bool isGap(int seq, int col) const;
        void calculateMaxLengths();

        /**
         * The following functions are for the iteration output order.
         */
        unsigned long getUniqueId(int seq);

        void debugPrintArray(){debugPrintSeqArray(&seqArray);}
        void debugPrintSeqArray(SeqArray* arrayToPrint);
        void debugPrintProfile1();
        void debugPrintProfile2();
        void debugPrintOutAlignInfo();
        void debugPrintAllNames();
        void debugPrintSequences();
        
        /* Attributes */

        /* Friends */
        class NameIterator;
        friend class NameIterator;
        
        class NameIterator 
        {
            private:
                Alignment *alignment;
                vector<string>::iterator i;
            public:
                void begin(Alignment *alignment);
                const string next();
                bool end();
        };
    private:
        /* Functions */
        
        void addSequencesToVector(vector<Sequence>* seqVector);
        int getSequenceLength(int index);
        void sortScores(vector<float>* scores, int f, int l);
        void swap(vector<float>* scores, int s1, int s2);
        bool keepPortionOfSeqArray(int beginRangeIndex, int endRangeIndex);
        
        void clearSeqArray();
        /* Attributes */
        int maxNames;
        int maxAlignmentLength;
        int lengthLongestSequence;
        int numSeqs;
        vector<int> outputIndex;
        vector<unsigned long> sequenceIds; // Mark change: To help with output order 
        vector<int> seqWeight;
        SeqArray seqArray;
        vector<string> names;
        vector<string> titles;
        vector<char> gapPenaltyMask1;
        vector<char> gapPenaltyMask2;
        vector<char> secStructMask1;
        vector<char> secStructMask2;
        string secStructName1;
        string secStructName2;
        vector<int> histogramColumnHeights; // NOTE July 13, for Qt
        int profile1NumSeqs;
        int gapPos1, gapPos2;
};
}
#endif

