/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef MYERSMILLERPROFILEALIGN_H
#define MYERSMILLERPROFILEALIGN_H

#include <vector>
#include "ProfileAlignAlgorithm.h"
#include "ProfileStandard.h"
#include "ProfileWithSub.h"
namespace clustalw
{

class MyersMillerProfileAlign : public ProfileAlignAlgorithm
{
    public:
  virtual ~MyersMillerProfileAlign(){};

    /* Functions */
        MyersMillerProfileAlign();
        virtual int profileAlign(Alignment* alnPtr, DistMatrix* distMat, 
                                 vector<int>* group, int* aligned);
    /* Attributes */
    
    private:
    /* Functions */
        void addGGaps(Alignment* alnPtr, SeqArray* seqArray);
        void addGGapsMask(vector<char>* mask,int len, vector<int>* path1, vector<int>* path2);
        int prfScore(int n, int m);
        int progTracepath();
        void progDel(int k);
        void progAdd(int k);
        void progAlign();
        int progDiff(int A, int B, int M, int N, int go1, int go2);
        int openPenalty1(int i, int j);
        int extPenalty1(int i, int j);
        int gapPenalty1(int i, int j, int k);
        int openPenalty2(int i, int j);
        int extPenalty2(int i, int j);
        int gapPenalty2(int i, int j, int k);    
    /* Attributes */
        ProfileWithSub* profileWithSub;
        ProfileStandard* profileStandard;
        int gapcoef1;
        int gapcoef2;
        int lencoef1;
        int lencoef2;
        vector<int> displ;
        vector<int> gS;
        vector<int> HH;
        vector<int> DD;
        vector<int> RR;
        vector<int> SS;
        vector<int> alnPath1;
        vector<int> alnPath2;
        int printPtr;
        int lastPrint;                        
        int matrix[32][32];
        vector<int> gaps;
        bool switchProfiles;
        const SeqArray* profile1;
        const SeqArray* profile2;
        int _gapPos1, _gapPos2;
        int alignmentLength;                   
};

}
#endif
