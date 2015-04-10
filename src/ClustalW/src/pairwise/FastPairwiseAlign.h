/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef FASTPAIRWISEALIGN_H
#define FASTPAIRWISEALIGN_H

#include "PairwiseAlignBase.h"

namespace clustalw
{

class FastPairwiseAlign : public PairwiseAlignBase
{
    public:
        /* Functions */
        FastPairwiseAlign();
	virtual ~FastPairwiseAlign(){};

        virtual void pairwiseAlign(Alignment *alignPtr, DistMatrix *distMat, int iStart, 
                                   int iEnd, int jStart, int jEnd);
        /* Attributes */

    private:
        /* Functions */
        void pairAlign(const vector<int>* seq, int l1, int l2); 
        void makePPtrs(vector<int>& tptr, vector<int>& pl, const vector<int>* seq, 
                       int length); 
        void makeNPtrs(vector<int>& tptr, vector<int>& pl, const vector<int>* seq, 
                       int length);
        void putFrag(int fs, int v1, int v2, int flen);
        int fragRelPos(int a1, int b1, int a2, int b2);
        void desQuickSort(vector<int>& array1, vector<int>& array2, int arraySize);

        /* Attributes */
        
        vector<int> displ;
        vector<int> zza;
        vector<int> zzb;
        vector<int> zzc;
        vector<int> zzd;        
        int next;
        int currFrag;
        int maxSoFar;
        int vatend;
        
        Array2D<int> accum;
        vector<int> diagIndex;
        vector<int> slopes;
        int _maxAlnLength;        
};
}
#endif
