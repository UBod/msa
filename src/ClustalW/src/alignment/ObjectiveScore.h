/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * This ObjectiveScore class is used to provide an objective function to score 
 * an alignment.
 * It is used with iteration to improve an alignment.
 */
 
#ifndef OBJECTIVESCORE_H
#define OBJECTIVESCORE_H
#include "../general/clustalw.h"
#include "../substitutionMatrix/globalmatrix.h" 
namespace clustalw
{

class Alignment;
typedef struct
{
    int first;
    int second; 
} Pair;  

class ObjectiveScore
{   
    public:
        ObjectiveScore();
        long getScore(const Alignment* alnToScore);    
    private:
        
        float scoreLetters(int seq1, int seq2);
        float scoreGaps(int seq1, int seq2);
        void calcNormalisedSeqWeights(const vector<int>* seqWeight, 
                                      vector<float>* normSeqWeight);
        long score;
        int matrix[NUMRES][NUMRES];
        const Alignment* alignToScore;
        long scale;
        int weightScale;
        int sagaGapEx, sagaGapOp;
        int gapPos1, gapPos2;
        static const int BOTHGAPS = 0;
        static const int NOGAPS = -1;
        static const int GAPINSEQB = 1;
        static const int GAPINSEQA = 2;
};

}
#endif
