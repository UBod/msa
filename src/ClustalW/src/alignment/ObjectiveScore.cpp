/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <sstream>
#include "ObjectiveScore.h"
#include "Alignment.h"
#include "../general/userparams.h"
#include "../general/debuglogObject.h"

namespace clustalw
{

ObjectiveScore::ObjectiveScore()
 : score(0),
   alignToScore(0),
   scale(100000),
   sagaGapEx(12),
   sagaGapOp(8),
   gapPos1(userParameters->getGapPos1()),
   gapPos2(userParameters->getGapPos2())
{

}

long ObjectiveScore::getScore(const Alignment* alnToScore)
{
    #if DEBUGFULL
        if(logObject && DEBUGLOG)
        {    
            logObject->logMsg("In getScore function");
        }
    #endif
        
    alignToScore = alnToScore;
    // If it doesnt point to any object return 0;
    if(!alignToScore)
    {
        return 0;
    }
    
    int maxRes = subMatrix->getAlnScoreMatrix(matrix);
    if (maxRes == 0)
    {
        utilityObject->error("Matrix for alignment scoring not found\n");
        return 0;
    }
    
    const vector<int>* seqWeight = alignToScore->getSeqWeights();
    vector<float> normalisedSeqWeights;
    
    calcNormalisedSeqWeights(seqWeight, &normalisedSeqWeights);
        
    int seq1, seq2;

    float w1, w2 = 1.0;
    float weight = 1.0;
    score = 0;
    float pwScoreLetters = 0, pwScoreGaps = 0, pwScore = 0;
    float scoreTotal = 0.0;
    int numSeqs = alignToScore->getNumSeqs();
    int sizeNormalSeqWeight = normalisedSeqWeights.size();
    
    for (seq1 = 1; seq1 <= numSeqs && seq1 <= sizeNormalSeqWeight; seq1++) 
    {
        w1 = normalisedSeqWeights[seq1 - 1];
        
        for (seq2 = seq1 + 1; seq2 <= numSeqs && seq2 <= sizeNormalSeqWeight; seq2++)
        {
            w2 = normalisedSeqWeights[seq2 - 1];
            weight = w1 * w2;
            
            pwScoreLetters = scoreLetters(seq1, seq2);
            pwScoreGaps = scoreGaps(seq1, seq2);
            pwScore = pwScoreLetters + pwScoreGaps;
                       
            scoreTotal += weight * pwScore;
            #if DEBUGFULL
                if(logObject && DEBUGLOG)
                {    
                    ostringstream outs;
                    outs << " weight = " << weight 
                         << " scoreLetters = " << pwScoreLetters << " scoreGaps = " 
                         << pwScoreGaps << " scorePair = " 
                         << pwScore << "\n"
                         << "scoreTotal = " << scoreTotal << "\n";
                    logObject->logMsg(outs.str());
                }
            #endif             
        }
    }

    score = static_cast<long>(scoreTotal);
    #if DEBUGFULL
        if(logObject && DEBUGLOG)
        {    
            ostringstream outs;
            outs << " score = " << score;
            logObject->logMsg(outs.str());
        }
    #endif
    
    utilityObject->info("Alignment Score %d\n", score);    
    return score;
}

float ObjectiveScore::scoreLetters(int seq1, int seq2)
{
    if(!alignToScore)
    {
        return 0;
    }
    const unsigned numColsSeq1 = alignToScore->getSeqLength(seq1);
    const unsigned numColsSeq2 = alignToScore->getSeqLength(seq2);
    
    if(numColsSeq1 != numColsSeq2)
    {
        return 0; // The sequences should be the same length after alignment.
    }
    
    float scoreLetters = 0;
    unsigned colStart = 1;
    bool gap1, gap2;
    
    for(unsigned col = 1; col < numColsSeq1; col++)
    {
        gap1 = alignToScore->isGap(seq1, col);
        gap2 = alignToScore->isGap(seq2, col);
        
        if(!gap1 || !gap2)
        {
            colStart = col;
            break;
        }
    }
 
    unsigned colEnd = numColsSeq1;
    
    for(unsigned col = numColsSeq1; col >= 1; --col)
    {
        gap1 = alignToScore->isGap(seq1, col);
        gap2 = alignToScore->isGap(seq2, col);
        
        if(!gap1 || !gap2)
        {
            colEnd = col;
            break;
        }
    }
    
    const SeqArray* seqArray = alignToScore->getSeqArray();
    int scoreMatch = 0;

        
    for(unsigned col = colStart; col <= colEnd; col++)
    {
        int res1 = (*seqArray)[seq1][col];
        int res2 = (*seqArray)[seq2][col];
        scoreMatch = matrix[res1][res2];
        scoreLetters += scoreMatch;
    }   
    
    return scoreLetters;
}

float ObjectiveScore::scoreGaps(int seq1, int seq2)
{
    if(!alignToScore)
    {
        return 0;
    }
    const unsigned numColsSeq1 = alignToScore->getSeqLength(seq1);
    const unsigned numColsSeq2 = alignToScore->getSeqLength(seq2);
    
    if(numColsSeq1 != numColsSeq2)
    {
        return 0; // The sequences should be the same length after alignment.
    }
    
    unsigned colStart = 1;
    bool gap1, gap2;
    
    for(unsigned col = 1; col < numColsSeq1; col++)
    {
        gap1 = alignToScore->isGap(seq1, col);
        gap2 = alignToScore->isGap(seq2, col);
        
        if(!gap1 || !gap2)
        {
            colStart = col;
            break;
        }
    }
 
    unsigned colEnd = numColsSeq1;
    
    for(unsigned col = numColsSeq1; col >= 1; --col)
    {
        gap1 = alignToScore->isGap(seq1, col);
        gap2 = alignToScore->isGap(seq2, col);
        
        if(!gap1 || !gap2)
        {
            colEnd = col;
            break;
        }
    }
 
    bool inGap1 = false;
    bool inGap2 = false;
    float gapOpen = userParameters->getGapOpen();
    float gapExtend = userParameters->getGapExtend();
    
    float scoreGaps = 0;
    for(unsigned col = colStart; col <= colEnd; col++)
    {
        gap1 = alignToScore->isGap(seq1, col);
        gap2 = alignToScore->isGap(seq2, col);
        
        if(gap1 && gap2)
        {
            continue;
        }
        if(gap1)
        {
            if(!inGap1)
            {
                // NOTE I left out the option of having different end gaps stuff.
                scoreGaps += gapOpen;
                inGap1 = true; // Opening a gap in seq1
            }
            else
            {
                // Already in a gap
                scoreGaps += gapExtend;
            }
            continue;
        }
        else if(gap2)
        {
            if(!inGap2)
            {
                // NOTE I left out the option of having different end gaps stuff.
                scoreGaps += gapOpen;
                inGap2 = true; // Opening a gap in seq2
            }
            else
            {
                // Already in a gap
                scoreGaps += gapExtend;
            }
            continue;        
        }
        inGap1 = inGap2 = false;    
    }       
    return scoreGaps;
} 

void ObjectiveScore::calcNormalisedSeqWeights(const vector<int>* seqWeight, 
                                              vector<float>* normSeqWeight)
{
    if(!seqWeight || !normSeqWeight)
    {
        return;
    }
    
    int sumWeights = 0;
    for(int i = 0; i < (int)seqWeight->size() - 1; i++)
    {       
        sumWeights += (*seqWeight)[i];
    }
    
    normSeqWeight->resize(seqWeight->size());
    for(int i = 0; i < (int)seqWeight->size() - 1; i++)
    {
        (*normSeqWeight)[i] = static_cast<float>((*seqWeight)[i]) /
                              static_cast<float>(sumWeights);
    #if DEBUGFULL
        if(logObject && DEBUGLOG)
        {    
            ostringstream outs;
            outs << " normSeqWeight[i] = " << (*normSeqWeight)[i];
            logObject->logMsg(outs.str());
        }
    #endif                           
    }
  
}                                      
        
}

