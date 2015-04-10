/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "DebugLog.h"
#include <sstream>
#include <iostream>
namespace clustalw
{

DebugLog::DebugLog(std::string _logFileName)
 : logFileName(_logFileName),
   logFile(0),
   numScores(0),
   sumSoFar(0.0),
   averageScore(0.0),
   minScore(0.0),
   maxScore(0.0)
{
    logFile = new std::ofstream();  
    logFile->open(logFileName.c_str(), ios::out);
    
    if(logFile->is_open())
    {
        std::cout << "Logging debug info to file: " << logFileName << std::endl;
    }
    else
    {
        std::cerr << "Could not open log file.\n";
    }    
}

DebugLog::~DebugLog()
{
    // Release the file!
    logFile->close();
    delete logFile;
}

void DebugLog::logMsg(std::string msg)
{
    if(logFile->is_open())
    {
        (*logFile) << msg << "\n";
    }
}

void DebugLog::logScore(float x)
{
    if(x < minScore)
    {
        minScore = x;
    }
    if(x > maxScore)
    {
        maxScore = x;
    }
        
    sumSoFar += x;
    numScores++;
}

void DebugLog::printScoreInfo()
{
    if(numScores > 0)
    {
        averageScore = sumSoFar / static_cast<float>(numScores);
        ostringstream outs;
        outs << "SCORE INFO--------------------------------------------------->"
             << " The score was calculated " << numScores << " times. The average = "
             << averageScore << "\n" << "The max score=" << maxScore << " The min score="
             << minScore << "\n";
        logMsg(outs.str());
    }
}

}
