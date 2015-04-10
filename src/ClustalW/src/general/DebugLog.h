/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * This class is used to log messages to a file that is specified when the object 
 * is created.
 * The file is closed when the object is destroyed. The user simply needs to 
 * create a DebugLog
 * object, and then call the logMsg function whenever they wish to write something to the
 * file. 
 */
 
#ifndef DEBUGLOG_H
#define DEBUGLOG_H

#include <string>
#include <fstream>

namespace clustalw
{

using namespace std;

class DebugLog
{
    public:
        DebugLog(std::string);
        ~DebugLog();
        void logMsg(std::string);
        void logScore(float x);
        void printScoreInfo();
    private:
        /* Attributes */
        std::string logFileName;
        std::ofstream* logFile;
        int numScores;
        float sumSoFar;
        float averageScore;
        float minScore;
        float maxScore;
        /* Functions */
        DebugLog();
        
};

}
#endif
