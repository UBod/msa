/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * The class CommandLineParser is used to parse the command line arguments. It then
 * sets some parameters, and calls the required functions.
 * To get it to parse the command line, and execute the required job, create a 
 * CommandLineParser object and pass the list of arguments to the constructor.
 */

#ifndef COMMANDLINEPARSER_H
#define COMMANDLINEPARSER_H

#include <vector>
#include <string>
#include "../Clustal.h"
#include "../general/clustalw.h"
#include "../general/userparams.h"
#include "../general/utils.h"
#include "../general/debuglogObject.h"
#include "../general/statsObject.h"

#include <Rcpp.h>
namespace clustalw
{

typedef std::vector<std::string> StringArray;

typedef struct {
    const char *str;
    int *flag;
    int type;
    StringArray* arg;
} CmdLineData;

class CommandLineParser
{
    public:
        /* Functions */
        CommandLineParser(StringArray* args, bool xmenus);
        ~CommandLineParser();

        void run(StringArray* args, bool xmenus,
        				   ClustalWInput *input, ClustalWOutput *output);

    private:
        /* Functions */
        
        /** separate argument options and values
         * 
         * put all the parameters and their values in the
         * vectors, and also to check that all parameters that require a value
         * have it.
         *
         * @param args contains parameters as value and optional option pairs
         * @param parameters used for storing parsed arguments
         * @param paramArgs used for storing values of parsed arguments
         * @return number of parsed parameters or -1 on error
         */
        int checkParam(StringArray* args, StringArray* params, StringArray* paramArg);
        void setOptionalParam(Rcpp::NumericMatrix substitutionMatrix);
        int findMatch(string probe, StringArray* list, int n);
        CmdLineData getCmdLineDataStruct(const char *str, int *flag, int type, StringArray* arg);
        void printCmdLineData(const CmdLineData& temp);
        string ConvertStringToLower(string strToConvert);
        void exitWithErrorMsg(string msg);
        void reportBadOptionAndExit(string option, string expectedType);
        void reportInvalidOptionAndExit(string option);

        
        /* Attributes */
        Clustal *clustalObj;
        
        static const int MAXARGS = 100;
        static const int NOARG = 0;
        static const int INTARG = 1;
        static const int FLTARG = 2;
        static const int STRARG = 3;
        static const int FILARG = 4;
        static const int OPTARG = 5;
        
        static const char *DUMMY_R_MATRIX;

        int setOptions;
        int setHelp;
        int setFullHelp;
        int setQuiet;
        int setInteractive;
        int setBatch;
        int setGapOpen;
        int setGapExtend;
        int setPWGapOpen;
        int setPWGapExtend;
        int setOutOrder;
        int setBootLabels;
        int setPWMatrix;
        int setMatrix;
        int setPWDNAMatrix;
        int setDNAMatrix;
        int setNegative;
        int setNoWeights;
        int setOutput;
        int setOutputTree;
        int setQuickTree;
        int setType;
        int setCase;
        int setSeqNo;
        int setSeqNoRange;
        int setRange;
        int setTransWeight;
        int setSeed;
        int setScore;
        int setWindow;
        int setKtuple;
        int setKimura;
        int setTopDiags;
        int setPairGap;
        int setTossGaps;
        int setNoPGap;
        int setNoHGap;
        int setNoVGap;
        int setHGapRes;
        int setUseEndGaps;
        int setMaxDiv;
        int setGapDist;
        int setDebug;
        int setOutfile;
        int setInfile;
        int setProfile1;
        int setProfile2;
        int setAlign;
        int setConvert;
        int setNewTree;
        int setUseTree;
        int setNewTree1;
        int setUseTree1;
        int setNewTree2;
        int setUseTree2;
        int setBootstrap;
        int setTree;
        int setProfile;
        int setSequences;
        int setSecStruct1;
        int setSecStruct2;
        int setSecStructOutput;
        int setHelixGap;
        int setStrandGap;
        int setLoopGap;
        int setTerminalGap;
        int setHelixEndIn;
        int setHelixEndOut;
        int setStrandEndIn;
        int setStrandEndOut;
        int profileType;
        int setDoIteration;
        int setNumIterations;
        int setTreeAlgorithm;
        int setMaxSeqLen;
        int setStatsFile;
        int setOutputPim;
        
        string userMatrixName;
        string pwUserMatrixName;
        string DNAUserMatrixName;
        string pwDNAUserMatrixName;
        
        CmdLineData cmdLineFile[4];
        CmdLineData cmdLineVerb[20];
        CmdLineData cmdLinePara[56];
        
        string clustalTreeName;
        string distTreeName;
        string phylipTreeName;
        string nexusTreeName;
        string p1TreeName;
        string p2TreeName;
        string pimName;
        
        StringArray* typeArg;
        StringArray* bootLabelsArg;
        StringArray* outOrderArg;
        StringArray* caseArg;
        StringArray* seqNoArg;
        StringArray* seqNoRangeArg;
        StringArray* scoreArg;
        StringArray* outputArg;
        StringArray* outputTreeArg;
        StringArray* outputSecStrArg;
        StringArray* cmdLineType;
        StringArray* clusterAlgorithm;
        StringArray* iterationArg;
        
        StringArray* params; // parameter names
        StringArray* paramArg; // parameter values
};

}
#endif
