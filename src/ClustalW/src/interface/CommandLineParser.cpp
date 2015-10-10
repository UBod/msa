/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <cctype>
#include <cstdio>
#include <iostream>
#include <fstream>
#include "CommandLineParser.h"
#include "../substitutionMatrix/globalmatrix.h"
#include "general/Utility.h"
#include "general/statsObject.h"

using namespace std;

namespace clustalw
{

const char *CommandLineParser::DUMMY_R_MATRIX ="dummyR.matrix";

CommandLineParser::CommandLineParser(StringArray* args, bool xmenus)
:   setOptions(-1),
    setHelp(-1),
    setFullHelp(-1),
    setQuiet(-1),
    setInteractive(-1),
    setBatch(-1),
    setGapOpen(-1),
    setGapExtend(-1),
    setPWGapOpen(-1),
    setPWGapExtend(-1),
    setOutOrder(-1),
    setBootLabels(-1),
    setPWMatrix(-1),
    setMatrix(-1),
    setPWDNAMatrix(-1),
    setDNAMatrix(-1),
    setNegative(-1),
    setNoWeights(-1),
    setOutput(-1),
    setOutputTree(-1),
    setQuickTree(-1),
    setType(-1),
    setCase(-1),
    setSeqNo(-1),
    setSeqNoRange(-1),
    setRange(-1),
    setTransWeight(-1),
    setSeed(-1),
    setScore(-1),
    setWindow(-1),
    setKtuple(-1),
    setKimura(-1),
    setTopDiags(-1),
    setPairGap(-1),
    setTossGaps(-1),
    setNoPGap(-1),
    setNoHGap(-1),
    setNoVGap(-1),
    setHGapRes(-1),
    setUseEndGaps(-1),
    setMaxDiv(-1),
    setGapDist(-1),
    setDebug(-1),
    setOutfile(-1),
    setInfile(-1),
    setProfile1(-1),
    setProfile2(-1),
    setAlign(-1),
    setConvert(-1),
    setNewTree(-1),
    setUseTree(-1),
    setNewTree1(-1),
    setUseTree1(-1),
    setNewTree2(-1),
    setUseTree2(-1),
    setBootstrap(-1),
    setTree(-1),
    setProfile(-1),
    setSequences(-1),
    setSecStruct1(-1),
    setSecStruct2(-1),
    setSecStructOutput(-1),
    setHelixGap(-1),
    setStrandGap(-1),
    setLoopGap(-1),
    setTerminalGap(-1),
    setHelixEndIn(-1),
    setHelixEndOut(-1),
    setStrandEndIn(-1),
    setStrandEndOut(-1),
    profileType(PROFILE),
    setDoIteration(-1),
    setNumIterations(-1),
    setTreeAlgorithm(-1),
    setMaxSeqLen(-1),
    setStatsFile(-1),
    setOutputPim(-1)
{
    int ctr=0;
    
    // The rest of the variables are arrays!
    try
    {
        clustalObj = new Clustal();
        
        // selecting the size prevents the resizing of the vector which is expensive. 
        typeArg = new StringArray(3);
        bootLabelsArg = new StringArray(3);
        outOrderArg = new StringArray(3);
        caseArg = new StringArray(3);
        seqNoArg = new StringArray(3);
        seqNoRangeArg = new StringArray(3);
        scoreArg = new StringArray(3);
        outputArg = new StringArray(8);
        outputTreeArg = new StringArray(5);
        outputSecStrArg = new StringArray(5);
        cmdLineType = new StringArray(6);
        clusterAlgorithm = new StringArray(3);
        iterationArg = new StringArray(4);
        
        params = new StringArray; // Wait until I need it!!!!!!!!!
        paramArg = new StringArray;
    }
    catch(const exception &ex)
    {
        cerr << ex.what() << endl;
        cerr << "Terminating program. Cannot continue" << std::endl;
        throw 1;
    }
                           
    (*typeArg)[0] = "protein";
    (*typeArg)[1] = "dna";
    (*typeArg)[2] = "";
        
    (*bootLabelsArg)[0] = "node";
    (*bootLabelsArg)[1] = "branch";
    (*bootLabelsArg)[2] = "";
    
    (*outOrderArg)[0] = "input";
    (*outOrderArg)[1] = "aligned";
    (*outOrderArg)[2] = "";
    
    (*caseArg)[0] = "lower";
    (*caseArg)[1] = "upper";
    (*caseArg)[2] = "";
    
    (*seqNoArg)[0] = "off";
    (*seqNoArg)[1] = "on";
    (*seqNoArg)[2] = "";
    
    (*seqNoRangeArg)[0] = "off";
    (*seqNoRangeArg)[1] = "on";
    (*seqNoRangeArg)[2] = "";
        
    (*scoreArg)[0] = "percent";
    (*scoreArg)[1] = "absolute";
    (*scoreArg)[2] = "";
    
    (*outputArg)[0] = "gcg";
    (*outputArg)[1] = "gde";
    (*outputArg)[2] = "pir";
    (*outputArg)[3] = "phylip";
    (*outputArg)[4] = "nexus";
    (*outputArg)[5] = "fasta";
    (*outputArg)[6] = "clustal";
    (*outputArg)[7] = "";
    
    (*outputTreeArg)[0] = "nj";
    (*outputTreeArg)[1] = "phylip";
    (*outputTreeArg)[2] = "dist";
    (*outputTreeArg)[3] = "nexus";
    (*outputTreeArg)[4] = "";
    
    (*outputSecStrArg)[0] = "structure";
    (*outputSecStrArg)[1] = "mask";
    (*outputSecStrArg)[2] = "both";
    (*outputSecStrArg)[3] = "none";
    (*outputSecStrArg)[4] = "";
    
    (*cmdLineType)[0] = " ";
    (*cmdLineType)[1] = "=n ";
    (*cmdLineType)[2] = "=f ";
    (*cmdLineType)[3] = "=string ";
    (*cmdLineType)[4] = "=filename ";
    (*cmdLineType)[5] = "";
    
    (*clusterAlgorithm)[0] = "nj";
    (*clusterAlgorithm)[1] = "upgma";
    (*clusterAlgorithm)[2] = "";

    (*iterationArg)[0] = "tree";
    (*iterationArg)[1] = "alignment";
    (*iterationArg)[2] = "none";
    (*iterationArg)[3] = "";
        
    userMatrixName = "";
    pwUserMatrixName = "";
    DNAUserMatrixName = "";
    pwDNAUserMatrixName = "";
    
    clustalTreeName = "";
    distTreeName = "";
    phylipTreeName = "";
    nexusTreeName = "";
    p1TreeName = "";
    p2TreeName = "";
    pimName = "";

    // NOTE there were only 3 params for the last one, so I put in NULL for the 4th.
    ctr=0;
    cmdLineFile[ctr++] = getCmdLineDataStruct("infile", &setInfile, FILARG, NULL);
    cmdLineFile[ctr++] = getCmdLineDataStruct("profile1", &setProfile1, FILARG, NULL);
    cmdLineFile[ctr++] = getCmdLineDataStruct("profile2", &setProfile2, FILARG, NULL);
    cmdLineFile[ctr++] = getCmdLineDataStruct("", NULL, -1, NULL);
    // FIXME: final ctr index is hardcoded in CommandLineParser
    
    ctr=0;
    cmdLineVerb[ctr++] = getCmdLineDataStruct("help", &setHelp, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("fullhelp", &setFullHelp, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("quiet", &setQuiet, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("check", &setHelp, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("options", &setOptions, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("align", &setAlign, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("newtree", &setNewTree, FILARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("usetree", &setUseTree, FILARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("newtree1", &setNewTree1, FILARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("usetree1", &setUseTree1, FILARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("newtree2", &setNewTree2, FILARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("usetree2", &setUseTree2, FILARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("bootstrap", &setBootstrap, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("tree", &setTree, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("quicktree", &setQuickTree, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("convert", &setConvert, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("interactive", &setInteractive, NOARG, NULL);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("batch", &setBatch, NOARG, NULL);
    // Mark change 16-feb-2007 I added options for doing LE and iteration
    cmdLineVerb[ctr++] = getCmdLineDataStruct("iteration", &setDoIteration, 
                                           OPTARG, iterationArg);
    cmdLineVerb[ctr++] = getCmdLineDataStruct("", NULL, -1, NULL);
    // FIXME: final ctr index is hardcoded in CommandLineParser.h
    
    // NOTE Start back here!!!!!!!!!!!!
    ctr=0;
    cmdLinePara[ctr++] = getCmdLineDataStruct("type", &setType, OPTARG, typeArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("profile", &setProfile, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("sequences", &setSequences, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("matrix", &setMatrix, FILARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("dnamatrix", &setDNAMatrix, FILARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("negative", &setNegative, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("noweights", &setNoWeights, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("gapopen", &setGapOpen, FLTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("gapext", &setGapExtend, FLTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("endgaps", &setUseEndGaps, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("nopgap", &setNoPGap, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("nohgap", &setNoHGap, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("novgap", &setNoVGap, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("hgapresidues", &setHGapRes, STRARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("maxdiv", &setMaxDiv, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("gapdist", &setGapDist, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("pwmatrix", &setPWMatrix, FILARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("pwdnamatrix", &setPWDNAMatrix, FILARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("pwgapopen", &setPWGapOpen, FLTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("pwgapext", &setPWGapExtend, FLTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("ktuple", &setKtuple, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("window", &setWindow, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("pairgap", &setPairGap, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("topdiags", &setTopDiags, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("score", &setScore, OPTARG, scoreArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("transweight", &setTransWeight, FLTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("seed", &setSeed, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("kimura", &setKimura, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("tossgaps", &setTossGaps, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("bootlabels", &setBootLabels, OPTARG,
                                             bootLabelsArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("debug", &setDebug, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("output", &setOutput, OPTARG, outputArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("outputtree", &setOutputTree, OPTARG,
                                             outputTreeArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("outfile", &setOutfile, FILARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("outorder", &setOutOrder, OPTARG, outOrderArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("case", &setCase, OPTARG, caseArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("seqnos", &setSeqNo, OPTARG, seqNoArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("seqno_range", &setSeqNoRange, OPTARG,
                                             seqNoRangeArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("range", &setRange, STRARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("nosecstr1", &setSecStruct1, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("nosecstr2", &setSecStruct2, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("secstrout", &setSecStructOutput, OPTARG,
                                             outputSecStrArg);
    cmdLinePara[ctr++] = getCmdLineDataStruct("helixgap", &setHelixGap, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("strandgap", &setStrandGap, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("loopgap", &setLoopGap, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("terminalgap", &setTerminalGap, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("helixendin", &setHelixEndIn, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("helixendout", &setHelixEndOut, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("strandendin", &setStrandEndIn, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("strandendout",&setStrandEndOut, INTARG, NULL);
    // NOTE these one was added to test the new LE scoring and iterations
    cmdLinePara[ctr++] = getCmdLineDataStruct("numiter",&setNumIterations, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("clustering", &setTreeAlgorithm, OPTARG,
                                           clusterAlgorithm);
    cmdLinePara[ctr++] = getCmdLineDataStruct("maxseqlen", &setMaxSeqLen, INTARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("stats", &setStatsFile, FILARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("pim", &setOutputPim, NOARG, NULL);
    cmdLinePara[ctr++] = getCmdLineDataStruct("", NULL, -1, NULL);
    // FIXME: final ctr index is hardcoded in CommandLineParser

    //parseParams(args, xmenus);
}

CommandLineParser::~CommandLineParser()
{
    // Free up menory used here!
    // NOTE any dynamically allocated memory (new) must be deleted.
    delete clustalObj;
    delete typeArg;
    delete bootLabelsArg;
    delete outOrderArg;
    delete caseArg;
    delete seqNoArg;
    delete seqNoRangeArg;
    delete scoreArg;
    delete outputArg;
    delete outputTreeArg;
    delete outputSecStrArg;
    delete cmdLineType;        
    delete params;
    delete paramArg;   
    delete clusterAlgorithm;
    delete iterationArg;
}

void CommandLineParser::run(StringArray* args, bool xmenus, ClustalWInput *input, ClustalWOutput *output)
{
	#if DEBUGFULL
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("Parsing Command line Parameters!\n");
        }
    #endif
        
    int i, j, temp;
    //int len;
    //static int cl_error_code = 0;
    //char path[FILENAMELEN];
    int numparams = 0;

    bool doAlign, doConvert, doAlignUseOldTree, doGuideTreeOnly, doTreeFromAlign, 
         doBootstrap, doProfileAlign, doSomething;

    /*if (!xmenus && userParameters->getDisplayInfo())
    {
        cout <<  std::endl << std::endl << std::endl;
        cout << " CLUSTAL " << userParameters->getRevisionLevel() 
             << " Multiple Sequence Alignments" << std::endl << std::endl << std::endl;
    }*/

    doAlign = doConvert = doAlignUseOldTree = doGuideTreeOnly = doTreeFromAlign = false;
    doBootstrap = doProfileAlign = doSomething = false; 
    
    numparams = checkParam(args, params, paramArg);
    
    if (numparams < 0) 
    {
    	throw 1;
    }

//**********************************************************************
//*** Note: This part of the code is to print out the options with  ****
//*** their expected value types/ranges                             ****
//**********************************************************************
    
    if(setHelp != -1) 
    {
        userParameters->setHelpFlag(true);
        if (xmenus) {
            // gui will display help
            // but parse rest of args anyway and don't return
        } else {
            clustalObj->getHelp('9');
            throw 1;
        }
    }

    if(setFullHelp != -1) 
    {
        userParameters->setFullHelpFlag(true);
        if (xmenus) {
            // gui will handle this
            // but parse rest of args anyway and don't return
        } else {
            clustalObj->getFullHelp();
            throw 1;
        }
    }

    if(setQuiet != -1) 
    {
        userParameters->setDisplayInfo(false);
        utilityObject->beQuiet(true);
    } else {
        userParameters->setDisplayInfo(true);
        utilityObject->beQuiet(false);
    }

        
    // need to check maxseqlen before reading input file
    if (setMaxSeqLen != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting max allowed sequence length.");
            }    
        #endif      
        temp = 0;
        if((*paramArg)[setMaxSeqLen].length() > 0)
        {
            if (sscanf((*paramArg)[setMaxSeqLen].c_str(),"%d", &temp) != 1) 
            {
                reportBadOptionAndExit("maxseqlen", "integer");
            }
        }
        if(temp > 0)
        {
            userParameters->setMaxAllowedSeqLength(temp);
        }
        else
        {
            cerr << "Cannot use a negative value for maximum sequence length. Using default"  << std::endl;
        }
        
    }

    if (setStatsFile != -1)
    {
        if((*paramArg)[setStatsFile].length() > 0) 
        {
            statsObject->setEnabled(true);
            statsObject->setStatsFile((*paramArg)[setStatsFile]);
        }
    }


    if (setOutputPim != -1)
    {
            userParameters->setOutputPim(true);
    }


    /*if(setDoIteration != -1)
    {
        userParameters->setDoRemoveFirstIteration(true);
    }*/
    
    if(setDoIteration != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting iteration parameter.");
            }
        #endif
        if((*paramArg)[setDoIteration].length() > 0) 
        { 
            temp = findMatch((*paramArg)[setDoIteration], iterationArg, 3);
            if(temp == 0)
            {
                userParameters->setDoRemoveFirstIteration(TREE);
            }
            else if(temp == 1)
            {
                userParameters->setDoRemoveFirstIteration(ALIGNMENT);
            }
            else if(temp == 2)
            {               
                userParameters->setDoRemoveFirstIteration(NONE);
            }
            else
            {
                cerr << "Unknown option for iteration. Setting to NONE"  << std::endl;
                userParameters->setDoRemoveFirstIteration(NONE);
            }
        }
    }    
    if(setOptions != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("Displaying options!\n");
            }
        #endif
                
        cout << "clustalw option list:-" << std::endl;
        
        for (i = 0; cmdLineVerb[i].str[0] != '\0'; i++) 
        {
            cout << "\t\t" << default_commandsep << cmdLineVerb[i].str 
                 << (*cmdLineType)[cmdLineVerb[i].type];

            if (cmdLineVerb[i].type == OPTARG) 
            {
                if (cmdLineVerb[i].arg !=  NULL)
                {     
                    cout << "=" << cmdLineVerb[i].arg->at(0);
                    for(j = 1; j < (int)cmdLineVerb[i].arg->size() - 1; j++)
                    {
                        cout << " OR " << cmdLineVerb[i].arg->at(j);
                    }
                }
            }
            cout << std::endl;
        }
        for (i = 0; cmdLineFile[i].str[0] != '\0'; i++) 
        {
            cout<<"\t\t" << default_commandsep << cmdLineFile[i].str 
                << (*cmdLineType)[cmdLineFile[i].type];
                
            if (cmdLineFile[i].type == OPTARG) 
            {
                if (cmdLineFile[i].arg !=  NULL)
                {
                    
                    cout << "=" << cmdLineFile[i].arg->at(0);
                    for(j = 1; j < (int)cmdLineFile[i].arg->size() - 1; j++)
                    {
                        cout << " OR " << cmdLineFile[i].arg->at(j);
                    }
                }
            }
            cout << std::endl;
        }
        for (i = 0; cmdLinePara[i].str[0] != '\0'; i++) 
        {
            cout <<"\t\t" << default_commandsep << cmdLinePara[i].str 
                 << (*cmdLineType)[cmdLinePara[i].type];
     
            if (cmdLinePara[i].type == OPTARG) 
            {
                if (cmdLinePara[i].arg !=  NULL)
                {
                    
                    cout << "=" << cmdLinePara[i].arg->at(0);
                    for(j = 1; j < (int)cmdLinePara[i].arg->size() - 1; j++)
                    {
                        cout << " OR " << cmdLinePara[i].arg->at(j);
                    }
                }
            }
            cout << std::endl;
        }
        throw 1;;
    }


//*****************************************************************************
//  Check to see if sequence type is explicitely stated..override ************
// the automatic checking (DNA or Protein).   /type=d or /type=p *************
//****************************************************************************

    if(setType != -1)
    {
        string msg;
        if(((*paramArg)[setType].length()) > 0) 
        {
            temp = findMatch((*paramArg)[setType], typeArg, 2);
            if(temp == 0) 
            { 
                userParameters->setDNAFlag(false);
                userParameters->setExplicitDNAFlag(true);
                /*msg = "Sequence type explicitly set to Protein";
		  cout << msg << std::endl;*/
            }
            else if(temp == 1) 
            {
                /*msg = "Sequence type explicitly set to DNA";
		  cout << msg << std::endl;*/
                userParameters->setDNAFlag(true);
                userParameters->setExplicitDNAFlag(true);
            }
            else
            {
                msg = "Unknown sequence type " + (*paramArg)[setType];
                cerr << std::endl << msg << endl;
            }
            #if DEBUGFULL 
                if(logObject && DEBUGLOG)
                {
                    logObject->logMsg(msg);
                }
            #endif            
        }
    }


//***************************************************************************
//   check to see if 1st parameter does not start with '/' i.e. look for an *
//   input file as first parameter.   The input file can also be specified  *
//   by /infile=fname.                                                      *
//***************************************************************************
// JULIE - moved to checkParam()
  //  if(paramstr[0] != '/') 
  //  {
  //      strcpy(seqName, params[0]);
  //  }


//*************************************************
//  Look for /infile=file.ext on the command line *
//*************************************************

    if(setInfile != -1) 
    {
        if((*paramArg)[setInfile].length() <= 0) 
        {
            exitWithErrorMsg("Bad sequence file name");
        }
        userParameters->setSeqName((*paramArg)[setInfile]);

        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                string msg = "Sequence file name (seqName) has been set to ";
                msg += userParameters->getSeqName();             
                logObject->logMsg(msg);
            }
        #endif        
    }
    
    // NOTE keep an eye on this part to see if it works.
    if(userParameters->getSeqName() != "") 
    {
        // NOTE I will need to cheack if it has been successful before setting
        // doSomething to true.
        int code;
        code = clustalObj->commandLineReadSeq(1, input);
        if(code == OK)
        {
            doSomething = true;
        }
        else
        {
            doSomething = false;
            throw 1;
        }
    }

    // NOTE call the other function to set the parameter values!!
    setOptionalParam((*input).substitutionMatrix);

//*********************************************************
// Look for /profile1=file.ext  AND  /profile2=file2.ext *
// You must give both file names OR neither.             *
//********************************************************


    if(setProfile1 != -1) 
    {
        if((*paramArg)[setProfile1].length() <= 0) 
        {
            exitWithErrorMsg("Bad profile 1 file name");
        }
                
        clustalObj->profile1Input((*paramArg)[setProfile1]);
    }

    if(setProfile2 != -1) 
    {
        if((*paramArg)[setProfile2].length() <= 0) 
        {
            exitWithErrorMsg("Bad profile 2 file name");
        }
        if(userParameters->getProfile1Empty()) 
        {
            exitWithErrorMsg("Only 1 profile file (profile 2) specified.");
        }
        
        clustalObj->profile2Input((*paramArg)[setProfile2]);
        doSomething = doProfileAlign = true; 
    }

//************************************************************************
// Look for /tree or /bootstrap or /align or /usetree ******************
//************************************************************************

    if (setBatch != -1)
    {
        userParameters->setInteractive(false);
    }

    if (setInteractive != -1)
    {
        userParameters->setInteractive(true);
        if (setQuiet!=-1) {
            cout << "interactive menu: overriding " << default_commandsep << "quiet" << std::endl;
            userParameters->setDisplayInfo(true);
            utilityObject->beQuiet(false);
        }
    }

    if (userParameters->getInteractive()) 
    {
        setTree = -1;
        setBootstrap = -1;
        setAlign = -1;
        setUseTree = -1;
        setUseTree1 = -1;
        setUseTree2 = -1;
        setNewTree = -1;
        setConvert = -1;
    }

    if(setTree != -1 )
    {
        if(userParameters->getEmpty()) 
        {
            exitWithErrorMsg("Cannot draw tree.  No input alignment file");
        }
        else
        { 
            doTreeFromAlign = true;
        }
    }

    if(setBootstrap != -1)
    {
        if(userParameters->getEmpty()) 
        {
            exitWithErrorMsg("Cannot bootstrap tree. No input alignment file");
        }
        else 
        {
            temp = 0;
            // Check if there is anything in the string!
            if((*paramArg)[setBootstrap].length() > 0)
            {
                if (sscanf((*paramArg)[setBootstrap].c_str(), "%d", &temp) != 1) 
                {
                    reportBadOptionAndExit("bootstrap", "integer");
                }
            }
            if(temp > 0)
            {    
                userParameters->setBootNumTrials(temp);
            }
            doBootstrap = true;
        }
    }

    if(setAlign != -1)
    {
        if(userParameters->getEmpty()) 
        {
            exitWithErrorMsg("Cannot align sequences.  No input file");
        }
        else
        { 
            doAlign = true;
        }
    }

    if(setConvert != -1)
    {
        if(userParameters->getEmpty()) 
        {
            exitWithErrorMsg("Cannot convert sequences.  No input file");
        }
        else
        { 
            doConvert = true;
        }
    }
 
    if(setUseTree != -1)
    {
        if(userParameters->getEmpty()) 
        {
            exitWithErrorMsg("Cannot align sequences.  No input file");
        }
        else  
        {
           if((*paramArg)[setUseTree].length() == 0) 
           {
               exitWithErrorMsg("Cannot align sequences.  No tree file specified");
           }
           else 
           {
               phylipTreeName = (*paramArg)[setUseTree];
           }
           userParameters->setUseTreeFile(true);
           doAlignUseOldTree = true;
        }
    }

    if(setNewTree != -1)
    {
        if(userParameters->getEmpty()) 
        {
            exitWithErrorMsg("Cannot align sequences.  No input file");
        }
        else  
        {
            if((*paramArg)[setNewTree].length() == 0) 
            {
                exitWithErrorMsg("Cannot align sequences.  No tree file specified");
            }
            else 
            {
                phylipTreeName = (*paramArg)[setNewTree];
            }
            userParameters->setNewTreeFile(true);
            doGuideTreeOnly = true;
        }
    }
 
    if(setUseTree1 != -1)
    {
        if(userParameters->getProfile1Empty()) 
        {
            exitWithErrorMsg("Cannot align profiles.  No input file");
        }
        else if(profileType == SEQUENCE) 
        {
            reportInvalidOptionAndExit("usetree1");
        }
        else  
        {
            if((*paramArg)[setUseTree1].length() == 0) 
            {
                exitWithErrorMsg("Cannot align profiles.  No tree file specified");
            }
            else 
            {
                p1TreeName = (*paramArg)[setUseTree1];
            }
            userParameters->setUseTree1File(true);
            doAlignUseOldTree = true;
        }
    }

    if(setNewTree1 != -1)
    {
        if(userParameters->getProfile1Empty()) 
        {
            exitWithErrorMsg("Cannot align profiles.  No input file");
        }
        else if(profileType == SEQUENCE) 
        {
            reportInvalidOptionAndExit("newtree1");
        }
        else  
        {
           if((*paramArg)[setNewTree1].length() == 0) 
           {
               exitWithErrorMsg("Cannot align profiles. No tree file specified");
           }
           else 
           {
                p1TreeName = (*paramArg)[setNewTree1];
           }
           userParameters->setNewTree1File(true);
        }
    }
 
    if(setUseTree2 != -1)
    {
        if(userParameters->getProfile2Empty()) 
        {
            exitWithErrorMsg("Cannot align profiles.  No input file");
        }
        else if(profileType == SEQUENCE) 
        {
            reportInvalidOptionAndExit("usetree2");
        }
        else  
        {
            if((*paramArg)[setUseTree2].length() == 0) 
            {
                exitWithErrorMsg("Cannot align profiles.  No tree file specified");
            }
            else 
            {
                p2TreeName = (*paramArg)[setUseTree2];
            }
            userParameters->setUseTree2File(true);
            doAlignUseOldTree = true;
        }
    }

    if(setNewTree2 != -1)
    {
        if(userParameters->getProfile2Empty()) 
        {
            exitWithErrorMsg("Cannot align profiles.  No input file");
        }
        else if(profileType == SEQUENCE) 
        {
            reportInvalidOptionAndExit("newtree2");
        }
        else  
        {
            if((*paramArg)[setNewTree2].length() == 0) 
            {
                exitWithErrorMsg("Cannot align profiles.  No tree file specified");
            }
            else 
            {
                p2TreeName = (*paramArg)[setNewTree2];
            }
            userParameters->setNewTree2File(true);
        }
    }


    if( (!doTreeFromAlign) && (!doBootstrap) && (!userParameters->getEmpty()) && (!doProfileAlign) && 
      (!doAlignUseOldTree) && (!doGuideTreeOnly) && (!doConvert))
    { 
        doAlign = true;
    }

//** ? /quicktree 
    if(setQuickTree != -1)
    {
        userParameters->setQuickPairAlign(true);
    }

    // NOTE 
    if(userParameters->getDNAFlag()) 
    {
        userParameters->setDNAParams();
    }
    else 
    {
        userParameters->setProtParams();
    }
    
    if(userParameters->getInteractive()) 
    {
        if (!xmenus)
        { 
            userParameters->setMenuFlag(true);
        }
        return;
    }


    if(!doSomething) 
    {
        exitWithErrorMsg("No input file(s) specified");
    }


    

//***************************************************************************
// Now do whatever has been requested ***************************************
//***************************************************************************
// NOTE This part is obviously not done yet! Functions not working in Clustal!!!!
    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("Now doing the requested task(s)");
        }
    #endif
        
    if(doProfileAlign) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Doing profile align");
            }
        #endif        
        if (profileType == PROFILE)
        {
            #if DEBUGFULL 
                if(logObject && DEBUGLOG)
                {
                    logObject->logMsg("        Calling Profile_align in clustal obj!");
                }
            #endif            
            clustalObj->profileAlign(&p1TreeName, &p2TreeName, output);
        }
        else
        {
            #if DEBUGFULL 
                if(logObject && DEBUGLOG)
                {
                logObject->logMsg("        Calling sequencesAlignToProfile in clustal obj!");
                }
            #endif
            clustalObj->sequencesAlignToProfile(&phylipTreeName, output);
        }
    }

    else if(doAlign)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Doing Alignment");
            }
        #endif
        //cout << "align: " << &phylipTreeName;
        clustalObj->align(&phylipTreeName, output);
    }

    else if(doConvert) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Doing filetype conversion");
            }
        #endif
        //cout << "align: " << &phylipTreeName;
        clustalObj->outputNow(output);
    }

    else if (doAlignUseOldTree)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Doing Alignment only");
            }
        #endif                 
        clustalObj->doAlignUseOldTree(&phylipTreeName, output);
    }

    else if(doGuideTreeOnly)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Doing tree only");
            }
        #endif               
        clustalObj->doGuideTreeOnly(&phylipTreeName);
    }

    else if(doTreeFromAlign)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {   
                logObject->logMsg("    Doing tree");
            }
        #endif                  
        clustalObj->phylogeneticTree(&phylipTreeName, &clustalTreeName, &distTreeName,
                                     &nexusTreeName, pimName);
    }

    else if(doBootstrap)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Doing Bootstrap");
            }
        #endif            
            clustalObj->bootstrapTree(&phylipTreeName, &clustalTreeName, &nexusTreeName);
    }

    //cout << std::endl;
    return;
}



int CommandLineParser::checkParam(StringArray* args, StringArray* params,
                                  StringArray* paramArg)
{
    int len, i, j, nameFirst, num;
    //int k;
    vector<int> match;
    match.resize(MAXARGS);
    bool name1 = false;
    
    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("Checking Parameters!");
        }
    #endif
        
    if(args->size() == 0) // NOTE I think this will work better!
    {
        cout << "The argument list is empty\n";
        return 0;
    }

    // AW: first arg is an input file if it doesnt start with commandsep
    if (VALID_COMMAND_SEP.find((*args)[0][0], 0) == string::npos)
    {
        name1 = true;
        params->push_back((*args)[0]); // Put the string at the back of this vector.
        paramArg->push_back(""); // Push a blank string onto the first element.
    }
    else // It is not a file name
    {
        params->push_back((*args)[0].substr(1));
    }
    
    for (i = 1; i < (int)args->size(); i++) 
    {
        params->push_back(""); // Empty string!

        for(j = 0; j < (int)(*args)[i].length() - 1; j++)
        {
            if(isprint((*args)[i][j + 1])) // Character printable?
            {
                // We start at j + 1 because each option should begin with commandsep
                (*params)[i].append((*args)[i].substr(j + 1, 1));
            }
        }
    }
    num = i;


    // AW:
    // params are now setup
    // extract paramArgs in next step

    
    if ((int)args->size() > MAXARGS)
    {
        cerr << "Error: too many command line arguments\n";
        return(-1);
    }
    /*
      special case - first parameter is input fileName
    */
    nameFirst = 0;
    if(name1 == true) // If name of file is first argument
    {
        userParameters->setSeqName((*params)[0]);
        /* Andreas Wilm (UCD) 2008-03-19:
           conversion nowadays unnecessary and makes trouble
           
           /@  JULIE
           convert to lower case now
           @/
           #ifndef UNIX
           if(logObject && DEBUGLOG)
           {
           logObject->logMsg("Converting seqName to lower case.\n");
           }         
           string temp = ConvertStringToLower(userParameters->getSeqName());
           userParameters->setSeqName(temp);
           #endif
        */
        nameFirst = 1;
    }
  
    // NOTE if name first we should start at the 2nd element in paramArg
    // This loop is used to set up the paramArg vector!
    for (i = nameFirst; i < num; i++) 
    {
        bool has_arg=false;
        paramArg->push_back(""); // Push a new empty string on.
        len = (*params)[i].length();
        for(j = 0; j < len; j++)
        {
            if((*params)[i][j] == '=') 
            {
                has_arg=true;
                (*paramArg)[i].assign((*params)[i].substr(j + 1, len - j -1));               
                // Trim off the bit from the '=' to the end, and put all in lower case!
                (*params)[i].assign(ConvertStringToLower((*params)[i].substr(0, j)));
                break;
            }
        }
        // Andreas Wilm (UCD): 2008-03-19:
        // this convert nonarg params to lowercase (-QuIcKtReE etc)
        if (!has_arg) {
            (*params)[i].assign(ConvertStringToLower((*params)[i]));          
        }
    }

    if(paramArg->size() != params->size())
    {
        cerr << "There is something wrong with arguments. Lengths different\n";
        return -1;
    }
    
    /*
      for each parameter given on the command line, first search the list of recognised
      optional parameters....
    */
    for (i = 0; i < num; i++) 
    {
        if ((i == 0) && (name1 == true))
        { 
            continue;
        }
        j = 0;
        match[i] = -1;
        for(;;) 
        {
            if (cmdLinePara[j].str[0] == '\0')
            { 
                // Think this means we have not found it!
                break;
            }
            if (!(*params)[i].compare(cmdLinePara[j].str))
            {
                match[i] = j; // Match has been found!
                *cmdLinePara[match[i]].flag = i;
                
                if ((cmdLinePara[match[i]].type != NOARG) && ((*paramArg)[i] == "")) 
                {
                    cerr <<  "Error: parameter required for " << default_commandsep << (*params)[i] << endl;
                    return -1;
                
                /* Andreas Wilm (UCD) 2008-03-19:
                 *  conversion nowadays unnecessary and breaks things
                 *  
                 * //  JULIE
                 * //    convert parameters to lower case now, unless the parameter is a fileName
                 * #ifdef UNIX
                 * else if (cmdLinePara[match[i]].type != FILARG && (*paramArg)[i] != "")
                 * #endif
                */

                } else if (cmdLinePara[match[i]].type != FILARG && (*paramArg)[i] != "") {
                    if ((*paramArg)[i] != "") 
                    {
                        // lowercase arg if not a filename to support mixed case
                        (*paramArg)[i].assign(ConvertStringToLower((*paramArg)[i]));
                    }
                }
                break;
            }
            j++;
        }
    }
    
    /*
      ....then the list of recognised input files,.... 
    */
    for (i = 0; i < num; i++) 
    {
        if ((i == 0) && (name1 == true)) 
        {
            continue;
        }
        if (match[i] != -1) 
        {
            continue;
        }
        j = 0;
        for(;;) 
        {
            if (cmdLineFile[j].str[0] == '\0') 
            {
                // Have not found a match!
                break;
            }
            if (!(*params)[i].compare(cmdLineFile[j].str)) 
            {
                match[i] = j;
                *cmdLineFile[match[i]].flag = i;
                if ((cmdLineFile[match[i]].type != NOARG) &&
                                    ((*paramArg)[i] == "")) 
                {
                    cerr << "Error: parameter required for " << default_commandsep << (*params)[i] << endl;
                    return -1;
                }
                break;
            }
            j++;
        }
    }
    
    /*
      ....and finally the recognised verbs. 
    */
    for (i = 0; i < num; i++) 
    {
        if ((i == 0) && (name1 == true)) 
        {
            continue;
        }
        if (match[i] != -1) 
        {
            continue;
        }
        j = 0;
        for(;;) 
        {
            if (cmdLineVerb[j].str[0] == '\0') 
            {
                // Havent found it!
                break;
            }
            if (!(*params)[i].compare(cmdLineVerb[j].str)) 
            {
                match[i] = j;
                *cmdLineVerb[match[i]].flag = i;
                if ((cmdLineVerb[match[i]].type != NOARG) && ((*paramArg)[i] == "")) 
                {
                    cerr << "Error: parameter required for " << default_commandsep << (*params)[i] << endl;
                    return -1;
                }
                break;
            }
            j++;
        }
    }

    /*
      check for any unrecognised parameters.
    */
    for (i = 0; i < num; i++) 
    {
        if (match[i] == -1) 
        {
            cerr << "Error: unknown option " << default_commandsep << (*params)[i] << endl;
            return -1;
        }
    }
    return(num);
}

void CommandLineParser::setOptionalParam(Rcpp::NumericMatrix substitutionMatrix)
{
    int temp;
    int _ktup, _windgap, _signif, _window = 0;
    string _matrixname;
    //int c, i;
    float ftemp;
    //char tstr[100];

    #if DEBUGFULL 
        if(logObject && DEBUGLOG)
        {
            logObject->logMsg("Setting optional parameters.");
        }
    #endif      
    //****************************************************************************
    //* look for parameters on command line  e.g. gap penalties, k-tuple etc.    *
    //****************************************************************************
  
    // Mark change 16-2-2007. 
    
    if(setNumIterations != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting num iterations parameter.");
            }    
        #endif      
        temp = 0;
        if((*paramArg)[setNumIterations].length() > 0)
        {
            if (sscanf((*paramArg)[setNumIterations].c_str(),"%d", &temp) != 1) 
            {
                reportBadOptionAndExit("numiter", "int");
                temp = 0;
            }
        }
        if(temp > 0)
        {
            userParameters->setNumIterations(temp);
        }
        else
        {
            exitWithErrorMsg("Cannot use a negative value for number of iterations.");
        }
    }

    
    //** ? /score=percent or /score=absolute *
    if(setScore != -1)
    {    
        if((*paramArg)[setScore].length() > 0) 
        {
            temp = findMatch((*paramArg)[setScore], scoreArg, 2);
            if(temp == 0)
            {
                userParameters->setPercent(true);
                #if DEBUGFULL 
                    if(logObject && DEBUGLOG)
                    {
                        logObject->logMsg("    Setting score parameter = percent");
                    }
                #endif                 
            }
            else if(temp == 1)
            {
                userParameters->setPercent(false);
                #if DEBUGFULL 
                    if(logObject && DEBUGLOG)
                    {
                        logObject->logMsg("    Setting score parameter = absolute");
                    }
                #endif                 
            }
            else
            {
                cerr << "\nUnknown SCORE type: " << (*paramArg)[setScore] << endl;
                #if DEBUGFULL 
                    if(logObject && DEBUGLOG)
                    {
                        logObject->logMsg("    problem setting score type!!!!!");
                    }
                #endif                
            }
        }
    }
    // NOTE I decided to stay with sscanf for getting the int. Options in c++ no better.
    //** ? /seed=n *
    if(setSeed != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting seed parameter.");
            }    
        #endif      
        temp = 0;
        if((*paramArg)[setSeed].length() > 0)
        {
            if (sscanf((*paramArg)[setSeed].c_str(),"%d",&temp) != 1) 
            {
                reportBadOptionAndExit("seed", "integer");
            }
        }
        if(temp > 0)
        {
            userParameters->setBootRanSeed(temp);
        }
        //cout<< "\ntemp = " << temp << "; seed = " << userParameters->getBootRanSeed()
        //    << ";\n";
    }
  
    if(setTreeAlgorithm != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting clustering algorithm parameter.");
            }
        #endif
        if((*paramArg)[setTreeAlgorithm].length() > 0) 
        { 
            temp = findMatch((*paramArg)[setTreeAlgorithm], clusterAlgorithm, 2);
            if(temp == 0)
            {
                userParameters->setClusterAlgorithm(NJ);
            }
            else if(temp == 1)
            {
                userParameters->setClusterAlgorithm(UPGMA);
            }
            else
            {
                cerr << "Unknown option for clustering algorithm. Using default\n";
                userParameters->setClusterAlgorithm(NJ);
            }
        }
    }
    
    
    //** ? /output=PIR, GCG, GDE or PHYLIP *
    if(setOutput != -1)
    {
        
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting output parameter.");
            }
        #endif          
        if((*paramArg)[setOutput].length() > 0) 
        {
            temp = findMatch((*paramArg)[setOutput], outputArg, 7);
            if (temp >= 0 && temp <= 6) 
            {
                userParameters->setOutputClustal(false);
                userParameters->setOutputGCG(false);
                userParameters->setOutputPhylip(false);
                userParameters->setOutputNbrf(false);
                userParameters->setOutputGde(false);
                userParameters->setOutputNexus(false);
                userParameters->setOutputFasta(false);
            }
            switch (temp) 
            {
                case 0: // GCG 
                    userParameters->setOutputGCG(true);
                    break;
                case 1: // GDE 
                    userParameters->setOutputGde(true);
                    break;
                case 2: // PIR 
                    userParameters->setOutputNbrf(true);
                    break;
                case 3: // PHYLIP
                    userParameters->setOutputPhylip(true);
                    break;
                case 4: // NEXUS
                    userParameters->setOutputNexus(true);
                    break;
                case 5: // FASTA
                    userParameters->setOutputFasta(true);
                    break;
                case 6: // CLUSTAL
                    userParameters->setOutputClustal(true);
                    break;
                default:
                    // FIXME AW: 1.83 behaves the same, but shouldnt
                    // we exit here?
                    exitWithErrorMsg("Unknown OUTPUT type: " + (*paramArg)[setOutput]);
            }
        }
    }
    //** ? /outputtree=NJ or PHYLIP or DIST or NEXUS 
    if(setOutputTree != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting outputtree parameter.");
            }   
        #endif           
        if((*paramArg)[setOutputTree].length() > 0) 
        {
            temp = findMatch((*paramArg)[setOutputTree], outputTreeArg, 4);
            switch (temp) 
            {
                case 0: // NJ 
                    userParameters->setOutputTreeClustal(true);
                    break;
                case 1: // PHYLIP
                    userParameters->setOutputTreePhylip(true);
                    break;
                case 2: // DIST
                    userParameters->setOutputTreeDistances(true);
                    break;
                case 3: // NEXUS 
                    userParameters->setOutputTreeNexus(true);
                    break;
                default:
                    cerr << "\nUnknown OUTPUT TREE type: " 
                         << (*paramArg)[setOutputTree] << endl;
            }
        }
    }
    
    //** ? /profile (sets type of second input file to profile) 
    if(setProfile != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting profileType = PROFILE.");
            }         
        #endif 
        profileType = PROFILE;
    }
  
    //** ? /sequences (sets type of second input file to list of sequences)
    if(setSequences != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting profileType = SEQUENCE.");
            }         
        #endif
        profileType = SEQUENCE;
    } 
  
  
    //** ? /ktuple=n
    if(setKtuple != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting ktup parameter.");
            }
        #endif         
        _ktup = 0;
        if((*paramArg)[setKtuple].length() > 0)
        {
            if (sscanf((*paramArg)[setKtuple].c_str(),"%d",&_ktup)!=1) 
            {
                reportBadOptionAndExit("ktuple", "integer");
                _ktup = 0;
            }
        }
        
        if(_ktup > 0) 
        {
            if(userParameters->getDNAFlag()) 
            {
                if(_ktup <= 4) 
                {
                    userParameters->setKtup(_ktup);
                    userParameters->setDNAKtup(_ktup);
                    userParameters->setWindowGap(_ktup + 4);
                    userParameters->setDNAWindowGap(_ktup + 4);
                } else {
                    // see comment in bug 185
                    cerr << "WARNING: Ignoring invalid ktuple of " << _ktup <<  " (must be <=4)" << std::endl;
                }
            }
            else 
            {
                if(_ktup <= 2) 
                {
                    userParameters->setKtup(_ktup);
                    userParameters->setAAKtup(_ktup);
                    userParameters->setWindowGap(_ktup + 3);
                    userParameters->setAAWindowGap(_ktup + 3);
                    // AW: why set setDNAWindowGap? we are in AA mode
                    // userParameters->setDNAWindowGap(_ktup + 4);
                } else {
                    // see comment in bug 185
                   cerr << "WARNING: Ignoring invalid ktuple of " << _ktup << " (must be <=2)" << std::endl;
                }
            }
        }
    }
    
    //** ? /pairgap=n 
    if(setPairGap != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting pairgap parameter.");
            }
        #endif         
        _windgap = 0;
        if((*paramArg)[setPairGap].length() > 0)
        {
            if (sscanf((*paramArg)[setPairGap].c_str(),"%d",&_windgap)!=1) 
            {
                reportBadOptionAndExit("pairgap", "integer");
                _windgap = 0;
            }
            if(_windgap > 0)
            {
                if(userParameters->getDNAFlag()) 
                {
                    if(_windgap > userParameters->getKtup()) 
                    {
                        userParameters->setWindowGap(_windgap);
                        userParameters->setDNAWindowGap(_windgap);
                    }
                }
                else 
                {
                    if(_windgap > userParameters->getKtup()) 
                    {
                        userParameters->setWindowGap(_windgap);
                        userParameters->setAAWindowGap(_windgap);
                    }
                }
            }
        }
    } 
  
    //** ? /topdiags=n  
    if(setTopDiags != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting topdiags parameter.");
            }
        #endif          
        _signif = 0;
        if((*paramArg)[setTopDiags].length() > 0)
        {
            if (sscanf((*paramArg)[setTopDiags].c_str(),"%d",&_signif)!=1) 
            {
                reportBadOptionAndExit("topdiags", "integer");
            }
        }
        if(_signif > 0)
        {
            if(userParameters->getDNAFlag()) 
            {
                if(_signif > userParameters->getKtup()) 
                {
                    userParameters->setSignif(_signif);
                    userParameters->setDNASignif(_signif);
                }
            }
            else 
            {
                if(_signif > userParameters->getKtup()) 
                {
                    userParameters->setSignif(_signif);
                    userParameters->setAASignif(_signif);
                }
            }
        }
    }
    

    //** ? /window=n 
    if(setWindow != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting window parameter.");
            }
        #endif        
        _window = 0;
        if((*paramArg)[setWindow].length() > 0)
        {
            if (sscanf((*paramArg)[setWindow].c_str(),"%d",&_window)!=1) 
            {
                reportBadOptionAndExit("window", "integer");
                _window = 0;
            }
        }
        if(_window > 0)
        {
            if(userParameters->getDNAFlag()) 
            {
                if(_window > userParameters->getKtup()) 
                {
                    userParameters->setWindow(_window);
                    userParameters->setDNAWindow(_window);
                }
            }
            else 
            {
                if(_window > userParameters->getKtup()) 
                {
                    userParameters->setWindow(_window);
                    userParameters->setAAWindow(_window);
                }
            }
        }
    }
  
    //** ? /kimura 
    if(setKimura != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting kimura=true");
            }
        #endif        
        userParameters->setKimura(true);
    }
  
    //** ? /tossgaps 
    if(setTossGaps != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting tossgaps=true");
            }
        #endif        
        userParameters->setTossGaps(true);
    }
    
    //** ? /negative  
    if(setNegative != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting useNegMatrix=true");
            }         
        #endif 
        userParameters->setUseNegMatrix(true);
    }
  
    //** ? /noweights
    if(setNoWeights!= -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting noweights=true");
            }       
        #endif 
        userParameters->setNoWeights(true);
    }
  
  
    //** ? /pwmatrix=ID (user's file) 
    if(setPWMatrix != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting pwmatrix parameter.");
            }        
        #endif
        temp = (*paramArg)[setPWMatrix].length();
        if(temp > 0) 
        {
            _matrixname = ConvertStringToLower((*paramArg)[setPWMatrix]);
            if (_matrixname.compare("blosum") == 0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 1, Protein, Pairwise);
            }
            else if (_matrixname.compare("pam") == 0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 2, Protein, Pairwise);
            }
            else if (_matrixname.compare("gonnet") == 0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 3, Protein, Pairwise);
            }
            else if (_matrixname.compare("id") == 0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 4, Protein, Pairwise);
            }
            else 
            {
                if (strcmp(DUMMY_R_MATRIX, (*paramArg)[setPWMatrix].c_str()) == 0) {
                	if(subMatrix->getUserMatFromR(substitutionMatrix, Protein, Pairwise))
					{
						subMatrix->setCurrentNameAndNum((*paramArg)[setPWMatrix], 5,
														 Protein, Pairwise);
						pwUserMatrixName = (*paramArg)[setPWMatrix];
					} else
	                	throw 1;
                } else {
					char hackTempName[FILENAMELEN + 1];
					strcpy(hackTempName, (*paramArg)[setPWMatrix].c_str());

					if(subMatrix->getUserMatFromFile(hackTempName, Protein, Pairwise))
					{
						subMatrix->setCurrentNameAndNum((*paramArg)[setPWMatrix], 5,
														 Protein, Pairwise);
						pwUserMatrixName = (*paramArg)[setPWMatrix];
					} else
	                	throw 1;
                }
            }
        }
    }

        
    //** ? /matrix=ID (user's file)
    if(setMatrix != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting matrix parameter.");
            }       
        #endif 
        temp = (*paramArg)[setMatrix].length();
        if(temp > 0) 
        {
        	_matrixname = ConvertStringToLower((*paramArg)[setMatrix]);

            if (_matrixname.compare("blosum")==0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 1, Protein, MultipleAlign);
            }
            else if (_matrixname.compare("pam")==0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 2, Protein, MultipleAlign);
            }
            else if (_matrixname.compare("gonnet")==0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 3, Protein, MultipleAlign);
            }
            else if (_matrixname.compare("id")==0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 4, Protein, MultipleAlign);
            }
            else 
            {
            	if (strcmp(DUMMY_R_MATRIX, (*paramArg)[setMatrix].c_str()) == 0) {
            		if(subMatrix->getUserMatSeriesFromR(substitutionMatrix))
					{
						subMatrix->setCurrentNameAndNum((*paramArg)[setMatrix], 4,
														 Protein, MultipleAlign);
						userMatrixName = (*paramArg)[setMatrix];
					} else throw 1;
            	} else {
					char hackTempName[FILENAMELEN + 1];
					strcpy(hackTempName, (*paramArg)[setMatrix].c_str());

					if(subMatrix->getUserMatSeriesFromFile(hackTempName))
					{
						subMatrix->setCurrentNameAndNum((*paramArg)[setMatrix], 4,
														 Protein, MultipleAlign);
						userMatrixName = (*paramArg)[setMatrix];
					} else
	                	throw 1;
            	}
            }
        }
    }

    //** ? /pwdnamatrix=ID (user's file)
    if(setPWDNAMatrix != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting pwdnamatrix parameter.");
            }       
        #endif 
        temp = (*paramArg)[setPWDNAMatrix].length();
        if(temp > 0) 
        {
            _matrixname = ConvertStringToLower((*paramArg)[setPWDNAMatrix]);

            if (_matrixname.compare("iub") == 0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 1, DNA, Pairwise);
            }
            else if (_matrixname.compare("clustalw") == 0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 2, DNA, Pairwise);
            }
            else 
            {
            	if (strcmp(DUMMY_R_MATRIX, (*paramArg)[setPWDNAMatrix].c_str()) == 0) {
            		if(subMatrix->getUserMatFromR(substitutionMatrix, DNA, Pairwise))
					{
						subMatrix->setCurrentNameAndNum((*paramArg)[setPWDNAMatrix], 3,
														 Protein, Pairwise);
						pwDNAUserMatrixName = (*paramArg)[setPWDNAMatrix];
					}
					else
						throw 1;
            	} else {
					char hackTempName[FILENAMELEN + 1];
					strcpy(hackTempName, (*paramArg)[setPWDNAMatrix].c_str());

					if(subMatrix->getUserMatFromFile(hackTempName, DNA, Pairwise))
					{
						subMatrix->setCurrentNameAndNum((*paramArg)[setPWDNAMatrix], 3,
														 Protein, Pairwise);
						pwDNAUserMatrixName = (*paramArg)[setPWDNAMatrix];
					}
					else
						throw 1;
            	}
            }
        }
    }

    //** ? /dnamatrix=ID (user's file)
    if(setDNAMatrix != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting dnamatrix parameter.");
            }        
        #endif 
        temp = (*paramArg)[setDNAMatrix].length();
        if(temp > 0) 
        {
            _matrixname = ConvertStringToLower((*paramArg)[setDNAMatrix]);

            if (_matrixname.compare("iub") == 0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 1, DNA, MultipleAlign);
            }
            else if (_matrixname.compare("clustalw") == 0) 
            {
                subMatrix->setCurrentNameAndNum(_matrixname, 2, DNA, MultipleAlign);
            }
            else 
            {
            	if (strcmp(DUMMY_R_MATRIX, (*paramArg)[setDNAMatrix].c_str()) == 0) {
            		if(subMatrix->getUserMatFromR(substitutionMatrix, DNA, MultipleAlign))
					{
						subMatrix->setCurrentNameAndNum((*paramArg)[setDNAMatrix], 3,
														 Protein, MultipleAlign);
						DNAUserMatrixName = (*paramArg)[setDNAMatrix];
					}
					else
						throw 1;
            	} else {
					char hackTempName[FILENAMELEN + 1];
					strcpy(hackTempName, (*paramArg)[setDNAMatrix].c_str());

					if(subMatrix->getUserMatFromFile(hackTempName, DNA, MultipleAlign))
					{
						subMatrix->setCurrentNameAndNum((*paramArg)[setDNAMatrix], 3,
														 Protein, MultipleAlign);
						DNAUserMatrixName = (*paramArg)[setDNAMatrix];
					}
					else
						throw 1;
            	}
            }
        }
    }
 
       
    //** ? /maxdiv= n
    if(setMaxDiv != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting maxdiv parameter.");
            }       
        #endif 
        temp = 0;
        if((*paramArg)[setMaxDiv].length() > 0)
        {
            if (sscanf((*paramArg)[setMaxDiv].c_str(),"%d",&temp)!=1) 
            {
                reportBadOptionAndExit("maxdiv", "integer");
                temp = 0;
            }
        }
        if (temp >= 0)
        {
            userParameters->setDivergenceCutoff(temp);
        }
    }

    //** ? /gapdist= n
    if(setGapDist != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting gapdist parameter.");
            }       
        #endif
         
        temp = 0;
        if((*paramArg)[setGapDist].length() > 0)
            if (sscanf((*paramArg)[setGapDist].c_str(),"%d",&temp)!=1) 
            {
                reportBadOptionAndExit("gapdist", "integer");
            }
        if (temp >= 0)
        {
            userParameters->setGapDist(temp);
        }
    }

    //** ? /debug= n 
    if(setDebug != -1) 
    {
        temp = 0;
        if((*paramArg)[setDebug].length() > 0)
            if (sscanf((*paramArg)[setDebug].c_str(),"%d",&temp)!=1) 
            {
                reportBadOptionAndExit("debug", "integer");
            }
        if (temp >= 0)
        {
            userParameters->setDebug(temp);
        }
    }

    //** ? /outfile= (user's file) 
    if(setOutfile != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting outfile parameter.");
            }        
        #endif
        if((*paramArg)[setOutfile].length() > 0) 
        {
            userParameters->setOutfileName((*paramArg)[setOutfile]);
        }
    }

    //*** ? /case= lower/upper 
    if(setCase != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting case parameter.");
            }       
        #endif 
        if((*paramArg)[setCase].length() > 0) 
        {
            temp = findMatch((*paramArg)[setCase], caseArg, 2);
            if(temp == 0) 
            {
                userParameters->setLowercase(true);
            }
            else if(temp == 1) 
            {
                userParameters->setLowercase(false);
            }
            else
            {
                cerr << "\nUnknown case " <<  (*paramArg)[setCase] << endl;
            }
        }
    }

    //*** ? /seqnos=off/on 
    if(setSeqNo != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting seqnos parameter.");
            }         
        #endif
        if((*paramArg)[setSeqNo].length() > 0) 
        {
            temp = findMatch((*paramArg)[setSeqNo], seqNoArg, 2);
            if(temp == 0) 
            {
                userParameters->setClSeqNumbers(false);
            }
            else if(temp == 1) 
            {
                userParameters->setClSeqNumbers(true);
            }
            else
            {
                cerr << "\nUnknown SEQNO option " << (*paramArg)[setSeqNo] << endl;
            }
        }
    }


    if(setSeqNoRange != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting seqno range parameter.");
            }        
        #endif
        if((*paramArg)[setSeqNoRange].length() > 0) 
        {
            temp = findMatch((*paramArg)[setSeqNoRange], seqNoRangeArg, 2);
            //cout << "\n comparing  "
            //     << "\nparamArg[setSeqNoRange]= " << (*paramArg)[setSeqNoRange]
            //     << "\n comparing \n ";

            if(temp == 0) 
            {
                userParameters->setSeqRange(false);
            }
            else if(temp == 1) 
            {
                userParameters->setSeqRange(true);
            }
            else
            {
                cerr << "\nUnknown Sequence range  option " 
                     << (*paramArg)[setSeqNoRange] << endl;
            }
        }
    }

    //*** ? /range=n:m
    if(setRange != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting range parameter.");
            }        
        #endif
        temp = 0;
        if((*paramArg)[setRange].length() > 0)
        {
            // NOTE I have made a big change here! Mark march 14th 2006. This was being done
            // in the Alignment output functions. 
            int iFirstRes = -1; 
            int iLastRes = -1;
            char ignore;
            
            if (sscanf((*paramArg)[setRange].c_str(), "%d%[ :,-]%d", &iFirstRes,
                        &ignore, &iLastRes) != 3) 
            {
                cerr << "setRange:  Syntax Error: Cannot set range, should be from:to \n";
            }
            else
            {
                userParameters->setRangeFrom(iFirstRes);
                userParameters->setRangeTo(iLastRes);
            }
        }
    }

    //*** ? /gapopen=n 
    if(setGapOpen != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting gapopen parameter.");
            }        
        #endif
        ftemp = 0.0;
        if((*paramArg)[setGapOpen].length() > 0)
        {
            if (sscanf((*paramArg)[setGapOpen].c_str(),"%f",&ftemp) != 1) 
            {
                reportBadOptionAndExit("gapopen", "real number");
                ftemp = 0.0;
            }
            if(ftemp >= 0.0)
            {
                if( userParameters->getDNAFlag()) 
                {
                    userParameters->setGapOpen(ftemp);
                    userParameters->setDNAGapOpen(ftemp);
                }
                else 
                {
                    userParameters->setGapOpen(ftemp);
                    userParameters->setProteinGapOpen(ftemp);
                }
            }
        }
    }

    //*** ? /gapext=n
    if(setGapExtend != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting gap extention parameter.");
            }        
        #endif
        
        ftemp = 0.0;
        if((*paramArg)[setGapExtend].length() > 0)
        {
            if (sscanf((*paramArg)[setGapExtend].c_str(),"%f",&ftemp) != 1) 
            {
                reportBadOptionAndExit("gapext", "real number");
               ftemp = 0.0;
            }
            if(ftemp >= 0)
            {
                if(userParameters->getDNAFlag()) 
                {
                    userParameters->setGapExtend(ftemp);
                    userParameters->setDNAGapExtend(ftemp);
                }
                else 
                {
                    userParameters->setGapExtend(ftemp);
                    userParameters->setProteinGapExtend(ftemp);
                }
            }
        }
    }

    //*** ? /transweight=n
    if(setTransWeight != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting transweight parameter.");
            }        
        #endif
        
        ftemp = 0.0;
        if((*paramArg)[setTransWeight].length() > 0)
        {
            if (sscanf((*paramArg)[setTransWeight].c_str(), "%f", &ftemp) != 1) 
            {
                reportBadOptionAndExit("transweight", "real number");
            }
        }
        userParameters->setTransitionWeight(ftemp);
    }

    //*** ? /pwgapopen=n 
    if(setPWGapOpen != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting pwgapopen parameter.");
            }        
        #endif
        
        ftemp = 0.0;
        if((*paramArg)[setPWGapOpen].length() > 0)
        {
            if (sscanf((*paramArg)[setPWGapOpen].c_str(), "%f", &ftemp) != 1) 
            {
                reportBadOptionAndExit("pwgapopen", "real number");
            }
        }
        if(ftemp >= 0.0)
        {
            if(userParameters->getDNAFlag()) 
            {
                userParameters->setPWGapOpen(ftemp);
                userParameters->setDNAPWGapOpenPenalty(ftemp);
            }
            else 
            {
                userParameters->setPWGapOpen(ftemp);
                userParameters->setProteinPWGapOpenPenalty(ftemp);
            }
        }
    }


    //*** ? /gapext=n 
    if(setPWGapExtend != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting pwgapext parameter.");
            }        
        #endif
        
        ftemp = 0.0;
        if((*paramArg)[setPWGapExtend].length() > 0)
        {
            if (sscanf((*paramArg)[setPWGapExtend].c_str(), "%f", &ftemp) != 1) 
            {
                reportBadOptionAndExit("pwgapext", "real number");
            }
        }
        if(ftemp >= 0)
        {
            if(userParameters->getDNAFlag()) 
            {
                userParameters->setPWGapExtend(ftemp);
                userParameters->setDNAPWGapExtendPenalty(ftemp);
            }
            else 
            {
                userParameters->setPWGapExtend(ftemp);
                userParameters->setProteinPWGapExtendPenalty(ftemp);
            }
        }
    }



    //*** ? /outorder=n 
    if(setOutOrder != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting outorder parameter.");
            }        
        #endif
        
        if((*paramArg)[setOutOrder].length() > 0)
        {
            temp = findMatch((*paramArg)[setOutOrder],outOrderArg,2);
        }

        if(temp == 0) 
        {    
            userParameters->setOutputOrder(INPUT);
        }
        else if(temp == 1) 
        {    
            userParameters->setOutputOrder(ALIGNED);
        }
        else
        {
            cerr << "\nUnknown OUTPUT ORDER type " <<  (*paramArg)[setOutOrder] << endl;
        }
    }

    //*** ? /bootlabels=n 
    if(setBootLabels != -1) 
    {

        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting bootlabels parameter.");
            }        
        #endif
        
        if((*paramArg)[setBootLabels].length() > 0)
        {
            temp = findMatch((*paramArg)[setBootLabels], bootLabelsArg, 2);
        }

        if(temp == 0)  
        {    
            userParameters->setBootstrapFormat(BS_NODE_LABELS);
        }
        else if(temp == 1)  
        {    
            userParameters->setBootstrapFormat(BS_BRANCH_LABELS);
        }
        else
        {
            cerr << "\nUnknown bootlabels type " << (*paramArg)[setBootLabels] << endl;
        }
    }

    //*** ? /endgaps 
    if(setUseEndGaps != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting useendgaps=true");
            }        
        #endif
        
        userParameters->setUseEndGaps(false);
    }

    //*** ? /nopgap 
    if(setNoPGap != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting noPrefPenalties=true");
            }        
        #endif
        
        userParameters->setNoPrefPenalties(true);
    }

    //*** ? /nohgap 
    if(setNoHGap != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting nohgap=true");
            }        
        #endif
        
        userParameters->setNoHydPenalties(true);
    }

    //*** ? /novgap 
    if(setNoVGap != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting novgap=false");
            }        
        #endif
        
        userParameters->setNoVarPenalties(false);
    }

    //*** ? /hgapresidues="string" 
    // NOTE I have made some big changes here. It looks as if there was an error here!
    if(setHGapRes != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting hgapresidues parameter.");
            }        
        #endif
        
        userParameters->setHydResidues((*paramArg)[setHGapRes]);
        
    }
              
    //*** ? /nosecstr1 
    if(setSecStruct1 != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting useSS1=false");
            }        
        #endif
        
        userParameters->setUseSS1(false);
    }

    //*** ? /nosecstr2 
    if(setSecStruct2 != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting useSS2=false");
            }        
        #endif
        
        userParameters->setUseSS2(false);
    }

    //*** ? /secstroutput
    if(setSecStructOutput != -1)
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting secstroutput parameter.");
            }        
        #endif
        
        if((*paramArg)[setSecStructOutput].length() > 0) 
        {
            temp = findMatch((*paramArg)[setSecStructOutput], outputSecStrArg, 4);
            if(temp >= 0 && temp <= 3)
            {
                userParameters->setOutputStructPenalties(temp);
            }
            else
            {
                cerr << "\nUnknown case " << (*paramArg)[setSecStructOutput] << endl;
            }
        }
    }


    //*** ? /helixgap= n
    if(setHelixGap != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting helixgap parameter.");
            }        
        #endif
        
        temp = 0;
        if((*paramArg)[setHelixGap].length() > 0)
        {
            if (sscanf((*paramArg)[setHelixGap].c_str(), "%d", &temp) != 1) 
            {
                reportBadOptionAndExit("helixgap", "integer");
            }
        }
        if (temp >= 1 && temp <= 9)
        {
            userParameters->setHelixPenalty(temp);
        }
    }
    
    //*** ? /strandgap= n 
    if(setStrandGap != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting strandgap parameter.");
            }         
        #endif
        
        temp = 0;
        if((*paramArg)[setStrandGap].length() > 0)
        {
            if (sscanf((*paramArg)[setStrandGap].c_str(), "%d", &temp) != 1) 
            {
                reportBadOptionAndExit("strandgap", "integer");
            }
        }
        if (temp >= 1 && temp <= 9)
        {
            userParameters->setStrandPenalty(temp);
        }
    }
    
    //*** ? /loopgap= n 
    if(setLoopGap != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting loopgap parameter.");
            }        
        #endif
        
        temp = 0;
        if((*paramArg)[setLoopGap].length() > 0)
        {
            if (sscanf((*paramArg)[setLoopGap].c_str(), "%d", &temp) != 1) 
            {
                reportBadOptionAndExit("loopgap", "integer");
            }
        }
        if (temp >= 1 && temp <= 9)
        {
            userParameters->setLoopPenalty(temp);
        }
    }

    //*** ? /terminalgap= n
    if(setTerminalGap != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting terminalgap parameter.");
            }        
        #endif
        temp = 0;
        if((*paramArg)[setTerminalGap].length() > 0)
        {
            if (sscanf((*paramArg)[setTerminalGap].c_str(), "%d", &temp) != 1) 
            {
                reportBadOptionAndExit("terminalgap", "integer");
                temp = 0;
            }
        }
        if (temp >= 1 && temp <= 9) 
        {
            userParameters->setHelixEndPenalty(temp);
            userParameters->setStrandEndPenalty(temp);
        }
    }
   
    //*** ? /helixendin= n
    if(setHelixEndIn != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting helixendin parameter.");
            }        
        #endif
        
        temp = 0;
        if((*paramArg)[setHelixEndIn].length() > 0)
        {
            if (sscanf((*paramArg)[setHelixEndIn].c_str(), "%d", &temp) != 1) 
            {
                reportBadOptionAndExit("helixendin", "integer");
                temp = 0;
            }
        }
        if (temp >= 0 && temp <= 3)
        {
            userParameters->setHelixEndMinus(temp);
        }
    }

    //*** ? /helixendout= n
    if(setHelixEndOut != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting helixendout parameter.");
            }        
        #endif
        
        temp = 0;
        if((*paramArg)[setHelixEndOut].length() > 0)
        {
            if (sscanf((*paramArg)[setHelixEndOut].c_str(), "%d", &temp) != 1) 
            {
                reportBadOptionAndExit("helixendout", "integer");
                temp = 0;
            }
        }
        if (temp >= 0 && temp <= 3)
        {
            userParameters->setHelixEndPlus(temp);
        }
    }

    //*** ? /strandendin= n 
    if(setStrandEndIn != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting strandendin parameter.");
            }        
        #endif
        
        temp = 0;
        if((*paramArg)[setStrandEndIn].length() > 0)
        {
            if (sscanf((*paramArg)[setStrandEndIn].c_str(), "%d", &temp) != 1) 
            {
                reportBadOptionAndExit("strandendin", "integer");
            }
        }
        if (temp >= 0 && temp <= 3)
        {
            userParameters->setStrandEndMinus(temp);
        }
    }

    //*** ? /strandendout= n 
    if(setStrandEndOut != -1) 
    {
        #if DEBUGFULL 
            if(logObject && DEBUGLOG)
            {
                logObject->logMsg("    Setting strandendout parameter.");
            }        
        #endif
        
        temp = 0;
        if((*paramArg)[setStrandEndOut].length() > 0)
        {
            if (sscanf((*paramArg)[setStrandEndOut].c_str(), "%d", &temp) != 1) 
            {
                reportBadOptionAndExit("strandendout", "integer");
            }
        }
        if (temp >= 0 && temp <= 3)
        {
            userParameters->setStrandEndPlus(temp);
        }
    }
}

int CommandLineParser::findMatch(string probe, StringArray* list, int n)
{
    int i, j, len;
    int count, match=0;
  
    len = probe.length();
    for (i = 0; i < len; i++) 
    {
        count = 0;
        for (j = 0; j < n; j++) 
        { 
            if (probe[i] == (*list)[j][i]) 
            {
                match = j;
                count++;
            }
        }
        if (count == 0)
        { 
            return((int)-1);
        }
        if (count == 1) 
        {
            return(match);
        }
    }
    return((int)-1);

}

/*
 * The function getCmdLineDataStruct is used to return a struct with the values
 * specified. It returns it by value. I cant return by reference or else the object will
 * be destroyed.
 */
CmdLineData CommandLineParser::getCmdLineDataStruct(const char *str, int *flag, int type, 
                                                      StringArray* arg)
{
    CmdLineData tempStruct = {str, flag, type, arg};
    return tempStruct;
}

void CommandLineParser::printCmdLineData(const CmdLineData& temp)
{
    std::cout << "The str is: " << temp.str << std::endl;
    std::cout << "The int* is: " << *(temp.flag) << std::endl;
    std::cout << "The type is: " << temp.type << std::endl;
    std::cout << "The StringArray is: " << std::endl;
    
    if(temp.arg == NULL)
    {
        std::cout << "    NULL" << std::endl;
    }
    else
    {
        cout << "The number of elements is " << temp.arg->size() << std::endl;
        for(int i = 0; i < (int)temp.arg->size(); i++)
        {
            cout << "The " << i << "th element is: " << temp.arg->at(i) << endl;
        }
    }
}

/*
 * Helper function to change string to lower case.
 *
 */
string CommandLineParser::ConvertStringToLower(string strToConvert)
{
   for(unsigned int i=0;i<strToConvert.length();i++)
   {
      strToConvert[i] = tolower(strToConvert[i]);
   }
   return strToConvert;
}




void  CommandLineParser::exitWithErrorMsg(string msg)
{
    cerr << "ERROR: " << msg << std::endl;
    throw 1;
}

void CommandLineParser::reportBadOptionAndExit(string option, string expectedType)
{
    string msg;
    msg = "Bad option for ";
    msg += default_commandsep;// has to be added separately
    msg += option + ": expected " +  expectedType;
    exitWithErrorMsg(msg);
}

void CommandLineParser::reportInvalidOptionAndExit(string option)
{
    string msg = "Invalid option ";
    msg += default_commandsep;// has to be added separately
    msg += option;
    exitWithErrorMsg(msg);
}


}

