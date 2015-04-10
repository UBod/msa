/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <iostream>
#include "alignment/Alignment.h"
#include "alignment/Sequence.h"
#include "general/clustalw.h"
#include "general/UserParameters.h"
#include "substitutionMatrix/SubMatrix.h"
#include "general/Utility.h"
#include "fileInput/FileReader.h"
// #include "interface/InteractiveMenu.h"
#include "interface/CommandLineParser.h"
#include "general/DebugLog.h"
#include "general/ClustalWResources.h"
#include "general/Stats.h"
#include <ctime>
namespace clustalw
{ 
    UserParameters* userParameters;
    Utility* utilityObject;
    SubMatrix *subMatrix; 
    DebugLog* logObject;
    Stats* statsObject;
}
using namespace std;
using namespace clustalw;

int main(int argc, char **argv)
{      
    userParameters = new UserParameters(false);
    utilityObject = new Utility();   
    subMatrix = new SubMatrix();
    statsObject = new Stats();
    ClustalWResources *resources = ClustalWResources::Instance();
    resources->setPathToExecutable(string(argv[0]));
    userParameters->setDisplayInfo(true);


    
    //userParameters->setDebug(5);       
    #if DEBUGFULL    
        if(DEBUGLOG)
        {
            cout << "debugging is on\n\n\n";
            logObject = new DebugLog("logfile.txt");
            logObject->logMsg("Loggin is on!");
        }
    #endif

    if (argc > 1)
    {    
        //time_t start, end;
        //double dif;
        //start = time (NULL);
        //userParameters->setDisplayInfo(false);        
        vector<string> args;
        for (int i = 1; i < argc; ++i)
        {
            args.push_back(argv[i]);
        }
        CommandLineParser cmdLineParser(&args, false);
        
        if (statsObject->isEnabled())
            statsObject->logCmdLine(argc,argv);
        
        //end = time (NULL);
        //dif = difftime(end, start);
        //cout << "It took " << dif << " seconds\n";
    }
    if (argc<=1 || userParameters->getInteractive())
    {
        // FIXME: additional parameters like infile are ignored!
       	//InteractiveMenu menu;
        userParameters->setMenuFlag(true);
        userParameters->setInteractive(true);
 	//menu.mainMenu();
    }
    delete userParameters;
    delete utilityObject;
    delete subMatrix;
    
    if(logObject)
    {
        delete logObject;
    }
    return 0;
}

