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
#include "interface/CommandLineParser.h"
#include "general/DebugLog.h"
#include "general/ClustalWResources.h"
#include "general/Stats.h"
#include <ctime>

#include <vector>
#include "RClustalWMain.h"
using namespace std;
using namespace clustalw;
using namespace Rcpp;

namespace clustalw
{
    UserParameters* userParameters;
    SubMatrix *subMatrix;
    Utility* utilityObject;
    DebugLog* logObject;
    Stats* statsObject;

RClustalWMain::RClustalWMain() {}

RClustalWMain::~RClustalWMain() {}

void RClustalWMain::run(UserArgs args, ClustalWInput *input, ClustalWOutput *output) {

	userParameters = new UserParameters(false);
	utilityObject = new Utility();
	subMatrix = new SubMatrix();
	statsObject = new Stats();
	ClustalWResources *resources = ClustalWResources::Instance();
	resources->setPathToExecutable(args.at(0));
	userParameters->setDisplayInfo(true);



	//userParameters->setDebug(5);
	#if DEBUGFULL
		if(DEBUGLOG) {
			cout << "debugging is on\n\n\n";
			logObject = new DebugLog("logfile.txt");
			logObject->logMsg("Loggin is on!");
		}
	#endif

	if (args.size() > 1) {
		//time_t start, end;
		//double dif;
		//start = time (NULL);
		//userParameters->setDisplayInfo(false);
		CommandLineParser cmdLineParser(&args, false);
		cmdLineParser.run(&args, false, input, output);
		//for (int i = 0; i < result.size(); i++) {
		//	utilityObject->info("[%s]", result[i].c_str());
		//}

		//if (statsObject->isEnabled())
		//	statsObject->logCmdLine(argc,argv);

		//end = time (NULL);
		//dif = difftime(end, start);
		//cout << "It took " << dif << " seconds\n";
	}

	delete userParameters;
	delete utilityObject;
	delete subMatrix;
	delete statsObject;

	if(logObject)
	{
		delete logObject;
	}
	return;
}
}
