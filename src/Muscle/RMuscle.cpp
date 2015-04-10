#include "RMuscle.h"
using namespace std;
#include <stdio.h>
#include <stdlib.h>
//#include "seqvect.h"
#include <R.h>
#include <Rinternals.h>
#include "seq.h"
#ifdef	WIN32
//#include <windows.h>	// for SetPriorityClass()
#include <io.h>			// for isatty()
#include <sstream>
#else
#include <unistd.h>		// for isatty()
#endif

using namespace Rcpp;

bool hasMuscleEntry(Rcpp::List params, const char* entry) {
	Rcpp::CharacterVector namesCV = params.names();
	int n = namesCV.size();
	vector<string> namesVector = Rcpp::as< std::vector< std::string > >(namesCV);
	for (int i = 0; i < n; i++) {
		const char* vec = namesVector[i].c_str();
		if (strcmp(vec, entry) == 0) {
			return !Rf_isNull(params[entry]);
		}
	}
	return false;
}

SEXP RMuscle(SEXP rInputSeqs,
             SEXP rCluster,
             SEXP rGapOpening,
             SEXP rGapExtension,
             SEXP rMaxiters,
             SEXP rSubstitutionMatrix,
             SEXP rType,
             SEXP rVerbose,
             SEXP rParams) {

	Rcpp::List retList;
	try {
		Rcpp::List params(rParams); // Get parameters in params.

		bool inputFlag = false;
		if (hasMuscleEntry(params, "inputSeqIsFileFlag")) {
			inputFlag = as<bool>(params["inputSeqIsFileFlag"]);
		}

		stringstream ss;
		MuscleInput input;
		Seq* seq = new Seq(); //valgrind

		if (inputFlag) {
			string inputFile = as<string>(rInputSeqs);
			ss << "-in " << inputFile;
		} else {
			ss << "-in noFile";

			//inputSeq, seqNames
			CharacterVector inputSeqs(rInputSeqs);
			vector<string> seqNames = as<vector<string> >(inputSeqs.attr("names"));

			for (int i = 0, n = inputSeqs.size(); i < n; i++) {
				seq->FromString(inputSeqs[i], seqNames[i].c_str()); //FIXME if string is too long, the string is omitted
				input.inputSeqs.AppendSeq(*seq);
			}
		}
		//cluster
		if (!Rf_isNull(rCluster)) {
			string cluster = as<string>(rCluster);
			ss << " -cluster " << cluster;
		}

		//gapOpen
		if (!Rf_isNull(rGapOpening)) {
			double gapOpening = as<double>(rGapOpening);
			ss << " -gapOpen " << gapOpening;
		}


		//gapExtend
		if (!Rf_isNull(rGapExtension)) {
			double gapExtension = as<double>(rGapExtension);
			ss << " -gapExtend " << gapExtension;
		}


		//maxiters
		if (!Rf_isNull(rMaxiters)) {
			int maxIters = as<int>(rMaxiters);
			ss << " -maxiters " << maxIters;
		}

		//substitutionMatrix TODO --> domuscle->DoMuscle --> ReadMX
		//vector<int> substitutionMatrix =
		//                 Rcpp::as<vector<int>(rSubstitutionMatrix);
		if (!Rf_isNull(rSubstitutionMatrix)) {
			NumericMatrix substitutionMatrix(rSubstitutionMatrix);
			int nrows = substitutionMatrix.nrow();
			int ncolumns = substitutionMatrix.ncol();
			for (int i = 0; i < ncolumns; i++) {
				for (int j = 0; j < nrows; j++) {
					input.substitutionMatrix[i][j] = substitutionMatrix(i,j);
				}
			}

			Function colnamesFun("colnames");
			input.colNames = as<vector<string> >(colnamesFun(substitutionMatrix));
			input.hasSubstitutionMatrix = true;
		} else {
			input.hasSubstitutionMatrix = false;
		}

		//type
		if (!Rf_isNull(rType)) {
			string type = as<string>(rType);
			ss << " -seqtype " << type;
		}


		//verbose
		bool verbose = false;
		if (!Rf_isNull(rVerbose)) {
			verbose = as<bool>(rVerbose);
			if (!verbose) {
				ss << " -quiet";
			}
		}

		//params

		if (hasMuscleEntry(params, "anchorspacing")) {
			int anchorspacing = as<int>(params["anchorspacing"]);
			ss << " -anchorspacing " << anchorspacing;
		}

		if (hasMuscleEntry(params, "center")) {
			double center = as<double>(params["center"]);
			ss << " -center " << center;
		}

		if (hasMuscleEntry(params, "cluster1")) {
			string cluster1 = as<string>(params["cluster1"]);
			ss << " -cluster1 " << cluster1;
		}

		if (hasMuscleEntry(params, "cluster2")) {
			string cluster2 = as<string>(params["cluster2"]);
			ss << " -cluster2 " << cluster2;
		}

		if (hasMuscleEntry(params, "diagbreak")) {
			int diagbreak = as<int>(params["diagbreak"]);
			ss << " -diagbreak " << diagbreak;
		}

		if (hasMuscleEntry(params, "diaglength")) {
			int diaglength = as<int>(params["diaglength"]);
			ss << " -diaglength " << diaglength;
		}

		if (hasMuscleEntry(params, "diagmargin")) {
			int diagmargin = as<int>(params["diagmargin"]);
			ss << " -diagmargin " << diagmargin;
		}

		if (hasMuscleEntry(params, "distance1")) {
			string distance1 = as<string>(params["distance1"]);
			ss << " -distance1 " << distance1;
		}

		if (hasMuscleEntry(params, "distance2")) {
			string distance2 = as<string>(params["distance2"]);
			ss << " -distance2 " << distance2;
		}

		if (hasMuscleEntry(params, "hydro")) {
			int hydro = as<int>(params["hydro"]);
			ss << " -hydro " << hydro;
		}

		if (hasMuscleEntry(params, "hydrofactor")) {
			double hydrofactor = as<double>(params["hydrofactor"]);
			ss << " -hydrofactor " << hydrofactor;
		}

		if (hasMuscleEntry(params, "in1")) {
			string in1 = as<string>(params["in1"]);
			ss << " -in1 " << in1;
		}

		if (hasMuscleEntry(params, "in2")) {
			string in2 = as<string>(params["in2"]);
			ss << " -in2 " << in2;
		}

		if (hasMuscleEntry(params, "maxhours")) {
			double maxhours = as<double>(params["maxhours"]);
			if (maxhours != -1) {
				ss << " -maxhours " << maxhours;
			}
		}

		if (hasMuscleEntry(params, "maxtrees")) {
			int maxtrees = as<int>(params["maxtrees"]);
			ss << " -maxtrees " << maxtrees;
		}

		if (hasMuscleEntry(params, "minbestcolscore")) {
			double minbestcolscore = as<double>(params["minbestcolscore"]);
			ss << " -minbestcolscore " << minbestcolscore;
		}

		if (hasMuscleEntry(params, "minsmoothscore")) {
			double minsmoothscore = as<double>(params["minsmoothscore"]);
			ss << " -minsmoothscore " << minsmoothscore;
		}

		if (hasMuscleEntry(params, "objscore")) {
			string objscore = as<string>(params["objscore"]);
			ss << " -objscore " << objscore;
		}

		if (hasMuscleEntry(params, "refinewindow")) {
			int refinewindow = as<int>(params["refinewindow"]);
			ss << " -refinewindow " << refinewindow;
		}

		if (hasMuscleEntry(params, "root1")) {
			string root1 = as<string>(params["root1"]);
			ss << " -root1 " << root1;
		}

		if (hasMuscleEntry(params, "root2")) {
			string root2 = as<string>(params["root2"]);
			ss << " -root2 " << root2;
		}

		if (hasMuscleEntry(params, "smoothscoreceil")) {
			double smoothscoreceil = as<double>(params["smoothscoreceil"]);
			ss << " -smoothscoreceil " << smoothscoreceil;
		}

		if (hasMuscleEntry(params, "smoothwindow")) {
			int smoothwindow = as<int>(params["smoothwindow"]);
			ss << " -smoothwindow " << smoothwindow;
		}

		if (hasMuscleEntry(params, "SUEFF")) {
			double SUEFF = as<double>(params["SUEFF"]);
			ss << " -SUEFF " << SUEFF;
		}

		if (hasMuscleEntry(params, "weight1")) {
			string weight1 = as<string>(params["weight1"]);
			ss << " -weight1 " << weight1;
		}

		if (hasMuscleEntry(params, "weight2")) {
			string weight2 = as<string>(params["weight2"]);
			ss << " -weight2 " << weight2;
		}

		if (hasMuscleEntry(params, "anchors")) {
			bool anchors = as<bool>(params["anchors"]);
			if (anchors) {
				ss << " -anchors ";
			}
		}

		if (hasMuscleEntry(params, "brenner")) {
			bool brenner = as<bool>(params["brenner"]);
			if (brenner) {
				ss << " -brenner ";
			}
		}

		if (hasMuscleEntry(params, "core")) {
			bool core = as<bool>(params["core"]);
			if (core) {
				ss << " -core ";
			}
		}

		if (hasMuscleEntry(params, "diags")) {
			bool diags = as<bool>(params["diags"]);
			if (diags) {
				ss << " -diags ";
			}
		}

		if (hasMuscleEntry(params, "diags1")) {
			bool diags1 = as<bool>(params["diags1"]);
			if (diags1) {
				ss << " -diags1 ";
			}
		}

		if (hasMuscleEntry(params, "diags2")) {
			bool diags2 = as<bool>(params["diags2"]);
			if (diags2) {
				ss << " -diags2 ";
			}
		}

		if (hasMuscleEntry(params, "dimer")) {
			bool dimer = as<bool>(params["dimer"]);
			if (dimer) {
				ss << " -dimer ";
			}
		}

		/*
		if (hasMuscleEntry(params, "group")) {
			bool group = as<bool>(params["group"]);
			//Rprintf("Group: %s\n", params["group"] ? "True" : "False");
			if (group) {
				ss << " -group ";
			}
		}
		*/

		if (hasMuscleEntry(params, "le")) {
			bool le = as<bool>(params["le"]);
			if (le) {
				ss << " -le ";
			}
		}

		if (hasMuscleEntry(params, "noanchors")) {
			bool noanchors = as<bool>(params["noanchors"]);
			if (noanchors) {
				ss << " -noanchors ";
			}
		}

		if (hasMuscleEntry(params, "nocore")) {
			bool nocore = as<bool>(params["nocore"]);
			if (nocore) {
				ss << " -nocore ";
			}
		}

        if (hasMuscleEntry(params, "profile")) {
            bool profile = as<bool>(params["profile"]);
            if (profile) {
                ss << " -profile ";
            }
        }

        if (hasMuscleEntry(params, "refine")) {
            bool refine = as<bool>(params["refine"]);
            if (refine) {
                ss << " -refine ";
            }
        }

		if (hasMuscleEntry(params, "refinew")) {
			bool refinew = as<bool>(params["refinew"]);
			if (refinew) {
				ss << " -refinew ";
			}
		}

		if (hasMuscleEntry(params, "sp")) {
			bool sp = as<bool>(params["sp"]);
			if (sp) {
				ss << " -sp ";
			}
		}

		if (hasMuscleEntry(params, "spn")) {
			bool spn = as<bool>(params["spn"]);
			if (spn) {
				ss << " -spn ";
			}
		}

		if (hasMuscleEntry(params, "spscore")) {
			bool spscore = as<bool>(params["spscore"]);
			if (spscore) {
				ss << " -spscore ";
			}
		}

		/*
		if (hasMuscleEntry(params, "stable")) {
			bool stable = as<bool>(params["stable"]);
			//Rprintf("Stable: %s\n", params["stable"] ? "True" : "False");
			if (stable) {
				ss << " -stable ";
			}
		}*/

		if (hasMuscleEntry(params, "sv")) {
			bool sv = as<bool>(params["sv"]);
			if (sv) {
				ss << " -sv ";
			}
		}

		/*
		if (hasMuscleEntry(params, "termgaps4")) {
			bool termgaps4 = as<bool>(params["termgaps4"]);
			//Rprintf("Termgaps4: %s\n", params["termgaps4"] ? "True" : "False");
			if (termgaps4) {
				ss << " -termgaps4 ";
			}
		}*/

		/*
		if (hasMuscleEntry(params, "termgapsfull")) {
			bool termgapsfull = as<bool>(params["termgapsfull"]);
			//Rprintf("Termgapsfull: %s\n", params["termgapsfull"] ? "True" : "False");
			if (termgapsfull) {
				ss << " -termgapsfull ";
			}
		}*/


		/*
		if (hasMuscleEntry(params, "termgapshalf")) {
			bool termgapshalf = as<bool>(params["termgapshalf"]);
			Rprintf("Termgapshalf: %s\n", params["termgapshalf"] ? "True" : "False");
			if (termgapshalf) {
				ss << " -termgapshalf ";
			}
		}*/


		/*
		if (hasMuscleEntry(params, "termgapshalflonger")) {
			bool termgapshalflonger = as<bool>(params["termgapshalflonger"]);
			//Rprintf("Termgapshalflonger: %s\n", params["termgapshalflonger"] ? "True" : "False");
			if (termgapshalflonger) {
				ss << " -termgapshalflonger ";
			}
		}*/

		if (hasMuscleEntry(params, "version")) {
			bool version = as<bool>(params["version"]);
			if (version) {
				ss << " -version ";
			}
		}

		//static params
		ss << " -clwstrict ";

		string s = ss.str();
		if (verbose) {
			Rprintf("params: %s\n", s.c_str());
		}
	
	
		//#if	WIN32
		// Multi-tasking does not work well in CPU-bound
		// console apps running under Win32.
		// Reducing the process priority allows GUI apps
		// to run responsively in parallel.
		//SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
		//#endif

		SetNewHandler();
		SetStartTime();
		ProcessArgStr(s.c_str());
		SetParams();
		SetLogFile();

		if (g_bVersion) {
			printf("%s\n", MUSCLE_LONG_VERSION);
		}

		if (!g_bQuiet)
			Credits();
	
		if (MissingCommand() && isatty(0)) {
			Usage();
		}

		MuscleOutput output;
	
		if (g_bCatchExceptions) {
			try	{
				Run(&input, &output);
			} catch (...) {
				OnException();
			}
		} else {
			Run(&input, &output);
		}

		//CleanupNewHandler(); //valgrind
		delete(seq);

		retList = Rcpp::List::create(Rcpp::Named("msa") = Rcpp::CharacterVector(output.msa.begin(), output.msa.end()));

	} catch(int i) {
		if (i == 0) {
			//Rprintf("MUSCLE finished successfully");
		} else {
			Rf_error("MUSCLE finished with errors");
		}
	} catch( std::exception &ex ) {
		forward_exception_to_r(ex);
	} catch(...) {
		Rf_error("MUSCLE finished by an unknown reason");
	}
	return retList;
}
