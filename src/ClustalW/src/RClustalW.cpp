#include "RClustalW.h"
#include "RClustalWMain.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

using namespace std;
using namespace Rcpp;

namespace clustalw
{
	RClustalWMain* rClustalWMain;
}
using namespace clustalw;

bool hasClustalWEntry(Rcpp::List params, const char* entry) {
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

bool fileExists (const char* name) {
    if (FILE *file = fopen(name, "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

SEXP RClustalW(SEXP rInputSeqs,
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
        struct ClustalWInput input;
        vector<string> args;
        args.push_back(".");
        Rcpp::List rparam(rParams); // Get parameters in params.

        //filename, if no file, then use "internalRsequence"
        //string file = "-INFILE=demoDNA.fasta"; //sequence from file
        bool inputFlag = false;
        if (hasClustalWEntry(rparam, "inputSeqIsFileFlag")) {
            inputFlag = as<bool>(rparam["inputSeqIsFileFlag"]);
        }

        if  (inputFlag) {
            string inputSeqFile = as<string>(rInputSeqs);
            //sequence from file
            string file = "-INFILE=" + inputSeqFile;
            args.push_back(file);
        } else {
            //sequence from R
            string file = "-INFILE=internalRsequence";
            args.push_back(file);

            //rInputSeqs
            vector<string> inputSeqs = Rcpp::as<vector<string> >(rInputSeqs);
            CharacterVector inputSeqsCV(rInputSeqs);
            vector<string> seqNames =
                    Rcpp::as<vector<string> >(inputSeqsCV.attr("names"));

            input.inputSeqs = inputSeqs;
            input.seqNames = seqNames;
        }

        //cluster
        if (!Rf_isNull(rCluster)) {
            string cluster = as<string>(rCluster);
            string clustering = "-CLUSTERING=" + cluster;
            args.push_back(clustering);
        }

        //gapOpening
        if (!Rf_isNull(rGapOpening)) {
            float gapOpening2 = as<float>(rGapOpening);
            stringstream gapOpening1;
            gapOpening1<<gapOpening2;
            string gapOpening = gapOpening1.str();
            string gapOpen = " GAPOPEN=" + gapOpening;
            args.push_back(std::string(gapOpen));
        }

        //gapExtension
        if (!Rf_isNull(rGapExtension)) {
            float gapExtension2 = as<float>(rGapExtension);
            stringstream gapExtension1;
            gapExtension1<<gapExtension2;
            string gapExtension = gapExtension1.str();
            string gapExt = " GAPEXT=" + gapExtension;
            args.push_back(std::string(gapExt));
        }

        //maxiters
        if (!Rf_isNull(rMaxiters)) {
            int maxiters2 = as<int>(rMaxiters);
            stringstream maxiters1;
            maxiters1<<maxiters2;
            string maxiters = maxiters1.str();
            string numiter = " NUMITER=" + maxiters;
            args.push_back(std::string(numiter));
        }

        //type
        if (!Rf_isNull(rType)) {
            string type = as<string>(rType);
            string autoType = "auto";
            //clustalW don't know type == "rna"
            if (autoType.compare(type) != 0 && type != "rna") {
                args.push_back(" TYPE=" + type);
            }
        }

        //rSubstitutionMatrix
        bool defaultFlag = false;
        if (hasClustalWEntry(rparam, "substitutionMatrixIsDefaultFlag")) {
            defaultFlag = as<bool>(rparam["substitutionMatrixIsDefaultFlag"]);
        }
        bool stringFlag = false;
        if (hasClustalWEntry(rparam, "substitutionMatrixIsStringFlag")) {
            stringFlag = as<bool>(rparam["substitutionMatrixIsStringFlag"]);
        }

        if (defaultFlag) {
            if (hasClustalWEntry(rparam, "pwdnamatrix")) {
		string pwDnaMat = as<string>(rparam["pwdnamatrix"]);
		Rprintf("use DNA substitution matrix: %s\n", pwDnaMat.c_str());
	    }
	    else {
		Rprintf("use default substitution matrix\n");
	    }
        } else if (stringFlag) {
            /*FIXME TODO change hardcoded rows and numbers*/
            string inputMatrixFile = as<string>(rSubstitutionMatrix);
            Rprintf("use matrix from file: %s\n", inputMatrixFile.c_str());
            string matrixFile = "-MATRIX=" + inputMatrixFile; //matrix from file
            args.push_back(matrixFile);
        } else {
            /*NumericMatrix substitutionMatrix(rSubstitutionMatrix);*/
            Rprintf("use user defined matrix\n");
            //FIXME TODO enable user matrix
            NumericMatrix substitutionMatrix(rSubstitutionMatrix);
            string dummyMatrix = "-MATRIX=dummyR.matrix";
            args.push_back(dummyMatrix);
            input.substitutionMatrix = substitutionMatrix;
        }

        //verbose
        bool verbose = false;
        if (!Rf_isNull(rVerbose)) {
            verbose = as<bool>(rVerbose);
            if (!verbose) {
                args.push_back("-QUIET");
            }
        }

        //all parameters in params

        //params$align
        if (hasClustalWEntry(rparam, "align")) {
            bool align = as<bool>(rparam["align"]);
            if (align) {
                args.push_back(" ALIGN");
            }
        }

        //params$bootlabels
        if (hasClustalWEntry(rparam, "bootlabels")) {
            string bootlabels = as<string>(rparam["bootlabels"]);
            args.push_back(" BOOTLABELS=" + bootlabels);
        }

        //params$bootstrap and params$bootstrapNo
        //DEACTIVATED, will be realized in later version
        /*if (hasClustalWEntry(rparam, "bootstrap")) {
            if (!hasClustalWEntry(rparam, "bootstrapNo")) {
                bool bootstrap = as<bool>(rparam["bootstrap"]);
                if (bootstrap) {
                    args.push_back(" BOOTSTRAP");
                }
            } else {
                //bootstrapNo can only be set if bootstrap==TRUE
                int bootstrapNo2 = as<int>(rparam["bootstrapNo"]);
                stringstream bootstrapNo1;
                bootstrapNo1<<bootstrapNo2;
                string bootstrapNo = bootstrapNo1.str();
                bool bootstrap = as<bool>(rparam["bootstrap"]);
                if (bootstrap) {
                    args.push_back(" BOOTSTRAP=" + bootstrapNo);
                }
            }
        }*/

        //params$case
        if (hasClustalWEntry(rparam, "case")) {
            string case1 = as<string>(rparam["case"]);
            args.push_back(" CASE=" + case1);
        }

        //params$check
        if (hasClustalWEntry(rparam, "check")) {
            bool check = as<bool>(rparam["check"]);
            if (check) {
                args.push_back(" CHECK:");
            }
        }

        //params$convert
        if (hasClustalWEntry(rparam, "convert")) {
            bool convert = as<bool>(rparam["convert"]);
            if (convert) {
                args.push_back(" CONVERT");
            }
        }

        //params$endgaps
        if (hasClustalWEntry(rparam, "endgaps")) {
            bool endgaps = as<bool>(rparam["endgaps"]);
            if (endgaps) {
                args.push_back(" ENDGAPS");
            }
        }

        //params$fullhelp
        if (hasClustalWEntry(rparam, "fullhelp")) {
            bool fullhelp = as<bool>(rparam["fullhelp"]);
            if (fullhelp) {
                args.push_back(" FULLHELP");
            }
        }

        //params$gapdist
        if (hasClustalWEntry(rparam, "gapdist")) {
            int gapdist2 = as<int>(rparam["gapdist"]);
            stringstream gapdist1;
            gapdist1<<gapdist2;
            string gapdist = gapdist1.str();
            args.push_back(std::string(" GAPDIST=" + gapdist));
        }

        //params$helixendin
        if (hasClustalWEntry(rparam, "helixendin")) {
            int helixendin2 = as<int>(rparam["helixendin"]);
            stringstream helixendin1;
            helixendin1<<helixendin2;
            string helixendin = helixendin1.str();
            args.push_back(std::string(" HELIXENDIN=" + helixendin));
        }

        //params$helixendout
        if (hasClustalWEntry(rparam, "helixendout")) {
            int helixendout2 = as<int>(rparam["helixendout"]);
            stringstream helixendout1;
            helixendout1<<helixendout2;
            string helixendout = helixendout1.str();
            args.push_back(std::string(" HELIXENDOUT=" + helixendout));
        }

        //params$helixgap
        if (hasClustalWEntry(rparam, "helixgap")) {
            int helixgap2 = as<int>(rparam["helixgap"]);
            stringstream helixgap1;
            helixgap1<<helixgap2;
            string helixgap = helixgap1.str();
            args.push_back(std::string(" HELIXGAP=" + helixgap));
        }


        //params$help
        if (hasClustalWEntry(rparam, "help")) {
            bool help = as<bool>(rparam["help"]);
            if (help) {
                args.push_back(" HELP");
            }
        }

        //params$hgapresidues
        if (hasClustalWEntry(rparam, "hgapresidues")) {
            string hgapresidues = as<string>(rparam["hgapresidues"]);
            args.push_back(" HGAPRESIDUES=" + hgapresidues);
        }

        //params$iteration
        if (hasClustalWEntry(rparam, "iteration")) {
            string iteration = as<string>(rparam["iteration"]);
            args.push_back(" ITERATION=" + iteration);
        }

        //params$kimura
        if (hasClustalWEntry(rparam, "kimura")) {
            bool kimura = as<bool>(rparam["kimura"]);
            if (kimura) {
                args.push_back(" KIMURA");
            }
        }

        //params$ktuple
        if (hasClustalWEntry(rparam, "ktuple")) {
            int ktuple2 = as<int>(rparam["ktuple"]);
            stringstream ktuple1;
            ktuple1<<ktuple2;
            string ktuple = ktuple1.str();
            args.push_back(std::string(" KTUPLE=" + ktuple));
        }

        //params$loopgap
        if (hasClustalWEntry(rparam, "loopgap")) {
            int loopgap2 = as<int>(rparam["loopgap"]);
            stringstream loopgap1;
            loopgap1<<loopgap2;
            string loopgap = loopgap1.str();
            args.push_back(std::string(" LOOPGAP=" + loopgap));
        }

        //params$maxdiv
        if (hasClustalWEntry(rparam, "maxdiv")) {
            int maxdiv2 = as<int>(rparam["maxdiv"]);
            stringstream maxdiv1;
            maxdiv1<<maxdiv2;
            string maxdiv = maxdiv1.str();
            args.push_back(std::string(" MAXDIV=" + maxdiv));
        }

        /* DEACTIVATED
        //params$maxseqlen
        if (hasClustalWEntry(rparam, "maxseqlen")) {
            int maxseqlen2 = as<int>(rparam["maxseqlen"]);
            stringstream maxseqlen1;
            maxseqlen1<<maxseqlen2;
            string maxseqlen = maxseqlen1.str();
            args.push_back(std::string(" MAXSEQLEN=" + maxseqlen));
        }*/

        //params$negative
        if (hasClustalWEntry(rparam, "negative")) {
            bool negative = as<bool>(rparam["negative"]);
            if (negative) {
                args.push_back(" NEGATIVE");
            }
        }

        //params$nohgap
        if (hasClustalWEntry(rparam, "nohgap")) {
            bool nohgap = as<bool>(rparam["nohgap"]);
            if (nohgap) {
                args.push_back(" NOHGAP");
            }
        }

        //params$nopgap
        if (hasClustalWEntry(rparam, "nopgap")) {
            bool nopgap = as<bool>(rparam["nopgap"]);
            if (nopgap) {
                args.push_back(" NOPGAP");
            }
        }

        //params$nosecstr1
        if (hasClustalWEntry(rparam, "nosecstr1")) {
            bool nosecstr1 = as<bool>(rparam["nosecstr1"]);
            if (nosecstr1) {
                args.push_back(" NOSECSTR1");
            }
        }

        //params$nosecstr2
        if (hasClustalWEntry(rparam, "nosecstr2")) {
            bool nosecstr2 = as<bool>(rparam["nosecstr2"]);
            if (nosecstr2) {
                args.push_back(" NOSECSTR2");
            }
        }

        //params$novgap
        if (hasClustalWEntry(rparam, "novgap")) {
            bool novgap = as<bool>(rparam["novgap"]);
            if (novgap) {
                args.push_back(" NOVGAP");
            }
        }

        //params$noweights
        if (hasClustalWEntry(rparam, "noweights")) {
            bool noweights = as<bool>(rparam["noweights"]);
            if (noweights) {
                args.push_back(" NOWEIGHTS");
            }
        }

        //params$options
        if (hasClustalWEntry(rparam, "options")) {
            bool options = as<bool>(rparam["options"]);
            if (options) {
                args.push_back(" OPTIONS");
            }
        }

        //params$outorder
        if (hasClustalWEntry(rparam, "outorder")) {
            string outorder = as<string>(rparam["outorder"]);
            args.push_back(" OUTORDER=" + outorder);
        }

        //params$output
        if (hasClustalWEntry(rparam, "output")) {
            string sOutput = as<string>(rparam["output"]);
            args.push_back(" OUTPUT=" + sOutput);
        }

        //params$outputtree
        if (hasClustalWEntry(rparam, "outputtree")) {
            string outputtree = as<string>(rparam["outputtree"]);
            args.push_back(" OUTPUTTREE=" + outputtree);
        }

        //params$pairgap
        if (hasClustalWEntry(rparam, "pairgap")) {
            int pairgap2 = as<int>(rparam["pairgap"]);
            stringstream pairgap1;
            pairgap1<<pairgap2;
            string pairgap = pairgap1.str();
            args.push_back(std::string(" PAIRGAP=" + pairgap));
        }

        //params$pim
        if (hasClustalWEntry(rparam, "pim")) {
            bool pim = as<bool>(rparam["pim"]);
            if (pim) {
                args.push_back(" PIM");
            }
        }

        //params$profile
        if (hasClustalWEntry(rparam, "profile")) {
            bool profile = as<bool>(rparam["profile"]);
            if (profile) {
                args.push_back(" PROFILE");
            }
        }

        //params$profile1
        if (hasClustalWEntry(rparam, "profile1")) {
            string profile1 = as<string>(rparam["profile1"]);
            args.push_back(" PROFILE1=" + profile1);
        }

        //params$profile2
        if (hasClustalWEntry(rparam, "profile2")) {
            string profile2 = as<string>(rparam["profile2"]);
            args.push_back(" PROFILE2=" + profile2);
        }

        //params$pwdnamatrix
        if (hasClustalWEntry(rparam, "pwdnamatrix")) {
            string pwdnamatrix = as<string>(rparam["pwdnamatrix"]);
            args.push_back(" PWDNAMATRIX=" + pwdnamatrix);
        }
	else if (hasClustalWEntry(rparam, "dnamatrix")) {
            string pwdnamatrix = as<string>(rparam["dnamatrix"]);
            args.push_back(" PWDNAMATRIX=" + pwdnamatrix);
	}

        //params$pwgapext
        if (hasClustalWEntry(rparam, "pwgapext")) {
            float pwgapext2 = as<float>(rparam["pwgapext"]);
            stringstream pwgapext1;
            pwgapext1<<pwgapext2;
            string pwgapext = pwgapext1.str();
            args.push_back(" PWGAPEXT=" + pwgapext);
        }

        //params$pwgapopen
        if (hasClustalWEntry(rparam, "pwgapopen")) {
            float pwgapopen2 = as<float>(rparam["pwgapopen"]);
            stringstream pwgapopen1;
            pwgapopen1<<pwgapopen2;
            string pwgapopen = pwgapopen1.str();
            args.push_back(std::string(" PWGAPOPEN=" + pwgapopen));
        }

        //params$pwmatrix
        if (hasClustalWEntry(rparam, "pwmatrix")) {
            string pwmatrix = as<string>(rparam["pwmatrix"]);
            args.push_back(" PWMATRIX=" + pwmatrix);
        }

        //params$quicktree
        if (hasClustalWEntry(rparam, "quicktree")) {
            bool quicktree = as<bool>(rparam["quicktree"]);
            if (quicktree) {
                args.push_back(" QUICKTREE");
            }
        }

        //params$range
        if (hasClustalWEntry(rparam, "range")) {
            std::vector<int> range = as<std::vector<int> >(rparam["range"]);
            int firstvalue2 = range.front();
            stringstream firstvalue1;
            firstvalue1<<firstvalue2;
            string firstvalue = firstvalue1.str();
            int secondvalue2 = range.back();
            stringstream secondvalue1;
            secondvalue1<<secondvalue2;
            string secondvalue = secondvalue1.str();
            string rangeVector= " RANGE=" + firstvalue + "," + secondvalue;
            if (!(firstvalue2 == -1 || secondvalue2 == -1)) {
                args.push_back(std::string(rangeVector));
            }
        }

        //params$score
        if (hasClustalWEntry(rparam, "score")) {
            string score = as<string>(rparam["score"]);
            args.push_back(" SCORE=" + score);
        }

        //params$secstrout
        if (hasClustalWEntry(rparam, "secstrout")) {
            string secstrout = as<string>(rparam["secstrout"]);
            args.push_back(" SECSTROUT=" + secstrout);
        }

        //params$seed
        if (hasClustalWEntry(rparam, "seed")) {
            int seed2 = as<int>(rparam["seed"]);
            stringstream seed1;
            seed1<<seed2;
            string seed = seed1.str();
            args.push_back(std::string(" SEED=" + seed));
        }

        //params$seqno_range
        if (hasClustalWEntry(rparam, "seqno_range")) {
            string seqno_range = as<string>(rparam["seqno_range"]);
            args.push_back(" SEQNO_RANGE=" + seqno_range);
        }

        //params$seqnos
        if (hasClustalWEntry(rparam, "seqnosFlag")) {
            bool isNull = as<bool>(rparam["seqnosFlag"]);
            if (!isNull) {
            //if (!Rf_isNull(rparam["seqnos"])){
                string seqnos = as<string>(rparam["seqnos"]);
                args.push_back(" SEQNOS=" + seqnos);
            }
        }

        //params$sequences
        if (hasClustalWEntry(rparam, "sequences")) {
            bool sequences = as<bool>(rparam["sequences"]);
            if (sequences) {
                args.push_back(" SEQUENCES");
            }
        }

        //params$strandendin
        if (hasClustalWEntry(rparam, "strandendin")) {
            int strandendin2 = as<int>(rparam["strandendin"]);
            stringstream strandendin1;
            strandendin1<<strandendin2;
            string strandendin = strandendin1.str();
            args.push_back(std::string(" STRANDENDIN=" + strandendin));
        }

        //params$strandendout
        if (hasClustalWEntry(rparam, "strandendout")) {
            int strandendout2 = as<int>(rparam["strandendout"]);
            stringstream strandendout1;
            strandendout1<<strandendout2;
            string strandendout = strandendout1.str();
            args.push_back(std::string(" STRANDENDOUT=" + strandendout));
        }

        //params$strandgap
        if (hasClustalWEntry(rparam, "strandgap")) {
            int strandgap2 = as<int>(rparam["strandgap"]);
            stringstream strandgap1;
            strandgap1<<strandgap2;
            string strandgap = strandgap1.str();
            args.push_back(std::string(" STRANDGAP=" + strandgap));
        }

        //params$terminalgap
        if (hasClustalWEntry(rparam, "terminalgap")) {
            int terminalgap2 = as<int>(rparam["terminalgap"]);
            stringstream terminalgap1;
            terminalgap1<<terminalgap2;
            string terminalgap = terminalgap1.str();
            args.push_back(std::string(" TERMINALGAP=" + terminalgap));
        }

        //params$topdiags
        if (hasClustalWEntry(rparam, "topdiags")) {
            int topdiags2 = as<int>(rparam["topdiags"]);
            stringstream topdiags1;
            topdiags1<<topdiags2;
            string topdiags = topdiags1.str();
            args.push_back(std::string(" TOPDIAGS=" + topdiags));
        }

        //params$tossgaps
        if (hasClustalWEntry(rparam, "tossgaps")) {
            bool tossgaps = as<bool>(rparam["tossgaps"]);
            if (tossgaps) {
                args.push_back(" TOSSGAPS");
            }
        }

        //params$transweight
        if (hasClustalWEntry(rparam, "transweight")) {
            float transweight2 = as<float>(rparam["transweight"]);
            stringstream transweight1;
            transweight1<<transweight2;
            string transweight = transweight1.str();
            args.push_back(std::string(" TRANSWEIGHT=" + transweight));
        }

        //params$tree
        if (hasClustalWEntry(rparam, "tree")) {
            bool tree = as<bool>(rparam["tree"]);
            if (tree) {
                args.push_back(" TREE");
            }
        }

        //params$window
        if (hasClustalWEntry(rparam, "window")) {
            int window2 = as<int>(rparam["window"]);
            stringstream window1;
            window1<<window2;
            string window = window1.str();
            args.push_back(std::string(" WINDOW=" + window));
        }

        if (verbose) {
            printf("params: ");
            for (int i = 0, n = args.size(); i < n; i++) {
                printf("%s", args[i].c_str());
                if (i < (n - 1)) {
                    printf(" ");
                }
            }
            printf("\n");
        }

        rClustalWMain = new RClustalWMain();

        ClustalWOutput output;
        rClustalWMain->run(args, &input, &output);

        delete(rClustalWMain);

        //return a more sophisticated object in later versions,
        //for now, we only return multiple sequence alignment
        retList = Rcpp::List::create(Rcpp::Named("msa") = Rcpp::CharacterVector(output.msa.begin(), output.msa.end()));

        if (fileExists("internalRsequence.aln")) {
        	remove("internalRsequence.aln");
        }
        if (fileExists("internalRsequence.dnd")) {
        	remove("internalRsequence.dnd");
        }

    } catch(int i) {
        if (i == 0) {
            Rprintf("ClustalW finished successfully");
        } else {
            Rf_error("ClustalW finished with errors");
        }
    } catch( std::exception &ex ) {
        forward_exception_to_r(ex);
    } catch(...) {
        Rf_error("ClustalW finished by an unknown reason");
    }
    return retList;
}
