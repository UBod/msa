#include "RMuscle.h"
//#include "muscle.h"
#include "textfile.h"
//#include "seqvect.h"
#include "distfunc.h"
#include "msa.h"
#include "tree.h"
#include "profile.h"
#include "timing.h"
#include <vector>
#include <string>
using namespace std;

static char g_strUseTreeWarning[] =
"\n******** WARNING ****************\n"
"\nYou specified the -usetree option.\n"
"Note that a good evolutionary tree may NOT be a good\n"
"guide tree for multiple alignment. For more details,\n"
"please refer to the user guide. To disable this\n"
"warning, use -usetree_nowarn <treefilename>.\n\n";

void DoMuscle(MuscleInput *msaInput, MuscleOutput *msaOutput)
	{
	SetOutputFileName(g_pstrOutFileName);
	SetInputFileName(g_pstrInFileName);

	SetMaxIters(g_uMaxIters);
	SetSeqWeightMethod(g_SeqWeight1);

	if (msaInput->inputSeqs.Length() == 0) {
		TextFile fileIn(g_pstrInFileName);
		msaInput->inputSeqs.FromFASTAFile(fileIn);
	}

	const unsigned uSeqCount = msaInput->inputSeqs.Length();

	if (0 == uSeqCount)
		Quit("No sequences in input file");

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType)
		{
	case SEQTYPE_Auto:
		Alpha = msaInput->inputSeqs.GuessAlpha();
		break;

	case SEQTYPE_Protein:
		Alpha = ALPHA_Amino;
		break;

	case SEQTYPE_DNA:
		Alpha = ALPHA_DNA;
		break;

	case SEQTYPE_RNA:
		Alpha = ALPHA_RNA;
		break;

	default:
		Quit("Invalid seq type");
		}
	SetAlpha(Alpha);
	msaInput->inputSeqs.FixAlpha();

	PTR_SCOREMATRIX UserMatrix = 0;
	if (0 != g_pstrMatrixFileName)
		{
		const char *FileName = g_pstrMatrixFileName;
		const char *Path = getenv("MUSCLE_MXPATH");
		if (Path != 0)
			{
			size_t n = strlen(Path) + 1 + strlen(FileName) + 1;
			char *NewFileName = new char[n];
			sprintf(NewFileName, "%s/%s", Path, FileName);
			FileName = NewFileName;
			}
		TextFile File(FileName);
		UserMatrix = ReadMx(File);
		g_Alpha = ALPHA_Amino;
		g_PPScore = PPSCORE_SP;
	} else if (msaInput->hasSubstitutionMatrix) {
		UserMatrix = ReadMxFromR(msaInput->colNames, msaInput->substitutionMatrix);
		g_Alpha = ALPHA_Amino;
		g_PPScore = PPSCORE_SP;
	}

	SetPPScore();

	if (0 != UserMatrix)
		g_ptrScoreMatrix = UserMatrix;

	unsigned uMaxL = 0;
	unsigned uTotL = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned L = msaInput->inputSeqs.GetSeq(uSeqIndex).Length();
		uTotL += L;
		if (L > uMaxL)
			uMaxL = L;
		}

	SetIter(1);
	g_bDiags = g_bDiags1;
	SetSeqStats(uSeqCount, uMaxL, uTotL/uSeqCount);

	SetMuscleSeqVect(v);

	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msaInput->inputSeqs.SetSeqId(uSeqIndex, uSeqIndex);

	if (0 == uSeqCount)
		Quit("Input file '%s' has no sequences", g_pstrInFileName);
	if (1 == uSeqCount)
		{
		TextFile fileOut(g_pstrOutFileName, true);
		msaInput->inputSeqs.ToFile(fileOut);
		return;
		}

	if (uSeqCount > 1)
		MHackStart(msaInput->inputSeqs);

// First iteration
	Tree GuideTree;
	if (0 != g_pstrUseTreeFileName)
		{
	// Discourage users...
		if (!g_bUseTreeNoWarn)
			fprintf(stderr, "%s", g_strUseTreeWarning);

	// Read tree from file
		TextFile TreeFile(g_pstrUseTreeFileName);
		GuideTree.FromFile(TreeFile);

	// Make sure tree is rooted
		if (!GuideTree.IsRooted())
			Quit("User tree must be rooted");

		if (GuideTree.GetLeafCount() != uSeqCount)
			Quit("User tree does not match input sequences");

		const unsigned uNodeCount = GuideTree.GetNodeCount();
		for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
			{
			if (!GuideTree.IsLeaf(uNodeIndex))
				continue;
			const char *LeafName = GuideTree.GetLeafName(uNodeIndex);
			unsigned uSeqIndex;
			bool SeqFound = msaInput->inputSeqs.FindName(LeafName, &uSeqIndex);
			if (!SeqFound)
				Quit("Label %s in tree does not match sequences", LeafName);
			unsigned uId = msaInput->inputSeqs.GetSeqIdFromName(LeafName);
			GuideTree.SetLeafId(uNodeIndex, uId);
			}
		}
	else
		TreeFromSeqVect(msaInput->inputSeqs, GuideTree, g_Cluster1, g_Distance1, g_Root1,
		  g_pstrDistMxFileName1);

	const char *Tree1 = ValueOpt("Tree1");
	if (0 != Tree1)
		{
		TextFile f(Tree1, true);
		GuideTree.ToFile(f);
		if (g_bClusterOnly)
			GuideTree.Clear(); //valgrind
			//GuideTree = 0;
			return;
		}

	SetMuscleTree(GuideTree);
	ValidateMuscleIds(GuideTree);

	MSA msa;
	ProgNode *ProgNodes = 0;
	int PNlength = (msaInput->inputSeqs).Length();
	if (g_bLow)
		ProgNodes = ProgressiveAlignE(msaInput->inputSeqs, GuideTree, msa);
	else
		ProgressiveAlign(msaInput->inputSeqs, GuideTree, msa);
	SetCurrentAlignment(msa);

	if (0 != g_pstrComputeWeightsFileName)
		{
		extern void OutWeights(const char *FileName, const MSA &msa);
		SetMSAWeightsMuscle(msa);
		OutWeights(g_pstrComputeWeightsFileName, msa);
		GuideTree.Clear();
		return;
		}

	ValidateMuscleIds(msa);

	if (1 == g_uMaxIters || 2 == uSeqCount)
		{
		//TextFile fileOut(g_pstrOutFileName, true);
		//MHackEnd(msa);
		//msa.ToFile(fileOut);
		msaOutput->msa = DoMuscleOutput(msa);
		GuideTree.Clear();
		return;
		}

	if (0 == g_pstrUseTreeFileName)
		{
		g_bDiags = g_bDiags2;
		SetIter(2);

		if (g_bLow)
			{
			if (0 != g_uMaxTreeRefineIters)
				RefineTreeE(msa, msaInput->inputSeqs, GuideTree, ProgNodes);
			}
		else
			RefineTree(msa, GuideTree);

		const char *Tree2 = ValueOpt("Tree2");
		if (0 != Tree2)
			{
			TextFile f(Tree2, true);
			GuideTree.ToFile(f);
			}
		}

	SetSeqWeightMethod(g_SeqWeight2);
	SetMuscleTree(GuideTree);

	if (g_bAnchors)
		RefineVert(msa, GuideTree, g_uMaxIters - 2);
	else
		RefineHoriz(msa, GuideTree, g_uMaxIters - 2, false, false);

#if	0
// Refining by subfamilies is disabled as it didn't give better
// results. I tried doing this before and after RefineHoriz.
// Should get back to this as it seems like this should work.
	RefineSubfams(msa, GuideTree, g_uMaxIters - 2);
#endif

	ValidateMuscleIds(msa);
	ValidateMuscleIds(GuideTree);

	//TextFile fileOut(g_pstrOutFileName, true);
	//MHackEnd(msa);
	//msa.ToFile(fileOut);
	msaOutput->msa = DoMuscleOutput(msa);

	//GuideTree = 0;
	extern void DeleteProgNode(ProgNode &Node);
	const unsigned uNodeCount = GuideTree.GetNodeCount();
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
			DeleteProgNode(ProgNodes[uNodeIndex]); //valgrind

	delete[] ProgNodes;

	GuideTree.Clear(); //valgrind
	return;
}

void Run(MuscleInput *msaInput, MuscleOutput *msaOutput)
	{
	SetStartTime();
	Log("Started %s\n", GetTimeAsStr());
	for (int i = 0; i < g_argc; ++i)
		Log("%s ", g_argv[i]);
	Log("\n");

#if	TIMING
	TICKS t1 = GetClockTicks();
#endif
	if (g_bRefine)
		Refine();
	else if (g_bRefineW)
		{
		extern void DoRefineW();
		DoRefineW();
		}
	else if (g_bProfDB)
		ProfDB();
	else if (g_bSW)
		Local();
	else if (0 != g_pstrSPFileName)
		DoSP();
	else if (g_bProfile)
		Profile();
	else if (g_bPPScore)
		PPScore();
	else if (g_bPAS)
		ProgAlignSubFams();
	else if (g_bMakeTree)
		{
		extern void DoMakeTree();
		DoMakeTree();
		}
	else
		DoMuscle(msaInput, msaOutput);

#if	TIMING
	extern TICKS g_ticksDP;
	extern TICKS g_ticksObjScore;
	TICKS t2 = GetClockTicks();
	TICKS TotalTicks = t2 - t1;
	TICKS ticksOther = TotalTicks - g_ticksDP - g_ticksObjScore;
	double dSecs = TicksToSecs(TotalTicks);
	double PctDP = (double) g_ticksDP*100.0/(double) TotalTicks;
	double PctOS = (double) g_ticksObjScore*100.0/(double) TotalTicks;
	double PctOther = (double) ticksOther*100.0/(double) TotalTicks;
	Log("                 Ticks     Secs    Pct\n");
	Log("          ============  =======  =====\n");
	Log("DP        %12ld  %7.2f  %5.1f%%\n",
	  (long) g_ticksDP, TicksToSecs(g_ticksDP), PctDP);
	Log("OS        %12ld  %7.2f  %5.1f%%\n",
	  (long) g_ticksObjScore, TicksToSecs(g_ticksObjScore), PctOS);
	Log("Other     %12ld  %7.2f  %5.1f%%\n",
	  (long) ticksOther, TicksToSecs(ticksOther), PctOther);
	Log("Total     %12ld  %7.2f  100.0%%\n", (long) TotalTicks, dSecs);
#endif

	ListDiagSavings();
	Log("Finished %s\n", GetTimeAsStr());

	}
