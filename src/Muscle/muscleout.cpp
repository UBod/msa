#include "muscle.h"
#include "msa.h"
#include "params.h"
#include "textfile.h"
#include <vector>
#include <string>
using namespace std;

static vector<string> DoOutput(MSA &msa)
	{
	bool AnyOutput = false;
	vector<string> result;

// Value options
	if (g_pstrFASTAOutFileName)
		{
		TextFile File(g_pstrFASTAOutFileName, true);
		msa.ToFASTAFile(File);
		AnyOutput = true;
		}

	if (g_pstrMSFOutFileName)
		{
		TextFile File(g_pstrMSFOutFileName, true);
		msa.ToMSFFile(File);
		AnyOutput = true;
		}

	if (g_pstrClwOutFileName)
		{
		TextFile File(g_pstrClwOutFileName, true);
		result = msa.ToAlnFile();
		AnyOutput = true;
		}
	if (g_pstrClwStrictOutFileName)
		{
		g_bClwStrict = true;
		TextFile File(g_pstrClwStrictOutFileName, true);
		result = msa.ToAlnFile();
		AnyOutput = true;
		}

	if (g_pstrHTMLOutFileName)
		{
		TextFile File(g_pstrHTMLOutFileName, true);
		msa.ToHTMLFile(File);
		AnyOutput = true;
		}

	if (g_pstrPHYIOutFileName)
		{
		TextFile File(g_pstrPHYIOutFileName, true);
		msa.ToPhyInterleavedFile(File);
		AnyOutput = true;
		}

	if (g_pstrPHYSOutFileName)
		{
		TextFile File(g_pstrPHYSOutFileName, true);
		msa.ToPhySequentialFile(File);
		AnyOutput = true;
		}

// Flag options, at most one used (because only one -out filename)
	TextFile fileOut(g_pstrOutFileName, true);
	if (g_bFASTA)
		{
		msa.ToFASTAFile(fileOut);
		AnyOutput = true;
		}
	else if (g_bMSF)
		{
		msa.ToMSFFile(fileOut);
		AnyOutput = true;
		}
	else if (g_bAln)
		{
		result = msa.ToAlnFile();
		AnyOutput = true;
		}
	else if (g_bHTML)
		{
		msa.ToHTMLFile(fileOut);
		AnyOutput = true;
		}
	else if (g_bPHYI)
		{
		msa.ToPhyInterleavedFile(fileOut);
		AnyOutput = true;
		}
	else if (g_bPHYS)
		{
		msa.ToPhySequentialFile(fileOut);
		AnyOutput = true;
		}

// If -out option was given but no flags, output as FASTA
	if (!AnyOutput) {
		msa.ToFASTAFile(fileOut);
	}
	//fileOut.Close();

	if (0 != g_pstrScoreFileName) {
		WriteScoreFile(msa);
	}
	return result;
}

vector<string> DoMuscleOutput(MSA &msa)
	{
	MHackEnd(msa);
	if (g_bStable)
		{
		MSA msaStable;
		Stabilize(msa, msaStable);
		msa.Clear();	// save memory
		return DoOutput(msaStable);
		}
	return DoOutput(msa);
	}
