#include "muscle.h"
#include "msa.h"
#include "seqvect.h"
#include "textfile.h"

#define MEMDEBUG	0

#if	MEMDEBUG
#include <crtdbg.h>
#endif

void MUSCLE(SeqVect &v, MSA &msaOut);

// Append msa2 at the end of msa1
void AppendMSA(MSA &msa1, const MSA &msa2)
	{
	const unsigned uSeqCount = msa1.GetSeqCount();

	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();

	const unsigned uColCountCat = uColCount1 + uColCount2;

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned uId = msa1.GetSeqId(uSeqIndex);
		unsigned uSeqIndex2;
		bool bFound = msa2.GetSeqIndex(uId, &uSeqIndex2);
		if (bFound)
			{
			for (unsigned uColIndex = 0; uColIndex < uColCount2; ++uColIndex)
				{
				const char c = msa2.GetChar(uSeqIndex2, uColIndex);
				msa1.SetChar(uSeqIndex, uColCount1 + uColIndex, c);
				}
			}
		else
			{
			for (unsigned uColIndex = 0; uColIndex < uColCount2; ++uColIndex)
				msa1.SetChar(uSeqIndex, uColCount1 + uColIndex, '-');
			}
		}
	}

static void SeqFromMSACols(const MSA &msa, unsigned uSeqIndex, unsigned uColFrom,
  unsigned uColTo, Seq &s)
	{
	s.Clear();
	s.SetName(msa.GetSeqName(uSeqIndex));
	s.SetId(msa.GetSeqId(uSeqIndex));
	for (unsigned uColIndex = uColFrom; uColIndex <= uColTo; ++uColIndex)
		{
		char c = msa.GetChar(uSeqIndex, uColIndex);
		if (!IsGapChar(c))
			s.AppendChar(c);
		}
	}

static void SeqVectFromMSACols(const MSA &msa, unsigned uColFrom, unsigned uColTo,
  SeqVect &v)
	{
	v.Clear();
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq s;
		SeqFromMSACols(msa, uSeqIndex, uColFrom, uColTo, s);
		v.AppendSeq(s);
		}
	}

void RefineW(const MSA &msaIn, MSA &msaOut)
	{
	const unsigned uSeqCount = msaIn.GetSeqCount();
	const unsigned uColCount = msaIn.GetColCount();

// Reserve same nr seqs, 20% more cols
	const unsigned uReserveColCount = (uColCount*120)/100;
	msaOut.SetSize(uSeqCount, uReserveColCount);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		msaOut.SetSeqName(uSeqIndex, msaIn.GetSeqName(uSeqIndex));
		msaOut.SetSeqId(uSeqIndex, msaIn.GetSeqId(uSeqIndex));
		}

	const unsigned uWindowCount = (uColCount + g_uRefineWindow - 1)/g_uRefineWindow;
	if (0 == g_uWindowTo)
		g_uWindowTo = uWindowCount - 1;

#if	MEMDEBUG
	_CrtSetBreakAlloc(1560);
#endif

	if (g_uWindowOffset > 0)
		{
		MSA msaTmp;
		MSAFromColRange(msaIn, 0, g_uWindowOffset, msaOut);
		}

	fprintf(stderr, "\n");
	for (unsigned uWindowIndex = g_uWindowFrom; uWindowIndex <= g_uWindowTo; ++uWindowIndex)
		{
		fprintf(stderr, "Window %d of %d    \r", uWindowIndex, uWindowCount);
		const unsigned uColFrom = g_uWindowOffset + uWindowIndex*g_uRefineWindow;
		unsigned uColTo = uColFrom + g_uRefineWindow - 1;
		if (uColTo >= uColCount)
			uColTo = uColCount - 1;
		assert(uColTo >= uColFrom);

		SeqVect v;
		SeqVectFromMSACols(msaIn, uColFrom, uColTo, v);

#if	MEMDEBUG
		_CrtMemState s1;
		_CrtMemCheckpoint(&s1);
#endif

		MSA msaTmp;
		MUSCLE(v, msaTmp);
		AppendMSA(msaOut, msaTmp);
		if (uWindowIndex == g_uSaveWindow)
			{
			MSA msaInTmp;
			unsigned uOutCols = msaOut.GetColCount();
			unsigned un = uColTo - uColFrom + 1;
			MSAFromColRange(msaIn, uColFrom, un, msaInTmp);

			char fn[256];
			sprintf(fn, "win%d_inaln.tmp", uWindowIndex);
			TextFile fIn(fn, true);
			msaInTmp.ToFile(fIn);

			sprintf(fn, "win%d_inseqs.tmp", uWindowIndex);
			TextFile fv(fn, true);
			v.ToFile(fv);

			sprintf(fn, "win%d_outaln.tmp", uWindowIndex);
			TextFile fOut(fn, true);
			msaTmp.ToFile(fOut);
			}

#if	MEMDEBUG
		void FreeDPMemSPN();
		FreeDPMemSPN();

		_CrtMemState s2;
		_CrtMemCheckpoint(&s2);

		_CrtMemState s;
		_CrtMemDifference(&s, &s1, &s2);

		_CrtMemDumpStatistics(&s);
		_CrtMemDumpAllObjectsSince(&s1);
		throw 1;
#endif
//#if	DEBUG
//		AssertMSAEqIgnoreCaseAndGaps(msaInTmp, msaTmp);
//#endif
		}
	fprintf(stderr, "\n");

//	AssertMSAEqIgnoreCaseAndGaps(msaIn, msaOut);//@@uncomment!
	}

void DoRefineW()
	{
	SetOutputFileName(g_pstrOutFileName);
	SetInputFileName(g_pstrInFileName);
	SetStartTime();

	SetMaxIters(g_uMaxIters);
	SetSeqWeightMethod(g_SeqWeight1);

	TextFile fileIn(g_pstrInFileName);
	MSA msa;
	msa.FromFile(fileIn);

	const unsigned uSeqCount = msa.GetSeqCount();
	if (0 == uSeqCount)
		Quit("No sequences in input file");

	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msa.SetSeqId(uSeqIndex, uSeqIndex);
	SetMuscleInputMSA(msa);

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType)
		{
	case SEQTYPE_Auto:
		Alpha = msa.GuessAlpha();
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
		Quit("Invalid SeqType");
		}
	SetAlpha(Alpha);
	msa.FixAlpha();

	if (ALPHA_DNA == Alpha || ALPHA_RNA == Alpha)
		SetPPScore(PPSCORE_SPN);

	MSA msaOut;
	RefineW(msa, msaOut);

//	ValidateMuscleIds(msa);

//	TextFile fileOut(g_pstrOutFileName, true);
//	msaOut.ToFile(fileOut);
	DoMuscleOutput(msaOut);
	}
