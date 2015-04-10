#include "muscle.h"
#include "msa.h"

/***
Compute Henikoff weights.
Steven Henikoff and Jorja G. Henikoff (1994), Position-based sequence weights.
J. Mol. Biol., 243(4):574-578.

Award each different residue an equal share of the weight, and then to divide up
that weight equally among the sequences sharing the same residue. So if in a
position of a multiple alignment, r different residues are represented, a residue
represented in only one sequence contributes a score of 1/r to that sequence, whereas a
residue represented in s sequences contributes a score of 1/rs to each of the s
sequences. For each sequence, the contributions from each position are summed to give
a sequence weight.

See also HenikoffWeightPB.
***/

void MSA::CalcHenikoffWeightsCol(unsigned uColIndex) const
	{
	const unsigned uSeqCount = GetSeqCount();

// Compute letter counts in this column
	unsigned uLetterCount[MAX_ALPHA];
	memset(uLetterCount, 0, sizeof(uLetterCount));
	unsigned uDifferentLetterCount = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned uLetter = GetLetterEx(uSeqIndex, uColIndex);
		if (uLetter >= 20)
			continue;
		unsigned uNewCount = uLetterCount[uLetter] + 1;
		uLetterCount[uLetter] = uNewCount;
		if (1 == uNewCount)
			++uDifferentLetterCount;
		}

// Compute weight contributions
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned uLetter = GetLetterEx(uSeqIndex, uColIndex);
		if (uLetter >= 20)
			continue;
		const unsigned uCount = uLetterCount[uLetter];
		unsigned uDenom = uCount*uDifferentLetterCount;
		if (uDenom == 0)
			continue;
		m_Weights[uSeqIndex] += (WEIGHT) (1.0/uDenom);
		}
	}

void MSA::SetHenikoffWeights() const
	{
	const unsigned uColCount = GetColCount();
	const unsigned uSeqCount = GetSeqCount();

	if (0 == uSeqCount)
		return;
	else if (1 == uSeqCount)
		{
		m_Weights[0] = (WEIGHT) 1.0;
		return;
		}
	else if (2 == uSeqCount)
		{
		m_Weights[0] = (WEIGHT) 0.5;
		m_Weights[1] = (WEIGHT) 0.5;
		return;
		}

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		m_Weights[uSeqIndex] = 0.0;

	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		CalcHenikoffWeightsCol(uColIndex);

// Set all-gap seqs weight to 0
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		if (IsGapSeq(uSeqIndex))
			m_Weights[uSeqIndex] = 0.0;

	Normalize(m_Weights, uSeqCount);
	}
