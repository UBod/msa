#include "muscle.h"
#include "pwpath.h"
#include "timing.h"
#include "textfile.h"
#include "msa.h"
#include "profile.h"

SCORE GlobalAlign(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	if (g_bDiags)
		return GlobalAlignDiags(PA, uLengthA, PB, uLengthB, Path);
	else
		return GlobalAlignNoDiags(PA, uLengthA, PB, uLengthB, Path);
	}

SCORE GlobalAlignNoDiags(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	switch (g_PPScore)
		{
	case PPSCORE_LE:
		return GlobalAlignLA(PA, uLengthA, PB, uLengthB, Path);

	case PPSCORE_SP:
		return GlobalAlignNS(PA, uLengthA, PB, uLengthB, Path);

	case PPSCORE_SV:
		return GlobalAlignSimple(PA, uLengthA, PB, uLengthB, Path);
		}
	return 0;
	}
