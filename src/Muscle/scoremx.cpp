#include "muscle.h"
#include "profile.h"

extern SCOREMATRIX VTML_LA;
extern SCOREMATRIX PAM200;
extern SCOREMATRIX PAM200NoCenter;
extern SCOREMATRIX VTML_SP;
extern SCOREMATRIX VTML_SPNoCenter;
extern SCOREMATRIX NUC_SP;

PTR_SCOREMATRIX g_ptrScoreMatrix;

void SetScoreMatrix()
	{
	switch (g_PPScore)
		{
	case PPSCORE_LE:
		g_ptrScoreMatrix = &VTML_LA;
		break;

	case PPSCORE_SP:
		if (g_bPrecompiledCenter)
			g_ptrScoreMatrix = &PAM200;
		else
			g_ptrScoreMatrix = &PAM200NoCenter;
		break;

	case PPSCORE_SV:
		if (g_bPrecompiledCenter)
			g_ptrScoreMatrix = &VTML_SP;
		else
			g_ptrScoreMatrix = &VTML_SPNoCenter;
		break;

	case PPSCORE_SPN:
		if (g_bPrecompiledCenter)
			g_ptrScoreMatrix = &NUC_SP;
		else
			Quit("SPN requires precompiled center");
		break;

	default:
		Quit("Invalid g_PPScore");
		}
	}
