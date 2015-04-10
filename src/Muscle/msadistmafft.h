#ifndef MSADistMAFFT_h
#define MSADistMAFFT_h

#include "msadist.h"
#include <math.h>

extern double PctIdToMAFFTDist(double dPctId);

class MSADistMAFFT : public MSADist
	{
public:
	virtual double ComputeDist(const MSA &msa, unsigned uSeqIndex1,
	  unsigned uSeqIndex2)
		{
		double dPctId = msa.GetPctIdentityPair(uSeqIndex1, uSeqIndex2);
		//if (dPctId < 0.05)
		//	dPctId = 0.05;
		//double dDist = -log(dPctId);
		//return dDist;
		return PctIdToMAFFTDist(dPctId);
		}
	};

#endif	// MSADistMAFFT_h
