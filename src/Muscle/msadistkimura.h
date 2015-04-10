#ifndef MSADistKimura_h
#define MSADistKimura_h

#include "msadist.h"

class MSADistKimura : public MSADist
	{
public:
	virtual double ComputeDist(const MSA &msa, unsigned uSeqIndex1,
	  unsigned uSeqIndex2);
	};

#endif	// MSADistKimura_h
