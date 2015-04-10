/*********************************************************************
 * Clustal Omega - Multiple sequence alignment
 *
 * Copyright (C) 2010 University College Dublin
 *
 * Clustal-Omega is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This file is part of Clustal-Omega.
 *
 ********************************************************************/

/*
 *  RCS $Id: hhalign_wrapper.h 266 2011-10-20 16:04:11Z andreas $
 */

extern void
SetDefaultHhalignPara(hhalign_para *prHhalignPara);


extern double
HHalignWrapper(mseq_t *mseq, int *piOrderLR,
               double *pdSeqWeights, int iNodeCount,
               hmm_light *prHMMList, int iHMMCount,
               int iProfProfSeparator, hhalign_para rHhalignPara);
void
SanitiseUnknown(mseq_t *mseq);
