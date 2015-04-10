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
 *  RCS $Id: mymain.h 97 2010-07-12 15:30:26Z andreas $
 */

typedef struct {
	int seqLength;
	char **inputSeqs;
	char **seqNames;
	char substitutionMatrix;
} ClustalOmegaInput;

typedef struct { //can be use for further result objects
	char **msa_v; //multiple sequence alignment
	int msa_c;
} ClustalOmegaOutput;


extern int
executeClustalOmega(int argc, char **argv, ClustalOmegaInput *msaInput, ClustalOmegaOutput *msaOutput);

