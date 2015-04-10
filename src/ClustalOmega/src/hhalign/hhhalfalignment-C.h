/* -*- mode: c; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */

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
 *  RCS $Id: hhhalfalignment-C.h 227 2011-03-28 17:03:09Z fabian $
 */

// hhfullalignment.C

#ifndef MAIN
#define MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <ctype.h>    // islower, isdigit etc
using std::ios;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
#include "util-C.h"     // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.h"     // list data structure
#include "hash.h"     // hash data structure
#include "hhdecl-C.h"      // constants, class 
#include "hhutil-C.h"      // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hhhmm.h"       // class HMM
#include "hhalignment.h" // class Alignment
#include "hhhit.h"
#endif

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Methods of class HalfAlignment
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
  


/////////////////////////////////////////////////////////////////////////////
// Constructor
HalfAlignment::HalfAlignment(int maxseqdis)
{
  n=0; 
  sname=seq=NULL; 
  nss_dssp = nss_pred = nss_conf = nsa_dssp = ncons= -1;
  h = new(int[maxseqdis]);   //h[k] = next position of sequence k to be written
  s = new(char*[maxseqdis]);  //s[k][h] = character in column h, sequence k of output alignment
  l = new(int*[maxseqdis]);   //counts non-gap residues: l[k][i] = index of last residue AT OR BEFORE match state i in seq k
  m = new(int*[maxseqdis]);   //counts positions:        m[k][i] = position of match state i in string seq[k]  
}

/////////////////////////////////////////////////////////////////////////////////////
// Destructor
HalfAlignment::~HalfAlignment()
{
  Unset();
  delete[] h; h = NULL;
  delete[] s; s = NULL;
  delete[] l; l = NULL;
  delete[] m; m = NULL;
}



/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Free memory in HalfAlignment arrays s[][], l[][], and m[][]
 */
void 
HalfAlignment::Unset()
{
   // Free memory for alignment characters and residue counts
  for (int k=0; k<n; k++) 
    {
      delete[] s[k]; s[k] = NULL;
      delete[] l[k]; l[k] = NULL;
      delete[] m[k]; m[k] = NULL;
    }
  n=0; 
  sname=seq=NULL; 
  nss_dssp = nss_pred = nss_conf = nsa_dssp = ncons= -1;
}


//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Prepare a2m/a3m alignment: 
 * Calculate l[k][i] (residue indices) and m[k][i] (position in seq[k]) 
*/
void 
HalfAlignment::Set(char* name, char** seq_in, char** sname_in, int n_in, int L_in, int n1, int n2, int n3, int n4, int nc, int L_in2/*<--FS*/)
{
    int i;     /* counts match states in seq[k] */
    int ll;    /* counts residues LEFT from or at current position in seq[k] */
    int mm;    /* counts postions in string seq[k] */
    int k;     /* counts sequences */
    char c;
    char warned=0;
    
    nss_dssp=n1; nss_pred=n2; nss_conf=n3; nsa_dssp=n4; ncons=nc;
    seq=seq_in;     /* flat copy of sequences */
    sname=sname_in; /* flat copy of sequence names */
    n=n_in;
    L=L_in;    
    pos=0;

    /* Allocate memory for alignment characters and residue counts */
    for (k=0; k<n; k++)  {
        s[k]=new char[LINELEN];
        l[k]=new int[L+10+L_in2/*<--FS*/]; 
        m[k]=new int[L+10+L_in2/*<--FS*/];
        if (!s[k] || !l[k] || !m[k]) MemoryError("space for formatting HMM-HMM alignment");
        h[k]=0; //starting positions in alignment = 0
    } /* k <= 0 < n (= n_in) */

  for (k=0; k<n; k++) {
      m[k][0]=0;  // 0'th match state (virtual) is begin state at 0
      //if k is consensus sequence 
      if (k==nc) {
          for (i=1; i<=L; i++) m[k][i]=l[k][i]=i; 
          m[k][L+1]=l[k][L+1]=L; 
          continue;
      }
      i=1; mm=1; ll=1;
      while ((c=seq[k][mm]))
          {
              if (MatchChr(c)==c)    //count match/delete states
                  {
                      l[k][i]=ll;
                      m[k][i]=mm;
                      i++;
                  }
              if (WordChr(c)) ll++;  //index of next residue
              mm++;
          }
      l[k][i]=ll-1; //set l[k][L+1] eq number of residues in seq k (-1 since there is no residue at L+1st match state)
      m[k][i]=mm;   //set m[k][L+1]

      if ((i-1)!=L && !warned) 
          {
              cerr<<"Warning: sequence "<<sname[k]<<" in HMM "<<name<<" has "<<i<<" match states but should have "<<L<<"\n"; 
              warned=1;
          }
  } /* k <= 0 < n (= n_in) */
  //DEBUG
  if (v>=5)
      {
          printf("  i chr   m   l\n");
          for(i=0;i<=L+1;i++) printf("%3i   %1c %3i %3i\n",i,seq[0][m[0][i]],m[0][i],l[0][i]);
          printf("\n");
      }

} /*** end HalfAlignment::Set() ***/


/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Fill in insert states following match state i (without inserting '.' to fill up)
 */
void 
HalfAlignment::AddInserts(int i)
{
  for (int k=0; k<n; k++)                        // for all sequences...
    for (int mm=m[k][i]+1; mm<m[k][i+1]; mm++)   // for all inserts between match state i and i+1...
      s[k][h[k]++]=seq[k][mm];                   // fill inserts into output alignment s[k]
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Fill up alignment with gaps '.' to generate flush end (all h[k] equal)
 */
void 
HalfAlignment::FillUpGaps()
{ 
  int k;      //counts sequences
  pos=0;

  // Determine max position h[k]
  for (k=0; k<n; k++) pos = imax(h[k],pos);
  
  // Fill in gaps up to pos
  for (k=0; k<n; k++) 
    {
      for (int hh=h[k]; hh<pos; hh++) s[k][hh]='.';
      h[k]=pos;
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Fill in insert states following match state i and fill up gaps with '.' 
 */
void 
HalfAlignment::AddInsertsAndFillUpGaps(int i) 
{ 
  AddInserts(i); 
  FillUpGaps(); 
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Add gap column '.'
 */
void 
HalfAlignment::AddChar(char c)  
{ 
  for (int k=0; k<n; k++) s[k][h[k]++]=c;                
  pos++;
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Add match state column i as is
 */
void 
HalfAlignment::AddColumn(int i) 
{ 
  for (int k=0; k<n; k++) s[k][h[k]++]=seq[k][m[k][i]];  
  pos++;
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Add match state column i as insert state
 */
void 
HalfAlignment::AddColumnAsInsert(int i) 
{ 
  char c; 
  for (int k=0; k<n; k++) 
    if ((c=seq[k][m[k][i]])!='-' && (c<'0' || c>'9')) 
      s[k][h[k]++]=InsertChr(c); 
  pos++;
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Remove all characters c from template sequences
 */
void 
HalfAlignment::RemoveChars(char c)
{ 
  int k,h,hh;
  for (k=0; k<n; k++)
    {
      for (h=hh=0; h<pos; h++)
	if (s[k][h]!=c) s[k][hh++]=s[k][h];
      s[k][++hh]='\0';
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transform alignment sequences from A3M to A2M (insert ".")
 */
void 
HalfAlignment::BuildFASTA()
{
  AddInserts(0); 
  FillUpGaps();
  for (int i=1; i<=L; i++)
    {
      AddColumn(i); 
      AddInserts(i); 
      FillUpGaps();
    }
  ToFASTA();
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transform alignment sequences from A3M to A2M (insert ".")
 */
void 
HalfAlignment::BuildA2M()
{
  AddInserts(0); 
  FillUpGaps();
  for (int i=1; i<=L; i++)
    {
      AddColumn(i); 
      AddInserts(i); 
      FillUpGaps();
    }
  AddChar('\0');
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transform alignment sequences from A3M to A2M (insert ".")
 */
void 
HalfAlignment::BuildA3M()
{
  AddInserts(0); 
  for (int i=1; i<=L; i++)
    {
      AddColumn(i); 
      AddInserts(i); 
    }
  AddChar('\0');
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Transform alignment sequences from A2M to FASTA ( lowercase to uppercase and '.' to '-')
 */
void 
HalfAlignment::ToFASTA()
{
  for (int k=0; k<n; k++)
    {
      uprstr(s[k]);
      strtr(s[k],".","-");
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Align query (HalfAlignment) to template (i.e. hit) match state structure
 */
void 
HalfAlignment::AlignToTemplate(Hit& hit)
{
  int i,j;
  int step;    // column of the HMM-HMM alignment (first:nstep, last:1)
  char state;

  if(0) {  //par.loc==0) { //////////////////////////////////////////// STILL NEEDED??
    // If in global mode: Add part of alignment before first MM state
    AddInserts(0); // Fill in insert states before first match state
    for (i=1; i<hit.i[hit.nsteps]; i++)
      {
	AddColumnAsInsert(i);
	AddInserts(i);
	if (par.outformat<=2) FillUpGaps();
      }
  }

  // Add endgaps (First state must be an MM state!!)
  for (j=1; j<hit.j[hit.nsteps]; j++)    
    {
      AddChar('-');
    }

  // Add alignment between first and last MM state
  for (step=hit.nsteps; step>=1; step--) 
  {
    state = hit.states[step];
    i = hit.i[step];

    switch(state)
      {
      case MM:  //MM pair state (both query and template in Match state)
	AddColumn(i);
	AddInserts(i);
	break;
      case DG: //D- state
      case MI: //MI state
	AddColumnAsInsert(i);
	AddInserts(i);
	break;
      case GD: //-D state
      case IM: //IM state
	AddChar('-');
	break;
      }
    if (par.outformat<=2) FillUpGaps();

  }

  if(0) { //par.loc==0) { //////////////////////////////////////////// STILL NEEDED??

    // If in global mode: Add part of alignment after last MM state
    for (i=hit.i[1]+1; i<=L; i++)    
      {
	AddColumnAsInsert(i);
	AddInserts(i);
	if (par.outformat==2) FillUpGaps();
      }
  }

  // Add endgaps 
  for (j=hit.j[1]+1; j<=hit.L; j++)    
    {
      AddChar('-');
    }

  // Add end-of-string character
  AddChar('\0');
}


/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Write the a2m/a3m alignment into alnfile 
 */
void 
HalfAlignment::Print(char* alnfile)
{
  int k;      //counts sequences
  int omitted=0; // counts number of sequences with no residues in match states
  FILE *outf;
  if (strcmp(alnfile,"stdout"))
    {
      if (par.append) outf=fopen(alnfile,"a"); else outf=fopen(alnfile,"w");
      if (!outf) OpenFileError(alnfile);
    } 
  else
    outf = stdout;
  if (v>=3) cout<<"Writing alignment to "<<alnfile<<"\n";

  for (k=0; k<n; k++)
    {
      // Print sequence only if it contains at least one residue in a match state
      if (1) //strpbrk(s[k],"ABCDEFGHIKLMNPQRSTUVWXYZ1234567890")) 
	{
	  fprintf(outf,">%s\n",sname[k]);
	  fprintf(outf,"%s\n",s[k]);
	} else {
	  omitted++;
	  if (v>=3) printf("%-14.14s contains no residue in match state. Omitting sequence\n",sname[k]);
	}
    }
  if (v>=2 && omitted) printf("Omitted %i sequences in %s which contained no residue in match state\n",omitted,alnfile);
  fclose(outf);
}


/** EOF hhhalfalignment-C.h **/
