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
 * RCS $Id: hhhit-C.h 284 2013-06-12 10:10:11Z fabian $
 */

// hhhit.C

#ifndef MAIN
#define MAIN
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <stdio.h>    // printf
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <time.h>     // clock
#include <ctype.h>    // islower, isdigit etc
#include "util-C.h"     // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "list.h"     // list data structure
#include "hash.h"     // hash data structure
#include "hhdecl-C.h"      // constants, class 
#include "hhutil-C.h"      // imax, fmax, iround, iceil, ifloor, strint, strscn, strcut, substr, uprstr, uprchr, Basename etc.
#include "hhhmm.h"       // class HMM
#include "hhalignment.h" // class Alignment
#include "hhhitlist.h"   // class HitList
#endif

#define CALCULATE_MAX6(max, var1, var2, var3, var4, var5, var6, varb) \
if (var1>var2) { max=var1; varb=STOP;} \
else           { max=var2; varb=MM;}; \
if (var3>max)  { max=var3; varb=GD;}; \
if (var4>max)  { max=var4; varb=IM;}; \
if (var5>max)  { max=var5; varb=DG;}; \
if (var6>max)  { max=var6; varb=MI;}; 

#define CALCULATE_SUM6(sum, var1, var2, var3, var4, var5, var6, varb) \
if (var1>var2) { sum=var1; varb=STOP;} \
else           { sum=var2; varb=MM;}; \
if (var3>sum)  { sum=var3; varb=GD;}; \
if (var4>sum)  { sum=var4; varb=IM;}; \
if (var5>sum)  { sum=var5; varb=DG;}; \
if (var6>sum)  { sum=var6; varb=MI;}; \
sum = var1 + var2 + var3 + var4 + var5 + var6;

#define CALCULATE_MAX4(max, var1, var2, var3, var4, varb) \
if (var1>var2) { max=var1; varb=STOP;} \
else           { max=var2; varb=MM;}; \
if (var3>max)  { max=var3; varb=MI;}; \
if (var4>max)  { max=var4; varb=IM;}; 

// Generate random number in [0,1[
#define frand() ((float) rand()/(RAND_MAX+1.0))


// Function declarations
inline float Score(float* qi, float* tj);
inline float ProbFwd(float* qi, float* tj);
inline float max2(const float& xMM, const float& xX, char& b); 
inline int pickprob2(const double& xMM, const double& xX, const int& state); 
inline int pickprob3_GD(const double& xMM, const double& xDG, const double& xGD); 
inline int pickprob3_IM(const double& xMM, const double& xMI, const double& xIM); 
inline int pickprob6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI); 
inline int pickmax2(const double& xMM, const double& xX, const int& state); 
inline int pickmax3_GD(const double& xMM, const double& xDG, const double& xGD); 
inline int pickmax3_IM(const double& xMM, const double& xMI, const double& xIM); 
inline int pickmax6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI); 
inline double Pvalue(double x, double a[]);
inline double Pvalue(float x, float lamda, float mu);
inline double logPvalue(float x, float lamda, float mu);
inline double logPvalue(float x, double a[]);
inline double Probab(Hit& hit);

//////////////////////////////////////////////////////////////////////////////
//// Constructor
//////////////////////////////////////////////////////////////////////////////
Hit::Hit()
{
  longname = name = file = dbfile = NULL;
  sname = NULL;
  seq = NULL;  
  bMM = bGD = bDG = bIM = bMI = NULL;
  self = 0;
  i = j = NULL;
  states = NULL;
  S = S_ss = P_posterior = NULL;
  Xcons = NULL;
  B_MM=B_MI=B_IM=B_DG=B_GD=NULL;
  F_MM=F_MI=F_IM=F_DG=F_GD=NULL;
  cell_off = NULL;
  scale = NULL;
  sum_of_probs=0.0; 
  Neff_HMM=0.0;
  score_ss = Pval = logPval = Eval = Probab = Pforward = 0.0;
  nss_conf = nfirst = i1 = i2 = j1 = j2 = matched_cols = ssm1 = ssm2 = 0;
}

//////////////////////////////////////////////////////////////////////////////
/**
 * @brief Free all allocated memory (to delete list of hits)
 */
void 
Hit::Delete()
{
  if (i){
    delete[] i; i = NULL;
  }
  if (j){
    delete[] j; j = NULL;
  }
  if (states){
    delete[] states; states = NULL;
  }
  if (S){
    delete[] S; S = NULL;
  }
  if (S_ss){
    delete[] S_ss; S_ss = NULL;
  }
  if (P_posterior){
    delete[] P_posterior; P_posterior = NULL;
  }
  if (Xcons){
    delete[] Xcons; Xcons = NULL;
  }
  //  delete[] l;    l = NULL;
  i = j = NULL;
  states = NULL;
  S = S_ss = P_posterior = NULL;
  Xcons = NULL;
  if (irep==1) // if irep>1 then longname etc point to the same memory locations as the first repeat. 
    {          // but these have already been deleted.
      // 	printf("Delete name = %s\n",name);//////////////////
      delete[] longname; longname = NULL;
      delete[] name; name = NULL;
      delete[] file; file = NULL;
      delete[] dbfile; dbfile = NULL;
      for (int k=0; k<n_display; k++) 
          {
              delete[] sname[k]; sname[k] = NULL;
              //delete[] seq[k]; seq[k] = NULL;  //disabled: MSAR-172
          }
      delete[] sname; sname = NULL;
      delete[] seq; seq = NULL;
    }
}


///////////////////////////////////////////////////////////////////////////
/**
 * @brief Allocate/delete memory for dynamic programming matrix
 */
void 
Hit::AllocateBacktraceMatrix(int Nq, int Nt)
{
  int i;
  bMM=new(char*[Nq]);
  bMI=new(char*[Nq]);
  bIM=new(char*[Nq]);
  bDG=new(char*[Nq]);
  bGD=new(char*[Nq]);
  cell_off=new(char*[Nq]);
  for (i=0; i<Nq; i++) 
    {
      bMM[i]=new(char[Nt]);
      bMI[i]=new(char[Nt]);
      bIM[i]=new(char[Nt]);
      bGD[i]=new(char[Nt]);
      bDG[i]=new(char[Nt]);
      cell_off[i]=new(char[Nt]);
      if (!bMM[i] || !bMI[i] || !bIM[i] || !bGD[i] || !bDG[i] || !cell_off[i]) 
	{
	  fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",i+1,Nq);
	  fprintf(stderr,"Suggestions:\n");
	  fprintf(stderr,"1. Cut query sequence into shorter segments\n");
	  fprintf(stderr,"2. Check stack size limit (Linux: ulimit -a)\n");
	  fprintf(stderr,"3. Run on a computer with bigger memory\n");
	  throw 3;
	} 
    }
}

/**
 * @brief
 */
void 
Hit::DeleteBacktraceMatrix(int Nq)
{
  int i;

  if (NULL != bMM){ /* FS, r259 -> r260 */
      for (i=0; i<Nq; i++)  {
          delete[] bMM[i]; bMM[i] = NULL;
          delete[] bMI[i]; bMI[i] = NULL;
          delete[] bIM[i]; bIM[i] = NULL;
          delete[] bGD[i]; bGD[i] = NULL;
          delete[] bDG[i]; bDG[i] = NULL;
          delete[] cell_off[i]; cell_off[i] = NULL;
      }
      delete[] bMM; bMM = NULL;
      delete[] bMI; bMI = NULL;
      delete[] bIM; bIM = NULL;
      delete[] bDG; bDG = NULL;
      delete[] bGD; bGD = NULL;
      delete[] cell_off; cell_off = NULL;
  }
}


///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Allocate/delete memory for Forward dynamic programming matrix
 */
void 
Hit::AllocateForwardMatrix(int Nq, int Nt)
{
  F_MM=new(double*[Nq]);
  F_MI=new(double*[Nq]);
  F_DG=new(double*[Nq]);
  F_IM=new(double*[Nq]);
  F_GD=new(double*[Nq]);
  scale=new(double[Nq+1]); // need Nq+3?
  for (int i=0; i<Nq; i++) 
    {
      F_MM[i] = new(double[Nt]);
      F_MI[i] = new(double[Nt]);
      F_DG[i] = new(double[Nt]);
      F_IM[i] = new(double[Nt]);
      F_GD[i] = new(double[Nt]);
      if (!F_MM[i] || !F_MI[i] || !F_IM[i] || !F_GD[i] || !F_DG[i]) 
	{
	  fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",i+1,Nq);
	  fprintf(stderr,"Suggestions:\n");
	  fprintf(stderr,"1. Cut query sequence into shorter segments\n");
	  fprintf(stderr,"2. Check stack size limit (Linux: ulimit -a)\n");
	  fprintf(stderr,"3. Run on a computer with bigger memory\n");
	  throw 3;
	} 

    }
}

/**
 * @brief
 */
void 
Hit::DeleteForwardMatrix(int Nq)
{

    if (NULL != F_MM){ /* FS, r259 -> r260 */
        for (int i=0; i<Nq; i++) {
            delete[] F_MM[i]; F_MM[i] = NULL;
            delete[] F_MI[i]; F_MI[i] = NULL;
            delete[] F_IM[i]; F_IM[i] = NULL;
            delete[] F_GD[i]; F_GD[i] = NULL;
            delete[] F_DG[i]; F_DG[i] = NULL;
        }
        delete[] F_MM; F_MM = NULL;
        delete[] F_MI; F_MI = NULL;
        delete[] F_IM; F_IM = NULL;
        delete[] F_DG; F_DG = NULL;
        delete[] F_GD; F_GD = NULL;
        delete[] scale; scale = NULL;
    }
}

/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Allocate/delete memory for Backward dynamic programming matrix (DO ONLY AFTER FORWARD MATRIX HAS BEEN ALLOCATED)
 */
void 
Hit::AllocateBackwardMatrix(int Nq, int Nt)
{
  B_MM=new(double*[Nq]);
  B_MI=F_MI; 
  B_DG=F_DG; 
  B_IM=F_IM; 
  B_GD=F_GD; 
  for (int i=0; i<Nq; i++) 
    {
      B_MM[i] = new(double[Nt]);
      if (!B_MM[i]) 
	{
	  fprintf(stderr,"Error: out of memory while allocating row %i (out of %i) for dynamic programming matrices \n",i+1,Nq);
	  fprintf(stderr,"Suggestions:\n");
	  fprintf(stderr,"1. Cut query sequence into shorter segments\n");
	  fprintf(stderr,"2. Check stack size limit (Linux: ulimit -a)\n");
	  fprintf(stderr,"3. Run on a computer with bigger memory\n");
	  throw 3;
	} 
    }
}

void Hit::DeleteBackwardMatrix(int Nq)
{

    if (NULL != B_MM){ /* FS, r259 -> r260 */
        for (int i=0; i<Nq; i++) {
            delete[] B_MM[i]; B_MM[i] = NULL;  /* is this all? FS */
        }
        delete[] B_MM; B_MM = NULL;
        B_MM=B_MI=B_IM=B_DG=B_GD=NULL;
    }
}



/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compare HMMs with one another and look for sub-optimal alignments that share no pair with previous ones
 * The function is called with q and t
 * If q and t are equal (self==1), only the upper right part of the matrix is calculated: j>=i+3
 */
void 
Hit::Viterbi(HMM& q, HMM& t, float** Sstruc)
{
  
    // Linear topology of query (and template) HMM:
    // 1. The HMM HMM has L+2 columns. Columns 1 to L contain 
    //    a match state, a delete state and an insert state each.
    // 2. The Start state is M0, the virtual match state in column i=0 (j=0). (Therefore X[k][0]=ANY)
    //    This column has only a match state and it has only a transitions to the next match state.
    // 3. The End state is M(L+1), the virtual match state in column i=L+1.(j=L+1) (Therefore X[k][L+1]=ANY)
    //    Column L has no transitions to the delete state: tr[L][M2D]=tr[L][D2D]=0.
    // 4. Transitions I->D and D->I are ignored, since they do not appear in PsiBlast alignments 
    //    (as long as the gap opening penalty d is higher than the best match score S(a,b)). 
    
    // Pairwise alignment of two HMMs:
    // 1. Pair-states for the alignment of two HMMs are 
    //    MM (Q:Match T:Match) , GD (Q:Gap T:Delete), IM (Q:Insert T:Match),  DG (Q:Delelte, T:Match) , MI (Q:Match T:Insert) 
    // 2. Transitions are allowed only between the MM-state and each of the four other states.
    
    // Saving space:
    // The best score ending in pair state XY sXY[i][j] is calculated from left to right (j=1->t.L) 
    // and top to bottom (i=1->q.L). To save space, only the last row of scores calculated is kept in memory.
    // (The backtracing matrices are kept entirely in memory [O(t.L*q.L)]).
    // When the calculation has proceeded up to the point where the scores for cell (i,j) are caculated,
    //    sXY[i-1][j'] = sXY[j']   for j'>=j (A below)  
    //    sXY[i][j']   = sXY[j']   for j'<j  (B below)
    //    sXY[i-1][j-1]= sXY_i_1_j_1         (C below) 
    //    sXY[i][j]    = sXY_i_j             (D below)
    //                   j-1   
    //                     j
    // i-1:               CAAAAAAAAAAAAAAAAAA
    //  i :   BBBBBBBBBBBBBD
    
    
    // Variable declarations
    //float sMM[MAXRES];          // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match) 
    //float sGD[MAXRES];          // sGD[i][j] = score of best alignment up to indices (i,j) ending in (Gap,Delete) 
    //float sDG[MAXRES];          // sDG[i][j] = score of best alignment up to indices (i,j) ending in (Delete,Gap)
    //float sIM[MAXRES];          // sIM[i][j] = score of best alignment up to indices (i,j) ending in (Ins,Match)
    //float sMI[MAXRES];          // sMI[i][j] = score of best alignment up to indices (i,j) ending in (Match,Ins) 
    float *sMM = new(float[par.maxResLen]);   // sMM[i][j] = score of best alignment up to indices (i,j) ending in (Match,Match) 
    float *sGD = new(float[par.maxResLen]);   // sGD[i][j] = score of best alignment up to indices (i,j) ending in (Gap,Delete) 
    float *sDG = new(float[par.maxResLen]);   // sDG[i][j] = score of best alignment up to indices (i,j) ending in (Delete,Gap)
    float *sIM = new(float[par.maxResLen]);   // sIM[i][j] = score of best alignment up to indices (i,j) ending in (Ins,Match)
    float *sMI = new(float[par.maxResLen]);   // sMI[i][j] = score of best alignment up to indices (i,j) ending in (Match,Ins) 
    float smin=(par.loc? 0:-FLT_MAX);  //used to distinguish between SW and NW algorithms in maximization         
    int i=0,j=0;      //query and template match state indices
    float sMM_i_j=0, sMI_i_j=0, sIM_i_j=0, sGD_i_j=0, sDG_i_j=0;
    float sMM_i_1_j_1=0, sMI_i_1_j_1=0, sIM_i_1_j_1=0, sGD_i_1_j_1=0, sDG_i_1_j_1=0;
    int jmin=0, jmax=0;

  // Reset crossed out cells?
  if(irep==1) InitializeForAlignment(q,t);

  // Initialization of top row, i.e. cells (0,j)
  for (j=0; j<=t.L; j++) 
    {
      sMM[j] = (self? 0 : -j*par.egt);
      sIM[j] = sMI[j] = sDG[j] = sGD[j] = -FLT_MAX; 
    }
  score=-INT_MAX; i2=j2=0; bMM[0][0]=STOP;

  // Viterbi algorithm
  for (i=1; i<=q.L; i++) // Loop through query positions i
    {
//       if (v>=5) printf("\n");

      
      if (self) 
	{
	  // If q is compared to itself, ignore cells below diagonal+SELFEXCL
	  jmin = i+SELFEXCL; 
	  jmax = t.L;
	  if (jmin>jmax) continue;
	}
      else
	{
	  // If q is compared to t, exclude regions where overlap of q with t < min_overlap residues
	  jmin=imax( 1, i+min_overlap-q.L);  // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq} 
	  jmax=imin(t.L,i-min_overlap+t.L);  // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt} 
	}      

      // Initialize cells
      if (jmin==1) 
	{
	  sMM_i_1_j_1 = -(i-1)*par.egq;  // initialize at (i-1,0)
	  sMM[0] = -i*par.egq;           // initialize at (i,0)
	  sIM_i_1_j_1 = sMI_i_1_j_1 = sDG_i_1_j_1 = sGD_i_1_j_1 = -FLT_MAX; // initialize at (i-1,jmin-1)
	} 
      else 
	{
	  // Initialize at (i-1,jmin-1) if lower left triagonal is excluded due to min overlap
	  sMM_i_1_j_1 = sMM[jmin-1];     // initialize at (i-1,jmin-1)
	  sIM_i_1_j_1 = sIM[jmin-1];     // initialize at (i-1,jmin-1)
	  sMI_i_1_j_1 = sMI[jmin-1];     // initialize at (i-1,jmin-1)
	  sDG_i_1_j_1 = sDG[jmin-1];     // initialize at (i-1,jmin-1)
	  sGD_i_1_j_1 = sGD[jmin-1];     // initialize at (i-1,jmin-1)
	  sMM[jmin-1] = -FLT_MAX;        // initialize at (i,jmin-1)
	}
      if (jmax<t.L) // initialize at (i-1,jmmax) if upper right triagonal is excluded due to min overlap
	sMM[jmax] = sIM[jmax] = sMI[jmax] = sDG[jmax] = sGD[jmax] = -FLT_MAX; 
      sIM[jmin-1] = sMI[jmin-1] = sDG[jmin-1] = sGD[jmin-1] = -FLT_MAX; // initialize at (i,jmin-1)
      
      for (j=jmin; j<=jmax; j++) // Loop through template positions j
	{

	  if (cell_off[i][j])
	    {
	      sMM_i_1_j_1 = sMM[j]; // sMM_i_1_j_1 (for j->j+1) = sMM(i-1,(j+1)-1) = sMM[j] 
	      sGD_i_1_j_1 = sGD[j];
	      sIM_i_1_j_1 = sIM[j];
	      sDG_i_1_j_1 = sDG[j];
	      sMI_i_1_j_1 = sMI[j];
	      sMM[j]=sMI[j]=sIM[j]=sDG[j]=sGD[j]=-FLT_MAX; // sMM[j] = sMM(i,j) is cell_off
	    }
	  else 
	    {
	      // Recursion relations
// 	      printf("S[%i][%i]=%4.1f  ",i,j,Score(q.p[i],t.p[j])); // DEBUG!!

	      CALCULATE_MAX6( sMM_i_j,
			      smin,
			      sMM_i_1_j_1 + q.tr[i-1][M2M] + t.tr[j-1][M2M], 
			      sGD_i_1_j_1 + q.tr[i-1][M2M] + t.tr[j-1][D2M],
			      sIM_i_1_j_1 + q.tr[i-1][I2M] + t.tr[j-1][M2M],
			      sDG_i_1_j_1 + q.tr[i-1][D2M] + t.tr[j-1][M2M],
			      sMI_i_1_j_1 + q.tr[i-1][M2M] + t.tr[j-1][I2M],
			      bMM[i][j]
			      );
 	      sMM_i_j += Score(q.p[i],t.p[j]) + ScoreSS(q,t,i,j) + par.shift 
		+ (Sstruc==NULL? 0: Sstruc[i][j]); 
	      

	      sGD_i_j = max2
	              (
		       sMM[j-1] + t.tr[j-1][M2D], // MM->GD gap opening in query 
		       sGD[j-1] + t.tr[j-1][D2D], // GD->GD gap extension in query 
		       bGD[i][j]
		       );
	      sIM_i_j = max2
 	              (
// 		       sMM[j-1] + q.tr[i][M2I] + t.tr[j-1][M2M] ,
		       sMM[j-1] + q.tr[i][M2I] + t.tr[j-1][M2M_GAPOPEN], // MM->IM gap opening in query 
		       sIM[j-1] + q.tr[i][I2I] + t.tr[j-1][M2M], // IM->IM gap extension in query 
		       bIM[i][j]
		       );
	      sDG_i_j = max2
	              (
// 		       sMM[j] + q.tr[i-1][M2D],
// 		       sDG[j] + q.tr[i-1][D2D], //gap extension (DD) in query
		       sMM[j] + q.tr[i-1][M2D] + t.tr[j][GAPOPEN], // MM->DG gap opening in template 
		       sDG[j] + q.tr[i-1][D2D] + t.tr[j][GAPEXTD], // DG->DG gap extension in template 
		       bDG[i][j]
		       );
	      sMI_i_j = max2
	              (
		       sMM[j] + q.tr[i-1][M2M] + t.tr[j][M2I], // MM->MI gap opening M2I in template 
		       sMI[j] + q.tr[i-1][M2M] + t.tr[j][I2I], // MI->MI gap extension I2I in template 
		       bMI[i][j]
		       );

	      sMM_i_1_j_1 = sMM[j];
	      sGD_i_1_j_1 = sGD[j];
	      sIM_i_1_j_1 = sIM[j];
	      sDG_i_1_j_1 = sDG[j];
	      sMI_i_1_j_1 = sMI[j];
	      sMM[j] = sMM_i_j;
	      sGD[j] = sGD_i_j;
	      sIM[j] = sIM_i_j;
	      sDG[j] = sDG_i_j;
	      sMI[j] = sMI_i_j;

          //if (isnan(sMM_i_j)||isinf(sMM_i_j)){
          //  printf("."); /* <DEBUG> FS*/
          //}
	      // Find maximum score; global alignment: maxize only over last row and last column
	      if(sMM_i_j>score && (par.loc || i==q.L)) { i2=i; j2=j; score=sMM_i_j; }

	    } // end if 
      //printf("i= %d\tj= %d\ti2= %d\tj2= %d\tsMM= %f\tscore= %f\n", i, j, i2, j2, sMM_i_j, score);
	} //end for j
      
      // if global alignment: look for best cell in last column
      if (!par.loc && sMM_i_j>score) { i2=i; j2=jmax; score=sMM_i_j; }
      
    } // end for i

  state=MM; // state with maximum score is MM state

  // If local alignment do length correction: -log(length)
  if (par.loc)
    {
      if (self)
	score=score-log(0.5*t.L*q.L/200.0/200.0)/LAMDA - 11.2; // offset of -11.2 to get approx same mean as for -global
      else 
	if (par.idummy==0 && q.lamda>0) //////////////////////////////////////////////
	  score=score-log(t.L*q.L/200.0/200.0)/q.lamda - 11.2; // offset of -11.2 to get approx same mean as for -global
      else if (par.idummy<=1) //////////////////////////////////////////////
	  score=score-log(t.L*q.L/200.0/200.0)/LAMDA - 11.2; // offset of -11.2 to get approx same mean as for -global
    }  
//   printf("Template=%-12.12s  i=%-4i j=%-4i score=%6.3f\n",t.name,i2,j2,score);

  delete[] sMM; sMM = NULL;
  delete[] sGD; sGD = NULL;
  delete[] sDG; sDG = NULL;
  delete[] sIM; sIM = NULL;
  delete[] sMI; sMI = NULL;

  return;

} /* this is the end of Hit::Viterbi() */



/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compare two HMMs with Forward Algorithm in lin-space (~ 2x faster than in log-space)
 */
int  
Hit::Forward(HMM& q, HMM& t, float** Pstruc)
{

    // Variable declarations
    int i,j;      // query and template match state indices
    double pmin=(par.loc? 1.0: 0.0);    // used to distinguish between SW and NW algorithms in maximization         
    double Cshift = pow(2.0,par.shift);   // score offset transformed into factor in lin-space
    double Pmax_i;                        // maximum of F_MM in row i
    double scale_prod=1.0;                // Prod_i=1^i (scale[i])
    int jmin;  
    
    // First alignment of this pair of HMMs?
    if(irep==1) 
        {
            q.tr[0][M2D] = q.tr[0][M2I] = 0.0;
            q.tr[0][I2M] = q.tr[0][I2I] = 0.0;
            q.tr[0][D2M] = q.tr[0][D2D] = 0.0;
            t.tr[0][M2M] = 1.0;
            t.tr[0][M2D] = t.tr[0][M2I] = 0.0;
            t.tr[0][I2M] = t.tr[0][I2I] = 0.0;
            t.tr[0][D2M] = t.tr[0][D2D] = 0.0;
            q.tr[q.L][M2M] = 1.0;
            q.tr[q.L][M2D] = q.tr[q.L][M2I] = 0.0;
            q.tr[q.L][I2M] = q.tr[q.L][I2I] = 0.0;
            q.tr[q.L][D2M] = 1.0;
            q.tr[q.L][D2D] = 0.0;
            t.tr[t.L][M2M] = 1.0;
            t.tr[t.L][M2D] = t.tr[t.L][M2I] = 0.0;
            t.tr[t.L][I2M] = t.tr[t.L][I2I] = 0.0;
            t.tr[t.L][D2M] = 1.0;
            t.tr[t.L][D2D] = 0.0;
            InitializeForAlignment(q,t);
        }	
    
    
    // Initialization of top row, i.e. cells (0,j)
    F_MM[1][0] = F_IM[1][0] = F_GD[1][0] =  F_MM[0][1] = F_MI[0][1] = F_DG[0][1] = 0.0;
    for (j=1; j<=t.L; j++) 
        {
            if (cell_off[1][j]) 
                F_MM[1][j] = F_MI[1][j] = F_DG[1][j] = F_IM[1][j] = F_GD[1][j] = 0.0;
            else 
                {
                    F_MM[1][j] = ProbFwd(q.p[1],t.p[j]) * fpow2(ScoreSS(q,t,1,j)) * Cshift * (Pstruc==NULL? 1: Pstruc[1][j]) ;
                    F_MI[1][j] = F_DG[1][j] = 0.0;
                    F_IM[1][j] = F_MM[1][j-1] * q.tr[1][M2I] * t.tr[j-1][M2M] + F_IM[1][j-1] * q.tr[1][I2I] * t.tr[j-1][M2M];
                    F_GD[1][j] = F_MM[1][j-1] * t.tr[j-1][M2D]                + F_GD[1][j-1] * t.tr[j-1][D2D];
                }
        }
    scale[0]=scale[1]=scale[2]=1.0;
    
    // Forward algorithm
    for (i=2; i<=q.L; i++) // Loop through query positions i
        {
            //       if (v>=5) printf("\n");
            
            if (self) jmin = imin(i+SELFEXCL+1,t.L); else jmin=1;
            
            if (scale_prod<DBL_MIN*100) scale_prod = 0.0; else scale_prod *= scale[i];
            
            // Initialize cells at (i,0)
            if (cell_off[i][jmin]) 
                F_MM[i][jmin] = F_MI[i][jmin] = F_DG[i][jmin] = F_IM[i][jmin] = F_GD[i][jmin] = 0.0;
            else 
                {
                    F_MM[i][jmin] = scale_prod * ProbFwd(q.p[i],t.p[jmin]) * fpow2(ScoreSS(q,t,i,jmin)) * Cshift * (Pstruc==NULL? 1: Pstruc[i][jmin]);
                    F_IM[i][jmin] = F_GD[i][jmin] = 0.0; 
                    F_MI[i][jmin] = scale[i] * (F_MM[i-1][jmin] * q.tr[i-1][M2M] * t.tr[jmin][M2I] + F_MI[i-1][jmin] * q.tr[i-1][M2M] * t.tr[jmin][I2I]);
                    F_DG[i][jmin] = scale[i] * (F_MM[i-1][jmin] * q.tr[i-1][M2D]                   + F_DG[i-1][jmin] * q.tr[i-1][D2D]);
                }
            Pmax_i=0;
            
            for (j=jmin+1; j<=t.L; j++) // Loop through template positions j
                {
                    // Recursion relations
                    //	      printf("S[%i][%i]=%4.1f  ",i,j,Score(q.p[i],t.p[j]));
                    
                    if (cell_off[i][j]) 
                        F_MM[i][j] = F_MI[i][j] = F_DG[i][j] = F_IM[i][j] = F_GD[i][j] = 0.0;
                    else
                        {
                            F_MM[i][j] = ProbFwd(q.p[i],t.p[j]) * fpow2(ScoreSS(q,t,i,j)) * Cshift * (Pstruc==NULL? 1: Pstruc[i][j]) * scale[i] *
                                ( pmin
                                  + F_MM[i-1][j-1] * q.tr[i-1][M2M] * t.tr[j-1][M2M] // BB -> MM (BB = Begin/Begin, for local alignment)
                                  + F_GD[i-1][j-1] * q.tr[i-1][M2M] * t.tr[j-1][D2M] // GD -> MM
                                  + F_IM[i-1][j-1] * q.tr[i-1][I2M] * t.tr[j-1][M2M] // IM -> MM
                                  + F_DG[i-1][j-1] * q.tr[i-1][D2M] * t.tr[j-1][M2M] // DG -> MM
                                  + F_MI[i-1][j-1] * q.tr[i-1][M2M] * t.tr[j-1][I2M] // MI -> MM
                                  );
                            F_GD[i][j] = 
                                ( F_MM[i][j-1] * t.tr[j-1][M2D]                    // GD -> MM
                                  + F_GD[i][j-1] * t.tr[j-1][D2D]                    // GD -> GD
                                  + (Pstruc==NULL? 0 : F_DG[i][j-1] * t.tr[j-1][M2D] * q.tr[i][D2M] ) // DG -> GD (only when structure scores given)
                                  );
                            F_IM[i][j] = 
                                ( F_MM[i][j-1] * q.tr[i][M2I] * t.tr[j-1][M2M]     // MM -> IM
                                  + F_IM[i][j-1] * q.tr[i][I2I] * t.tr[j-1][M2M]     // IM -> IM
                                  + (Pstruc==NULL? 0 : F_MI[i][j-1] * q.tr[i][M2I] * t.tr[j-1][I2M] ) // MI -> IM (only when structure scores given)
                                  );
                            F_DG[i][j] = scale[i] * 
                                ( F_MM[i-1][j] * q.tr[i-1][M2D]                    // DG -> MM
                                  + F_DG[i-1][j] * q.tr[i-1][D2D]                    // DG -> DG
                                  ) ;
                            F_MI[i][j] = scale[i] * 
                                ( F_MM[i-1][j] * q.tr[i-1][M2M] * t.tr[j][M2I]     // MI -> MM 
                                  + F_MI[i-1][j] * q.tr[i-1][M2M] * t.tr[j][I2I]     // MI -> MI
                                  );
                            
                            if(F_MM[i][j]>Pmax_i) Pmax_i=F_MM[i][j];
                            
                        } // end else  	  
                    
                } //end for j
            
            pmin *= scale[i];
            scale[i+1] = 1.0/(Pmax_i+1.0);
            //      scale[i+1] = 1.0;
            
        } // end for i
    
    // Calculate P_forward * Product_{i=1}^{Lq+1}(scale[i])
    if (par.loc) 
        {
            Pforward = 1.0; // alignment contains no residues (see Mueckstein, Stadler et al.)
            for (i=1; i<=q.L; i++) // Loop through query positions i
                {
                    for (j=1; j<=t.L; j++) // Loop through template positions j
                        Pforward += F_MM[i][j];
                    Pforward *= scale[i+1];
                }
        }
    else  // global alignment
        {
            Pforward = 0.0;
            for (i=1; i<q.L; i++) {
                Pforward = (Pforward + F_MM[i][t.L]) * scale[i+1];
            }
            for (j=1; j<=t.L; j++) {
                Pforward += F_MM[q.L][j];
            }
            Pforward *= scale[q.L+1];
        }
    
    // Calculate log2(P_forward)
    score = log2(Pforward)-10.0f;
    for (i=1; i<=q.L+1; i++) score -= log2(scale[i]);
    //   state = MM;
    
    if (par.loc) 
        {
            if (self)
                score=score-log(0.5*t.L*q.L)/LAMDA+14.; // +14.0 to get approx same mean as for -global
            else 
                score=score-log(t.L*q.L)/LAMDA+14.; // +14.0 to get approx same mean as for -global
        }
    
    // Debugging output
    if (v>=6) 
        {
            const int i0=0, i1=q.L;
            const int j0=0, j1=t.L;
            scale_prod=1;
            printf("\nFwd      scale     ");
            for (j=j0; j<=j1; j++) printf("%3i     ",j);
            printf("\n");
            for (i=i0; i<=i1; i++) 
                {
                    scale_prod *= scale[i];
                    printf("%3i: %9.3G ",i,1/scale_prod);
                    for (j=j0; j<=j1; j++) 
                        printf("%7.4f ",(F_MM[i][j]+F_MI[i][j]+F_IM[i][j]+F_DG[i][j]+F_GD[i][j]));
                    printf("\n");
                    // 	  printf(" MM  %9.5f ",1/scale[i]);
                    // 	  for (j=j0; j<=j1; j++) 
                    // 	    printf("%7.4f ",F_MM[i][j]);
                    // 	  printf("\n");
                }
        }
    //   printf("Template=%-12.12s  score=%6.3f i2=%i  j2=%i \n",t.name,score,i2,j2);

    /* check for NaN and or infinities, FS, r241 -> r243 */
    if (isnan(score) || isinf(score) || isnan(Pforward) || isinf(Pforward) ){
        fprintf(stderr, "%s:%s:%d: Forward score is %g, Pforward is %g\n",
                __FUNCTION__, __FILE__, __LINE__, score, Pforward);
        return FAILURE;
    }
    /* alignment might fail if no useful characters in sequence, FS, r259 -> r260 */
    if ( (q.L <= 0) || (t.L <= 0) ){
        fprintf(stderr, "%s:%s:%d: length(s) of profile(s) invalid (q.L=%d/t.L=%d)\n",
                __FUNCTION__, __FILE__, __LINE__, q.L, t.L);
        return FAILURE;
    }
    i = q.L-1; j = t.L-1; /* FS, r241 -> r243 */
    if (isinf(F_MM[i][j]+F_MI[i][j]+F_IM[i][j]+F_DG[i][j]+F_GD[i][j])){
        fprintf(stderr, "%s:%s:%d: F_MM[i][j]=%g, F_IM[i][j]=%g, F_MI[i][j]=%g, F_DG[i][j]=%g, F_GD[i][j]=%g (i=%d,j=%d)\n", 
                __FUNCTION__, __FILE__, __LINE__, F_MM[i][j], F_MI[i][j], F_IM[i][j], F_DG[i][j], F_GD[i][j], i, j);
        return FAILURE;
    }
    return OK;

} /* this is the end of Hit::Forward() */





/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compare two HMMs with Backward Algorithm (in lin-space, 2x faster), for use in MAC alignment 
 */
int 
Hit::Backward(HMM& q, HMM& t)
{
    
    // Variable declarations
    int i,j;      // query and template match state indices
    double pmin=(par.loc? 1.0: 0.0);    // used to distinguish between SW and NW algorithms in maximization         
    double Cshift = pow(2.0,par.shift);   // score offset transformed into factor in lin-space
    double scale_prod=scale[q.L+1];
    int jmin;
    //double dMaxB = -1.0;
    
    // Initialization of top row, i.e. cells (0,j)
    for (j=t.L; j>=1; j--) 
        {
            if (cell_off[q.L][j]) 
                B_MM[q.L][j] = 0.0;
            else 
                B_MM[q.L][j] = scale[q.L+1];
            //dMaxB = dMaxB>B_MM[q.L][j]?dMaxB:B_MM[q.L][j];
            B_IM[q.L][j] = B_MI[q.L][j] = B_DG[q.L][j] = B_GD[q.L][j] = 0.0;
        }
    if (par.loc) pmin = scale[q.L+1]; // transform pmin (for local alignment) to scale of present (i'th) row 
    
    // Backward algorithm
    for (i=q.L-1; i>=1; i--) // Loop through query positions i
        {
            //       if (v>=5) printf("\n");
            
            if (self) jmin = imin(i+SELFEXCL,t.L); else jmin=1; // jmin = i+SELFEXCL and not (i+SELFEXCL+1) to set matrix element at boundary to zero
            
            // Initialize cells at (i,t.L+1)
            scale_prod *= scale[i+1];
            if (cell_off[i][t.L]) 
                B_MM[i][t.L] = 0.0;  
            else 
                B_MM[i][t.L] = scale_prod; 
            //if (isnan(B_MM[i][t.L])||isinf(B_MM[i][t.L])){
            //  printf("."); /* <DEBUG> FS*/
            //}
            //dMaxB = dMaxB>B_MM[i][t.L]?dMaxB:B_MM[i][t.L];
            B_IM[i][t.L] = B_MI[i][t.L] = B_DG[i][t.L] = B_GD[i][t.L] = 0.0; 
            pmin *= scale[i+1]; // transform pmin (for local alignment) to scale of present (i'th) row 
            
            for (j=t.L-1; j>=jmin; j--) // Loop through template positions j
                {
                    // Recursion relations
                    //	      printf("S[%i][%i]=%4.1f  ",i,j,Score(q.p[i],t.p[j]));
                    if (cell_off[i][j]) 
                        B_MM[i][j] = B_GD[i][j] = B_IM[i][j] = B_DG[i][j] = B_MI[i][j] = 0.0;  
                    else 
                        {
                            double pmatch = B_MM[i+1][j+1] * ProbFwd(q.p[i+1],t.p[j+1]) * fpow2(ScoreSS(q,t,i+1,j+1)) * Cshift * scale[i+1];
                            //if (isnan(pmatch)||isinf(pmatch)){
                            //  printf("."); /* <DEBUG> FS*/
                            //}
                            B_MM[i][j] =  
                                (
                                 + pmin                                                    // MM -> EE (End/End, for local alignment)
                                 + pmatch       * q.tr[i][M2M] * t.tr[j][M2M]              // MM -> MM
                                 + B_GD[i][j+1]                * t.tr[j][M2D]              // MM -> GD (q.tr[i][M2M] is already contained in GD->MM)
                                 + B_IM[i][j+1] * q.tr[i][M2I] * t.tr[j][M2M]              // MM -> IM
                                 + B_DG[i+1][j] * q.tr[i][M2D]                * scale[i+1] // MM -> DG (t.tr[j][M2M] is already contained in DG->MM)
                                 + B_MI[i+1][j] * q.tr[i][M2M] * t.tr[j][M2I] * scale[i+1] // MM -> MI
                                 );
                            //if (isnan(B_MM[i][j])||isinf(B_MM[i][j])){
                            //  printf("."); /* <DEBUG> FS*/
                            //}
                            //dMaxB = dMaxB>B_MM[i][j]?dMaxB:B_MM[i][j];

                            B_GD[i][j] = 
                                (
                                 + pmatch       * q.tr[i][M2M] * t.tr[j][D2M]              // GD -> MM 
                                 + B_GD[i][j+1]                * t.tr[j][D2D]              // DG -> DG   
                                 );
                            B_IM[i][j] = 
                                (
                                 + pmatch       * q.tr[i][I2M] * t.tr[j][M2M]              // IM -> MM
                                 + B_IM[i][j+1] * q.tr[i][I2I] * t.tr[j][M2M]              // IM -> IM
                                 );
                            B_DG[i][j] =  
                                (
                                 + pmatch       * q.tr[i][D2M] * t.tr[j][M2M]              // DG -> MM
                                 + B_DG[i+1][j] * q.tr[i][D2D]                * scale[i+1] // DG -> DG
                                 //   	         + B_GD[i][j+1] * q.tr[i][D2M] * t.tr[j][M2D]              // DG -> GD
                                 );
                            B_MI[i][j] = 
                                (
                                 + pmatch       * q.tr[i][M2M] * t.tr[j][I2M]              // MI -> MM       
                                 + B_MI[i+1][j] * q.tr[i][M2M] * t.tr[j][I2I] * scale[i+1] // MI -> MI
                                 // 	         + B_IM[i][j+1] * q.tr[i][M2I] * t.tr[j][I2M]              // MI -> IM    
                                 );
                            
                        } // end else	      
                    
                } //end for j
            
        } // end for i
    
    // Debugging output
    if (v>=6)
        {
            const int i0=0, i1=q.L;
            const int j0=0, j1=t.L;
            double scale_prod[q.L+2];
            scale_prod[q.L] = scale[q.L+1];
            for (i=q.L-1; i>=1; i--) scale_prod[i] = scale_prod[i+1] * scale[i+1];
            
            printf("\nBwd      scale     ");
            for (j=j0; j<=j1; j++) printf("%3i     ",j);
            printf("\n");
            for (i=i0; i<=i1; i++) 
                {
                    printf("%3i: %9.3G ",i,1/scale_prod[i]);
                    for (j=j0; j<=j1; j++)
                        printf("%7.4f ",(B_MM[i][j]+B_MI[i][j]+B_IM[i][j]+B_DG[i][j]+B_GD[i][j]) * (ProbFwd(q.p[i],t.p[j])*fpow2(ScoreSS(q,t,i,j)) * Cshift));
                    printf("\n");
                    
                    // 	  printf("MM   %9.5f ",1/scale[i]);
                    // 	  for (j=j0; j<=j1; j++)
                    // 	    printf("%7.4f ",B_MM[i][j] * (ProbFwd(q.p[i],t.p[j])*fpow2(ScoreSS(q,t,i,j)) * Cshift));
                    // 	  printf("\n");
                }
            printf("\nPost     scale     ");
            for (j=j0; j<=j1; j++) printf("%3i     ",j);
            printf("\n");
            for (i=i0; i<=i1; i++) 
                {
                    printf("%3i: %9.3G ",i,1/scale_prod[i]);
                    for (j=j0; j<=j1; j++) 
                        printf("%7.4f ",B_MM[i][j]*F_MM[i][j]/Pforward);
                    printf("\n");
                }
            printf("\n");
        }
    
    if (v>=4) printf("\nForward total probability ratio: %8.3G\n",Pforward);
    
    // Calculate Posterior matrix and overwrite Backward matrix with it
    for (i=1; i<=q.L; i++) {
        for (j=1; j<=t.L; j++) { 
            B_MM[i][j] *= F_MM[i][j]/Pforward;
            //if (isnan(B_MM[i][j]) || isinf(B_MM[i][j])){
            //  printf("."); /* <DEBUG> FS*/
            //}
            //dMaxB = dMaxB>B_MM[i][j]?dMaxB:B_MM[i][j];
        }
    }

    //printf("Max-B_MM = %f\n", dMaxB);

    /* check for NaN and or infinities, FS, r241 -> r243 */
    if (isnan(score) || isinf(score)){
        fprintf(stderr, "%s:%s:%d: Backward score is %g\n",
                __FUNCTION__, __FILE__, __LINE__, score);
        return FAILURE;
    }
    i = j = 1;
    if (isinf(B_MM[i][j]+B_MI[i][j]+B_IM[i][j]+B_DG[i][j]+B_GD[i][j])){
        fprintf(stderr, "%s:%s:%d: B_MM[1][1]=%g, B_IM[1][1]=%g, B_MI[1][1]=%g, B_DG[1][1]=%g, B_GD[1][1]=%g\n", 
                __FUNCTION__, __FILE__, __LINE__, B_MM[i][j], B_MI[i][j], B_IM[i][j], B_DG[i][j], B_GD[i][j]);
        for (i = 1; (i < q.L) && isinf(B_MM[i][1]); i++);
        i--;
        for (j = 1; (j < t.L) && isinf(B_MM[i][j]); j++);
        j--;
        fprintf(stderr, "%s:%s:%d: B_MM[%d][%d]=%g, B_MM[%d][%d]=%g, B_MM[%d][%d]=%g\n",
                __FUNCTION__, __FILE__, __LINE__, 
                i, j, B_MM[i][j], i+1, 1, B_MM[i+1][1], i, j+1, B_MM[i][j+1]);
        return FAILURE;
    }
    return OK;
    
} /* this is the end of Hit::Backward() */



/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Maximum Accuracy alignment 
 */
void 
Hit::MACAlignment(HMM& q, HMM& t)
{
  // Use Forward and Backward matrices to find that alignment which 
  // maximizes the expected number of correctly aligned pairs of residues (mact=0)
  // or, more generally, which maximizes the expectation value of the number of 
  // correctly aligned pairs minus (mact x number of aligned pairs)
  // "Correctly aligned" can be based on posterior probabilities calculated with
  // a local or a global version of the Forward-Backward algorithm.

  int i,j;           // query and template match state indices
  int jmin,jmax;     // range of dynamic programming for j
  double** S=F_MI;    // define alias for new score matrix
  double score_MAC;   // score of the best MAC alignment

  // Initialization of top row, i.e. cells (0,j)
  for (j=0; j<=t.L; j++) S[0][j] = 0.0;
  score_MAC=-INT_MAX; i2=j2=0; bMM[0][0]=STOP;

  // Dynamic programming 
  for (i=1; i<=q.L; i++) // Loop through query positions i
    {
      
      if (self) 
	{
	  // If q is compared to itself, ignore cells below diagonal+SELFEXCL
	  jmin = i+SELFEXCL; 
	  jmax = t.L;
	  if (jmin>jmax) continue;
	}
      else
	{
	  // If q is compared to t, exclude regions where overlap of q with t < min_overlap residues
	  jmin=imax( 1, i+min_overlap-q.L);  // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq} 
	  jmax=imin(t.L,i-min_overlap+t.L);  // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt} 
	}      

      // Initialize cells
      S[i][jmin-1] = 0.0;
      if (jmax<t.L) S[i-1][jmax] = 0.0; // initialize at (i-1,jmax) if upper right triagonal is excluded due to min overlap
      
      for (j=jmin; j<=jmax; j++) // Loop through template positions j
	{

	  if (cell_off[i][j]) 
	    S[i][j] = -FLT_MIN;
	  else 
	    {
	      // Recursion
	     
	      // NOT the state before the first MM state)
	      CALCULATE_MAX4(
		 S[i][j],
		 B_MM[i][j] - par.mact,  // STOP signifies the first MM state, NOT the state before the first MM state (as in Viterbi)
		 S[i-1][j-1] + B_MM[i][j] - par.mact, // B_MM[i][j] contains posterior probability
		 S[i-1][j] - 0.5*par.mact,  // gap penalty prevents alignments such as this: XX--xxXX
		 S[i][j-1] - 0.5*par.mact,  //                                               YYyy--YY  
		 bMM[i][j]   // backtracing matrix
		 );

// 	      if (i==6 && j==8) 
// 		printf("i=%i  j=%i  S[i][j]=%8.3f  MM:%7.3f  MI:%7.3f  IM:%7.3f  b:%i\n",i,j,S[i][j],S[i-1][j-1]+B_MM[i][j]-par.mact,S[i-1][j],S[i][j-1],bMM[i][j]);
	      
	      // Find maximum score; global alignment: maximize only over last row and last column
	      if(S[i][j]>score_MAC && (par.loc || i==q.L)) { i2=i; j2=j; score_MAC=S[i][j]; }	      
	      
	    } // end if 
	  
	} //end for j
      
	  // if global alignment: look for best cell in last column
      if (!par.loc && S[i][jmax]>score_MAC) { i2=i; j2=jmax; score_MAC=S[i][jmax]; }
      
    } // end for i
  
  // DEBUG
  if (v>=5) 
    {
      printf("\nScore  ");
      for (j=0; j<=t.L; j++) printf("%3i   ",j);
      printf("\n");
      for (i=0; i<=q.L; i++) 
	{
	  printf("%2i:    ",i);
 	  for (j=0; j<=t.L; j++) 
	    printf("%5.2f ",S[i][j]);
	  printf("\n");
	}
      printf("\n");
      printf("Template=%-12.12s  i=%-4i j=%-4i score=%6.3f\n",t.name,i2,j2,score);
    }  

  return;

} /* this is the end of Hit::MACAlignment() */


/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Trace back alignment of two profiles based on matrices bXX[][]
 */
void 
Hit::Backtrace(HMM& q, HMM& t)
{
  // Trace back trough the matrices bXY[i][j] until first match state is found (STOP-state)

  int step;      // counts steps in path through 5-layered dynamic programming matrix
  int i,j;       // query and template match state indices

  InitializeBacktrace(q,t);
  
  // Make sure that backtracing stops when t:M1 or q:M1 is reached (Start state), e.g. sMM[i][1], or sIM[i][1] (M:MM, B:IM)
  for (i=0; i<=q.L; i++) bMM[i][1]=bGD[i][1]=bIM[i][1] = STOP;
  for (j=1; j<=t.L; j++) bMM[1][j]=bDG[1][j]=bMI[1][j] = STOP;
  

  // Back-tracing loop
  matched_cols=0; // for each MACTH (or STOP) state matched_col is incremented by 1
  step=0;         // steps through the matrix correspond to alignment columns (from 1 to nsteps)
  //  state=MM;       // state with maximum score must be MM state  // already set at the end of Viterbi()
  i=i2; j=j2;     // last aligned pair is (i2,j2)
  while (state)   // while (state!=STOP)  because STOP=0
    {
      step++;
      states[step] = state;
      this->i[step] = i;
      this->j[step] = j;
      // Exclude cells in direct neighbourhood from all further alignments
      for (int ii=imax(i-2,1); ii<=imin(i+2,q.L); ii++)
	cell_off[ii][j]=1;     
      for (int jj=imax(j-2,1); jj<=imin(j+2,t.L); jj++)
	cell_off[i][jj]=1;     
      
      switch (state)
	{
	case MM: // current state is MM, previous state is bMM[i][j]
	  matched_cols++; 
	  state = bMM[i--][j--];
	  break;	      
	case GD: // current state is GD
	  switch (bGD[i][j--])
	    {
	    case STOP: state = STOP; break; // current state does not have predecessor
	    case MM:   state = MM;   break; // previous state is Match state
	    }                               // default: previous state is same state (GD)
	  break;	      
	case IM: 
	  switch (bIM[i][j--]) 
	    {
	    case STOP: state = STOP; break; // current state does not have predecessor
	    case MM:   state = MM;   break; // previous state is Match state
	    }                               // default: previous state is same state (IM)
	  break;	      
	case DG:
	  switch (bDG[i--][j])
	    {
	    case STOP: state = STOP; break; // current state does not have predecessor
	    case MM:   state = MM;   break; // previous state is Match state
	    }                               // default: previous state is same state (DG)
	  break;	      
	case MI:
	  switch (bMI[i--][j])
	    {
	    case STOP: state = STOP; break; // current state does not have predecessor
	    case MM:   state = MM;   break; // previous state is Match state
		}                               // default: previous state is same state (MI)
	  break;
	default:
	  fprintf(stderr,"Error: unallowed state value %i occurred during backtracing at step %i, (i,j)=(%i,%i)\n",state,step,i,j);
	  state=0;
	  v=4;
	  break;
	} //end switch (state)
    } //end while (state)
 
  i1 = this->i[step];
  j1 = this->j[step];
  states[step] = MM;  // first state (STOP state) is set to MM state
  nsteps=step; 
  
  // Allocate new space for alignment scores
  if (t.Xcons) Xcons = new( char[q.L+2]); // for template consensus sequence aligned to query
  S    = new( float[nsteps+1]);
  S_ss = new( float[nsteps+1]);
  if (!S_ss) MemoryError("space for HMM-HMM alignments");

  // Add contribution from secondary structure score, record score along alignment,
  // and record template consensus sequence in master-slave-alignment to query sequence
  score_ss=0.0f;
  int ssm=ssm1+ssm2;
  for (step=1; step<=nsteps; step++)
    {
      switch(states[step])
	{
	case MM: 
	  i = this->i[step];
	  j = this->j[step];
	  S[step] = Score(q.p[i],t.p[j]);
	  S_ss[step] = ScoreSS(q,t,i,j,ssm);
	  score_ss += S_ss[step];
	  if (Xcons) Xcons[i]=t.Xcons[j]; //record database consensus sequence
	  break;
	case MI: //if gap in template  
	case DG:   
	  if (Xcons) Xcons[this->i[step]]=GAP; //(no break hereafter)
	default: //if gap in T or Q
	  S[step]=S_ss[step]=0.0f;
	  break;
	}
    }
  if (ssm2>=1) score-=score_ss;    // subtract SS score added during alignment!!!!
  if (Xcons) 
    {
      for (i=0; i<i1; i++) Xcons[i]=ENDGAP; // set end gap code at beginning and end of template consensus sequence
      for (i=i2+1; i<=q.L+1; i++) Xcons[i]=ENDGAP;
    }
  
  // Add contribution from correlation of neighboring columns to score
  float Scorr=0;
  if (nsteps)
    {
      for (step=2; step<=nsteps; step++) Scorr+=S[step]*S[step-1];
      for (step=3; step<=nsteps; step++) Scorr+=S[step]*S[step-2];
      for (step=4; step<=nsteps; step++) Scorr+=S[step]*S[step-3];
      for (step=5; step<=nsteps; step++) Scorr+=S[step]*S[step-4];
      score+=par.corr*Scorr;
    }
  
  // Set score, P-value etc.
  score_sort = score_aass = -score;
  logPval=0; Pval=1;
  if (t.mu)
    {
      logPvalt=logPvalue(score,t.lamda,t.mu); 
      Pvalt=Pvalue(score,t.lamda,t.mu); 
    }
  else { logPvalt=0; Pvalt=1;}
  //   printf("%-10.10s lamda=%-9f  score=%-9f  logPval=%-9g\n",name,t.lamda,score,logPvalt);
  

  //DEBUG: Print out Viterbi path
  if (v>=4) 
    {
      printf("NAME=%7.7s score=%7.3f  score_ss=%7.3f\n",name,score,score_ss);
      printf("step  Q T    i    j  state   score    T Q cf ss-score\n");
      for (step=nsteps; step>=1; step--)
	{
	  switch(states[step])
	    {
	    case MM: 
	      printf("%4i  %1c %1c ",step,q.seq[q.nfirst][this->i[step]],seq[nfirst][this->j[step]]); 
	      break;
	    case GD: 
	    case IM: 
	      printf("%4i  - %1c ",step,seq[nfirst][this->j[step]]); 
	      break;
	    case DG:
	    case MI: 
	      printf("%4i  %1c - ",step,q.seq[q.nfirst][this->i[step]]); 
	      break;
	    }
	  printf("%4i %4i     %2i %7.2f    ",this->i[step],this->j[step],(int)states[step],S[step]); 
	  printf("%c %c %1i %7.2f\n",i2ss(t.ss_dssp[this->j[step]]),i2ss(q.ss_pred[this->i[step]]),q.ss_conf[this->i[step]]-1,S_ss[step]); 
	}
    }

 return;

} /* this is the end of Hit::Backtrace() */



/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief GLOBAL stochastic trace back through the forward matrix of probability ratios
 */
void 
Hit::StochasticBacktrace(HMM& q, HMM& t, char maximize)
{
  int step;        // counts steps in path through 5-layered dynamic programming matrix
  int i,j;         // query and template match state indices
//  float pmin=(par.loc? 1.0: 0.0);    // used to distinguish between SW and NW algorithms in maximization         
  const float pmin=0;
  double* scale_cum = new(double[q.L+2]);
  

  scale_cum[0]=1;
  for (i=1; i<=q.L+1; i++) scale_cum[i] = scale_cum[i-1]*scale[i];

  // Select start cell for GLOBAL alignment
  // (Implementing this in a local version would make this method work for local backtracing as well)
  if (maximize) 
    {
      double F_max=0;
      for (i=q.L-1; i>=1; i--) 
	if (F_MM[i][t.L]/scale_cum[i]>F_max) {i2=i; j2=t.L; F_max=F_MM[i][t.L]/scale_cum[i];}
      for (j=t.L; j>=1; j--) 
	if (F_MM[q.L][j]/scale_cum[q.L]>F_max) {i2=q.L; j2=j; F_max=F_MM[q.L][j]/scale_cum[q.L];}
    }
  else 
    {
//      float sumF[q.L+t.L];
      double* sumF=new(double[q.L+t.L]);
      sumF[0]=0.0;
      for (j=1; j<=t.L; j++)        sumF[j] = sumF[j-1] + F_MM[q.L][j]/scale_cum[q.L];;
      for (j=t.L+1; j<t.L+q.L; j++) sumF[j] = sumF[j-1] + F_MM[j-t.L][t.L]/scale_cum[j-t.L];;
      float x = sumF[t.L+q.L-1]*frand(); // generate random number between 0 and sumF[t.L+q.L-1]
      for (j=1; j<t.L+q.L; j++) 
	if (x<sumF[j]) break;
      if (j<=t.L) {i2=q.L; j2=j;} else {i2=j-t.L; j2=t.L;}
      delete[] sumF; sumF = NULL;
    }

  InitializeBacktrace(q,t);

  int (*pick2)(const double&, const double&, const int&);
  int (*pick3_GD)(const double&, const double&, const double&);
  int (*pick3_IM)(const double&, const double&, const double&);
  int (*pick6)(const double&, const double&, const double&, const double&, const double&, const double&);
  if (maximize) 
    {
      pick2 = &pickmax2;
      pick3_GD = &pickmax3_GD;      
      pick3_IM = &pickmax3_IM;
      pick6 = &pickmax6;
    }
  else 
    {
      pick2 = &pickprob2;
      pick3_GD = &pickprob3_GD;
      pick3_IM = &pickprob3_IM;
      pick6 = &pickprob6;
    }

  // Back-tracing loop
  matched_cols=0;     // for each MACTH (or STOP) state matched_col is incremented by 1
  step=0;             // steps through the matrix correspond to alignment columns (from 1 to nsteps)
  state = MM;
  i=i2; j=j2;    // start at end of query and template
  while (state)  // while not reached STOP state or upper or left border 
    {
      step++;
      states[step] = state;
      this->i[step] = i;
      this->j[step] = j;

      switch (state)
	{
	  
	case MM: // current state is MM, previous state is state
// 	  fprintf(stderr,"%4i  %1c %1c %4i %4i     MM %7.2f\n",step,q.seq[q.nfirst][i],seq[nfirst][j],i,j,Score(q.p[i],t.p[j])); 
// 	  printf("0:%7.3f   MM:%7.3f   GD:%7.3f   IM:%7.3f   DG:%7.3f   MI:%7.3f \n",
// 		        pmin*scale_cum[i-1],
// 		        F_MM[i-1][j-1] * q.tr[i-1][M2M] * t.tr[j-1][M2M], 
// 			F_GD[i-1][j-1] * q.tr[i-1][M2M] * t.tr[j-1][D2M],
// 			F_IM[i-1][j-1] * q.tr[i-1][I2M] * t.tr[j-1][M2M],
// 			F_DG[i-1][j-1] * q.tr[i-1][D2M] * t.tr[j-1][M2M],
// 		        F_MI[i-1][j-1] * q.tr[i-1][M2M] * t.tr[j-1][I2M]);
	  matched_cols++; 
	  if (j>1 && i>1)
	    state = (*pick6)( 
			pmin*scale_cum[i-1],
			F_MM[i-1][j-1] * q.tr[i-1][M2M] * t.tr[j-1][M2M], 
			F_GD[i-1][j-1] * q.tr[i-1][M2M] * t.tr[j-1][D2M],
			F_IM[i-1][j-1] * q.tr[i-1][I2M] * t.tr[j-1][M2M],
			F_DG[i-1][j-1] * q.tr[i-1][D2M] * t.tr[j-1][M2M],
			F_MI[i-1][j-1] * q.tr[i-1][M2M] * t.tr[j-1][I2M]
			);
	  else state=0;	  
	  i--; j--;
	  break;	      
	case GD: // current state is GD
// 	  fprintf(stderr,"%4i  - %1c %4i %4i     GD %7.2f\n",step,q.seq[q.nfirst][j],i,j,Score(q.p[i],t.p[j])); 
	  if (j>1) 
	    state = (*pick3_GD)(
			F_MM[i][j-1] * t.tr[j-1][M2D],
                        F_DG[i][j-1] * t.tr[j-1][M2D] * q.tr[i][D2M],   // DG -> GD
			F_GD[i][j-1] * t.tr[j-1][D2D]                   // gap extension (DD) in template
			);
	  else state=0;	  
	  j--;
	  break;	      
	case IM: 
// 	  fprintf(stderr,"%4i  - %1c %4i %4i     IM %7.2f\n",step,q.seq[q.nfirst][j],i,j,Score(q.p[i],t.p[j])); 
	  if (j>1) 
	    state = (*pick3_IM)(
			F_MM[i][j-1] * q.tr[i][M2I] * t.tr[j-1][M2M_GAPOPEN],
			F_MI[i][j-1] * q.tr[i][M2I] * t.tr[j-1][I2M],  // MI -> IM
			F_IM[i][j-1] * q.tr[i][I2I] * t.tr[j-1][M2M]   // gap extension (II) in query
			); 
	  else state=0;	  
	  j--;
	  break;	      
	case DG:
// 	  fprintf(stderr,"%4i  %1c - %4i %4i     DG %7.2f\n",step,q.seq[q.nfirst][i],i,j,Score(q.p[i],t.p[j])); 
	  if (i>1) 
	    state = (*pick2)(
			F_MM[i-1][j] * q.tr[i-1][M2D] * t.tr[j][GAPOPEN],
			F_DG[i-1][j] * q.tr[i-1][D2D] * t.tr[j][GAPEXTD], //gap extension (DD) in query
			DG
			);
	  else state=0;	  
	  i--; 
	  break;	      
	case MI:
// 	  fprintf(stderr,"%4i  %1c - %4i %4i     MI %7.2f\n",step,q.seq[q.nfirst][i],i,j,Score(q.p[i],t.p[j])); 
	  if (i>1) 
	    state = (*pick2)(
			F_MM[i-1][j] * q.tr[i-1][M2M] * t.tr[j][M2I],
			F_MI[i-1][j] * q.tr[i-1][M2M] * t.tr[j][I2I], //gap extension (II) in template
			MI
			);
	  else state=0;
	  i--; 
	  break;

	} //end switch (state)

    } //end while (state)
 
  i1 = this->i[step];
  j1 = this->j[step];
  states[step] = MM;  // first state (STOP state) is set to MM state
  nsteps=step; 

  // Allocate new space for alignment scores
  if (t.Xcons) Xcons = new( char[q.L+2]); // for template consensus sequence aligned to query
  S    = new( float[nsteps+1]);
  S_ss = new( float[nsteps+1]);
  if (!S_ss) MemoryError("space for HMM-HMM alignments");

  // Add contribution from secondary structure score, record score along alignment,
  // and record template consensus sequence in master-slave-alignment to query sequence
  score_ss=0.0f;
  int ssm=ssm1+ssm2;
  for (step=1; step<=nsteps; step++)
    {
      switch(states[step])
	{
	case MM: 
	  i = this->i[step];
	  j = this->j[step];
	  S[step] = Score(q.p[i],t.p[j]);
	  S_ss[step] = ScoreSS(q,t,i,j,ssm);
	  score_ss += S_ss[step];
	  if (Xcons) Xcons[i]=t.Xcons[j]; //record database consensus sequence
	  break;
	case MI: //if gap in template  
	case DG:   
	  if (Xcons) Xcons[this->i[step]]=GAP; //(no break hereafter)
	default: //if gap in T or Q
	  S[step]=S_ss[step]=0.0f;
	  break;
	}
    }
  if (ssm2>=1) score-=score_ss;    // subtract SS score added during alignment!!!!
  if (Xcons) 
    {
      for (i=0; i<i1; i++) Xcons[i]=ENDGAP; // set end gap code at beginning and end of template consensus sequence
      for (i=i2+1; i<=q.L+1; i++) Xcons[i]=ENDGAP;
    }

  delete[] scale_cum; scale_cum = NULL;

  return;
}





/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Trace back alignment of two profiles based on matrices bXX[][]
 */
void 
Hit::BacktraceMAC(HMM& q, HMM& t)
{
  // Trace back trough the matrix b[i][j] until STOP state is found

  char** b=bMM;  // define alias for backtracing matrix
  int step;      // counts steps in path through 5-layered dynamic programming matrix
  int i,j;       // query and template match state indices

  InitializeBacktrace(q,t);
  
  // Make sure that backtracing stops when t:M1 or q:M1 is reached (Start state), e.g. sMM[i][1], or sIM[i][1] (M:MM, B:IM)
  for (i=0; i<=q.L; i++) b[i][1] = STOP;
  for (j=1; j<=t.L; j++) b[1][j] = STOP;
  

  // Back-tracing loop
  // In contrast to the Viterbi-Backtracing, STOP signifies the first Match-Match state, NOT the state before the first MM state
  matched_cols=1; // for each MACTH (or STOP) state matched_col is incremented by 1
  state=MM;       // lowest state with maximum score must be match-match state 
  step=0;         // steps through the matrix correspond to alignment columns (from 1 to nsteps)
  i=i2; j=j2;     // last aligned pair is (i2,j2)
  while (state!=STOP) 
    {
      step++;
      states[step] = state = b[i][j];
      this->i[step] = i;
      this->j[step] = j;
      // Exclude cells in direct neighbourhood from all further alignments
      for (int ii=imax(i-2,1); ii<=imin(i+2,q.L); ii++)
	cell_off[ii][j]=1;     
      for (int jj=imax(j-2,1); jj<=imin(j+2,t.L); jj++)
	cell_off[i][jj]=1;     
      if (state==MM) matched_cols++; 

      switch (state)
	{
	case MM: i--; j--; break;
	case IM: j--; break;
	case MI: i--; break;
	case STOP: break;
	default:
	  fprintf(stderr,"Error: unallowed state value %i occurred during backtracing at step %i, (i,j)=(%i,%i)\n",state,step,i,j);
	  state=0;
	  v=4;
	  break;
	} //end switch (state)
    } //end while (state)
 
  i1 = this->i[step];
  j1 = this->j[step];
  states[step] = MM;  // first state (STOP state) is set to MM state
  nsteps=step; 
    
  // Allocate new space for alignment scores
  if (t.Xcons) Xcons = new( char[q.L+2]); // for template consensus sequence aligned to query
  S    = new( float[nsteps+1]);
  S_ss = new( float[nsteps+1]);
  P_posterior = new( float[nsteps+1]);
  if (!P_posterior) MemoryError("space for HMM-HMM alignments");

  // Add contribution from secondary structure score, record score along alignment,
  // and record template consensus sequence in master-slave-alignment to query sequence
  score_ss=0.0f;
  sum_of_probs=0.0;       // number of identical residues in query and template sequence
  int ssm=ssm1+ssm2;
//   printf("Hit=%s\n",name); /////////////////////////////////////////////////////////////
  for (step=1; step<=nsteps; step++)
    {
      switch(states[step])
	{
	case MM: 
	  i = this->i[step];
	  j = this->j[step];
	  S[step] = Score(q.p[i],t.p[j]);
	  S_ss[step] = ScoreSS(q,t,i,j,ssm);
	  score_ss += S_ss[step];
	  P_posterior[step] = B_MM[this->i[step]][this->j[step]];
	  // Add probability to sum of probs if no dssp states given or dssp states exist and state is resolved in 3D structure
	  if (t.nss_dssp<0 || t.ss_dssp[j]>0) sum_of_probs += P_posterior[step]; 
// 	  printf("j=%-3i  dssp=%1i  P=%4.2f  sum=%6.2f\n",j,t.ss_dssp[j],P_posterior[step],sum_of_probs); //////////////////////////
	  if (Xcons) Xcons[i]=t.Xcons[j]; //record database consensus sequence
	  break;
	case MI: //if gap in template  
	case DG:   
	  if (Xcons) Xcons[this->i[step]]=GAP; //(no break hereafter)
	default: //if gap in T or Q
	  S[step] = S_ss[step] = P_posterior[step] = 0.0;
	  break;
	}
    }
//   printf("\n"); /////////////////////////////////////////////////////////////
  if (ssm2>=1) score-=score_ss;    // subtract SS score added during alignment!!!!
  if (Xcons) 
    {
      for (i=0; i<i1; i++) Xcons[i]=ENDGAP; // set end gap code at beginning and end of template consensus sequence
      for (i=i2+1; i<=q.L+1; i++) Xcons[i]=ENDGAP;
    }

  // Add contribution from correlation of neighboring columns to score
  float Scorr=0;
  if (nsteps)
    {
      for (step=1; step<=nsteps-1; step++) Scorr+=S[step]*S[step+1];
      for (step=1; step<=nsteps-2; step++) Scorr+=S[step]*S[step+2];
      for (step=1; step<=nsteps-3; step++) Scorr+=S[step]*S[step+3];
      for (step=1; step<=nsteps-4; step++) Scorr+=S[step]*S[step+4];
      score+=par.corr*Scorr;
    }
  
  // Set score, P-value etc.
  score_sort = score_aass = -score;
  logPval=0; Pval=1;
  if (t.mu)
    {
      logPvalt=logPvalue(score,t.lamda,t.mu); 
      Pvalt=Pvalue(score,t.lamda,t.mu); 
    }
  else { logPvalt=0; Pvalt=1;}
//   printf("%-10.10s lamda=%-9f  score=%-9f  logPval=%-9g\n",name,t.lamda,score,logPvalt);


  //DEBUG: Print out MAC alignment path
  if (v>=4) 
    {
      float sum_post=0.0;
      printf("NAME=%7.7s score=%7.3f  score_ss=%7.3f\n",name,score,score_ss);
      printf("step  Q T    i    j  state   score    T Q cf ss-score   P_post Sum_post\n");
      for (step=nsteps; step>=1; step--)
	{
	  switch(states[step])
	    {
	    case MM: 
	      sum_post+=P_posterior[step];
	      printf("%4i  %1c %1c ",step,q.seq[q.nfirst][this->i[step]],seq[nfirst][this->j[step]]); 
	      break;
	    case IM: 
	      printf("%4i  - %1c ",step,seq[nfirst][this->j[step]]); 
	      break;
	    case MI: 
	      printf("%4i  %1c - ",step,q.seq[q.nfirst][this->i[step]]); 
	      break;
	    }
	  printf("%4i %4i     %2i %7.1f    ",this->i[step],this->j[step],(int)states[step],S[step]); 
	  printf("%c %c  %1i  %7.1f  ",i2ss(t.ss_dssp[this->j[step]]),i2ss(q.ss_pred[this->i[step]]),q.ss_conf[this->i[step]]-1,S_ss[step]); 
	  printf("%7.5f  %7.2f\n",P_posterior[step],sum_post); 
	}
    }

 return;
}



/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Functions that calculate probabilities
 */
void 
Hit::InitializeForAlignment(HMM& q, HMM& t)
{
  int i,j;

  // SS scoring during (ssm2>0) or after (ssm1>0) alignment? Query SS known or Template SS known?
  switch (par.ssm) 
    {
    case 0:
      ssm1=0;
      ssm2=0;
      break;
    case 1:
      ssm2=0;  // SS scoring after alignment
      if (t.nss_dssp>=0 && q.nss_pred>=0) ssm1=1;
      else if (q.nss_dssp>=0 && t.nss_pred>=0) ssm1=2;    
      else if (q.nss_pred>=0 && t.nss_pred>=0) ssm1=3;
      else ssm1=0;
      break;
    case 2:
      ssm1=0;  // SS scoring during alignment
      if (t.nss_dssp>=0 && q.nss_pred>=0) ssm2=1;
      else if (q.nss_dssp>=0 && t.nss_pred>=0) ssm2=2;   
      else if (q.nss_pred>=0 && t.nss_pred>=0) ssm2=3;
      else ssm2=0;
      break;
    case 3:
      ssm2=0;  // SS scoring after alignment
      if (q.nss_pred>=0 && t.nss_pred>=0) ssm1=3; else ssm1=0;  
      break;
    case 4:
      ssm1=0;  // SS scoring during alignment
      if (q.nss_pred>=0 && t.nss_pred>=0) ssm2=3; else ssm2=0;
      break;
      //     case 5:
      //       ssm2=0;  // SS scoring after alignment
      //       if (q.nss_dssp>=0 && t.nss_dssp>=0) ssm1=4; else ssm1=0;  
      //       break;
      //     case 6:
      //       ssm1=0;  // SS scoring during alignment
      //       if (q.nss_dssp>=0 && t.nss_dssp>=0) ssm2=4; else ssm2=0;
      //       break;
    }

  if (self)  
    {
      // Cross out cells in lower diagonal for self-comparison?
      for (i=1; i<=q.L; i++) 
	{
	  int jmax = imin(i+SELFEXCL,t.L);
	  for (j=1; j<=jmax; j++) 
	    cell_off[i][j]=1;   // cross out cell near diagonal
	  for (j=jmax+1; j<=t.L+1; j++)  
	    cell_off[i][j]=0;   // no other cells crossed out yet
	}
    }
  else 
    // Compare two different HMMs Q and T
    {
      // Activate all cells in dynamic programming matrix
      for (i=1; i<=q.L; i++) 
	for (j=1; j<=t.L; j++) 
	  cell_off[i][j]=0;   // no other cells crossed out yet

      // Cross out cells that are excluded by the minimum-overlap criterion
      if (par.min_overlap==0) 
	min_overlap = imin(60, (int)(0.333f*imin(q.L,t.L))+1); // automatic minimum overlap
      else 
	min_overlap = imin(par.min_overlap, (int)(0.8f*imin(q.L,t.L)));

      for (i=0; i<min_overlap; i++) 
	for (j=i-min_overlap+t.L+1; j<=t.L; j++) // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt} 
	  cell_off[i][j]=1;
      for (i=q.L-min_overlap+1; i<=q.L; i++) 
	for (j=1; j<i+min_overlap-q.L; j++)      // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq} 
	  cell_off[i][j]=1;
    }

  // Cross out rows which are contained in range given by exclstr ("3-57,238-314")
  if (par.exclstr) 
    {
      char* ptr=par.exclstr;
      int i0, i1;
      while (1) 
	{
	  i0 = abs(strint(ptr));
	  i1 = abs(strint(ptr));
	  if (!ptr) break;
	  for (i=i0; i<=imin(i1,q.L); i++) 
	    for (j=1; j<=t.L; j++) 
	      cell_off[i][j]=1; 
	}
    }
}
	
/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Allocate memory for data of new alignment (sequence names, alignment, scores,...)
 */
void 
Hit::InitializeBacktrace(HMM& q, HMM& t)
{
    if (irep==1) //if this is the first single repeat repeat hit with this template
        {
            //Copy information about template profile to hit and reset template pointers to avoid destruction
            longname=new(char[strlen(t.longname)+1])();
            name    =new(char[strlen(t.name)+1])();
            file    =new(char[strlen(t.file)+1])();
            if (!file) {
                MemoryError("space for alignments with database HMMs. \nNote that all alignments have to be kept in memory");
            }
            strcpy(longname,t.longname);
            strcpy(name,t.name);
            strcpy(fam ,t.fam);
            strcpy(sfam ,t.sfam);
            strcpy(fold ,t.fold);
            strcpy(cl ,t.cl);
            strcpy(file,t.file);
            sname=new(char*[t.n_display])();   // Call Compare only once with irep=1
            seq  =new(char*[t.n_display])();   // Call Compare only once with irep=1
            if (!sname || !seq) {
                MemoryError("space for alignments with database HMMs.\nNote that all sequences for display have to be kept in memory");
            }

            for (int k=0; k<t.n_display; k++)	{
                if (NULL != t.sname){
                    sname[k]=t.sname[k]; t.sname[k]=NULL;
                }
                else {
                    sname[k]=NULL;
                }
                seq[k]  =t.seq[k];   t.seq[k]=NULL;
            }
            
            n_display=t.n_display; t.n_display=0;
            ncons  = t.ncons;
            nfirst = t.nfirst;
            nss_dssp = t.nss_dssp;
            nsa_dssp = t.nsa_dssp;
            nss_pred = t.nss_pred;
            nss_conf = t.nss_conf;
            L = t.L;
            Neff_HMM = t.Neff_HMM;
            Eval   = 1.0;
            Pval   = 1.0;
            Pvalt  = 1.0;
            logPval = 0.0;
            logPvalt= 0.0;
            Probab = 1.0;
        }    
    
    // Allocate new space
    this->i = new( int[i2+j2+2])();
    this->j = new( int[i2+j2+2])();
    states  = new( char[i2+j2+2])();
    S = S_ss = P_posterior = NULL; // set to NULL to avoid deleting data from irep=1 when hit with irep=2 is removed 
    Xcons = NULL;
} /* this is the end of Hit::InitializeBacktrace() */

/////////////////////////////////////////////////////////////////////////////////////
// Some score functions 
/////////////////////////////////////////////////////////////////////////////////////


/**
 * @brief Calculate score between columns i and j of two HMMs (query and template)
 */
inline float 
Score(float* qi, float* tj)
{
//   if (par.columnscore==9)
//     return (tj[0] *qi[0] +tj[1] *qi[1] +tj[2] *qi[2] +tj[3] *qi[3] +tj[4]*qi[4]
//            +tj[5] *qi[5] +tj[6] *qi[6] +tj[7] *qi[7] +tj[8] *qi[8] +tj[9]*qi[9]
//            +tj[10]*qi[10]+tj[11]*qi[11]+tj[12]*qi[12]+tj[13]*qi[13]+tj[14]*qi[14]
//            +tj[15]*qi[15]+tj[16]*qi[16]+tj[17]*qi[17]+tj[18]*qi[18]+tj[19]*qi[19]);
//   else
  return fast_log2(
          tj[0] *qi[0] +tj[1] *qi[1] +tj[2] *qi[2] +tj[3] *qi[3] +tj[4] *qi[4]
         +tj[5] *qi[5] +tj[6] *qi[6] +tj[7] *qi[7] +tj[8] *qi[8] +tj[9] *qi[9]
         +tj[10]*qi[10]+tj[11]*qi[11]+tj[12]*qi[12]+tj[13]*qi[13]+tj[14]*qi[14]
         +tj[15]*qi[15]+tj[16]*qi[16]+tj[17]*qi[17]+tj[18]*qi[18]+tj[19]*qi[19]
	  );
}

/**
 * @brief Calculate score between columns i and j of two HMMs (query and template)
 */
inline float 
ProbFwd(float* qi, float* tj)
{
  return  tj[0] *qi[0] +tj[1] *qi[1] +tj[2] *qi[2] +tj[3] *qi[3] +tj[4] *qi[4]
         +tj[5] *qi[5] +tj[6] *qi[6] +tj[7] *qi[7] +tj[8] *qi[8] +tj[9] *qi[9]
         +tj[10]*qi[10]+tj[11]*qi[11]+tj[12]*qi[12]+tj[13]*qi[13]+tj[14]*qi[14]
         +tj[15]*qi[15]+tj[16]*qi[16]+tj[17]*qi[17]+tj[18]*qi[18]+tj[19]*qi[19];
}


/**
 * @brief Calculate secondary structure score between columns i and j of two HMMs (query and template)
 */
inline float 
Hit::ScoreSS(HMM& q, HMM& t, int i, int j, int ssm)
{
  switch (ssm) //SS scoring during alignment 
    {
    case 0: // no SS scoring during alignment 
      return 0.0;
    case 1: // t has dssp information, q has psipred information 
      return par.ssw * S73[ (int)t.ss_dssp[j]][ (int)q.ss_pred[i]][ (int)q.ss_conf[i]];
    case 2: // q has dssp information, t has psipred information 
      return par.ssw * S73[ (int)q.ss_dssp[i]][ (int)t.ss_pred[j]][ (int)t.ss_conf[j]];
    case 3: // q has dssp information, t has psipred information 
      return par.ssw * S33[ (int)q.ss_pred[i]][ (int)q.ss_conf[i]][ (int)t.ss_pred[j]][ (int)t.ss_conf[j]];
//     case 4: // q has dssp information, t has dssp information 
//       return par.ssw*S77[ (int)t.ss_dssp[j]][ (int)t.ss_conf[j]];
    }
  return 0.0;
}

/**
 * @brief Calculate secondary structure score between columns i and j of two HMMs (query and template)
 */
inline float 
Hit::ScoreSS(HMM& q, HMM& t, int i, int j)
{
  return ScoreSS(q,t,i,j,ssm2);
}


/**
 * @brief Calculate score between columns i and j of two HMMs (query and template)
 */
inline float 
Hit::ScoreTot(HMM& q, HMM& t, int i, int j)
{
  return Score(q.p[i],t.p[j]) + ScoreSS(q,t,i,j) + par.shift;
}

/*
 * Calculate score between columns i and j of two HMMs (query and template)
 */
inline float 
Hit::ScoreAA(HMM& q, HMM& t, int i, int j)
{
  return Score(q.p[i],t.p[j]);
}


/////////////////////////////////////////////////////////////////////////////////////
/*
 * Function for Viterbi()
 */
inline float 
max2(const float& xMM, const float& xX, char& b) 
{
  if (xMM>xX) { b=MM; return xMM;} else { b=SAME;  return xX;}
}


/////////////////////////////////////////////////////////////////////////////////////
/*
 * Functions for StochasticBacktrace()
 */

inline int 
pickprob2(const double& xMM, const double& xX, const int& state) 
{
  if ( (xMM+xX)*frand() < xMM) return MM; else return state; 
}

inline int 
pickprob3_GD(const double& xMM, const double& xDG, const double& xGD) 
{
  double x = (xMM+xDG+xGD)*frand();
  if ( x<xMM) return MM; 
  else if ( x<xMM+xDG) return DG; 
  else return GD;
}

inline int 
pickprob3_IM(const double& xMM, const double& xMI, const double& xIM) 
{
  double x = (xMM+xMI+xIM)*frand();
  if ( x<xMM) return MM; 
  else if ( x<xMM+xMI) return MI; 
  else return IM;
}

inline int 
pickprob6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI) 
{
  double x = (x0+xMM+xGD+xIM+xDG+xMI)*frand();
  x-=xMM; if (x<0) return MM; 
  x-=x0;  if (x<0) return STOP; 
  x-=xGD; if (x<0) return GD;
  x-=xIM; if (x<0) return IM;
  if (x < xDG) return DG; else return MI;
}

inline int 
pickmax2(const double& xMM, const double& xX, const int& state) 
{
  if (xMM > xX) return MM; else return state; 
}

inline int 
pickmax3_GD(const double& xMM, const double& xDG, const double& xGD) 
{
  char state;
  double x;
  if ( xMM>xDG) {state=MM; x=xMM;} 
  else          {state=DG; x=xDG;}
  if ( xGD>x)   {state=GD; x=xGD;}
  return state;
}

inline int 
pickmax3_IM(const double& xMM, const double& xMI, const double& xIM) 
{
  char state;
  double x;
  if ( xMM>xMI) {state=MM; x=xMM;}
  else          {state=MI; x=xMI;}
  if ( xIM>x)   {state=IM; x=xIM;}
  return state;
}

inline int 
pickmax6(const double& x0, const double& xMM, const double& xGD, const double& xIM, const double& xDG, const double& xMI) 
{
  char state;
  double x;
  if ( x0 >xMM) {state=STOP; x=x0;} 
  else          {state=MM; x=xMM;}
  if ( xGD>x)   {state=GD; x=xGD;}
  if ( xIM>x)   {state=IM; x=xIM;}
  if ( xDG>x)   {state=DG; x=xDG;}
  if ( xMI>x)   {state=MI; x=xMI;}
  return state;
}


/////////////////////////////////////////////////////////////////////////////////////
//// Functions that calculate P-values and probabilities 
/////////////////////////////////////////////////////////////////////////////////////


//// Evaluate the CUMULATIVE extreme value distribution at point x
//// p(s)ds = lamda * exp{ -exp[-lamda*(s-mu)] - lamda*(s-mu) } ds = exp( -exp(-x) - x) dx = p(x) dx
//// => P(s>S) = integral_-inf^inf {p(x) dx}  = 1 - exp{ -exp[-lamda*(S-mu)] }
inline double 
Pvalue(double x, double a[])
{
  //a[0]=lamda, a[1]=mu
  double h = a[0]*(x-a[1]);
  return (h>10)? exp(-h) : double(1.0)-exp( -exp(-h));
}

inline double 
Pvalue(float x, float lamda, float mu)
{
  double h = lamda*(x-mu);
  return (h>10)? exp(-h) : (double(1.0)-exp( -exp(-h)));
}

inline double 
logPvalue(float x, float lamda, float mu)
{
  double h = lamda*(x-mu);
  return (h>10)? -h : (h<-2.5)? -exp(-exp(-h)): log( ( double(1.0) - exp(-exp(-h)) ) );
}

inline double 
logPvalue(float x, double a[])
{
  double h = a[0]*(x-a[1]);
  return (h>10)? -h : (h<-2.5)? -exp(-exp(-h)): log( ( double(1.0) - exp(-exp(-h)) ) );
}

// Calculate probability of true positive : p_TP(score)/( p_TP(score)+p_FP(score) )
// TP: same superfamily OR MAXSUB score >=0.1
inline double 
Probab(Hit& hit)
{
  double s=-hit.score_aass;
  double t;
  if (s>200) return 100.0; 
  if (par.loc) 
    {
      if (par.ssm && (hit.ssm1 || hit.ssm2) && par.ssw>0) 
	{
	  // local with SS
	  const double a=sqrt(6000.0);
	  const double b=2.0*2.5;
	  const double c=sqrt(0.12);
	  const double d=2.0*32.0;
	  t = a*exp(-s/b) + c*exp(-s/d);
	}
      else
	{
	  // local no SS
	  const double a=sqrt(4000.0);
	  const double b=2.0*2.5;
	  const double c=sqrt(0.15);
	  const double d=2.0*34.0;
	  t = a*exp(-s/b) + c*exp(-s/d);
	}
    }
  else
    {
      if ( (par.ssm>0) && (par.ssw>0) ) /* FIXME: was '&', should be '&&' (or not?) */
	{
	  // global with SS
	  const double a=sqrt(4000.0);
	  const double b=2.0*3.0;
	  const double c=sqrt(0.13);
	  const double d=2.0*34.0;
	  t = a*exp(-s/b) + c*exp(-s/d);
	}
      else
	{
	  // global no SS
	  const double a=sqrt(6000.0);
	  const double b=2.0*2.5;
	  const double c=sqrt(0.10);
	  const double d=2.0*37.0;
	  t = a*exp(-s/b) + c*exp(-s/d);
	}

    }

  return 100.0/(1.0+t*t);
}

// #define Weff(Neff) (1.0+par.neffa*(Neff-1.0)+(par.neffb-4.0*par.neffa)/16.0*(Neff-1.0)*(Neff-1.0))

// /////////////////////////////////////////////////////////////////////////////////////
// // Merge HMM with next aligned HMM  
// /////////////////////////////////////////////////////////////////////////////////////
// void Hit::MergeHMM(HMM& Q, HMM& t, float wk[])
// {
//   int i,j;    // position in query and target
//   int a;      // amino acid
//   int step;   // alignment position (step=1 is end)
//   float Weff_M, Weff_D, Weff_I;
//   for (step=nsteps; step>=2; step--) // iterate only to one before last alignment column
//     {
//       i = this->i[step];
//       j = this->j[step];
//       switch(states[step])
// 	{
// 	case MM: 
// 	  Weff_M = Weff(t.Neff_M[j]-1.0);
// 	  Weff_D = Weff(t.Neff_D[j]-1.0);
// 	  Weff_I = Weff(t.Neff_I[j]-1.0);
// 	  for (a=0; a<20; a++) Q.f[i][a] += t.f[j][a]*wk[j]*Weff_M;
// 	  switch(states[step-1])
// 	    {
// 	    case MM:  // MM->MM
// 	      Q.tr_lin[i][M2M]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2D]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2I]+= t.tr_lin[j][M2I]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][D2M]+= t.tr_lin[j][D2M]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2D]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      break; 
// 	    case MI: // MM->MI
// 	      Q.tr_lin[i][M2D]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2D]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2M]+= t.tr_lin[j][M2I]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2M]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2D]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      break;
// 	    case DG: // MM->DG
// 	      Q.tr_lin[i][M2D]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2D]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2M]+= t.tr_lin[j][M2I]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2M]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2D]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      break;
// 	    case IM: // MM->IM
// 	      Q.tr_lin[i][M2I]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2M]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2I]+= t.tr_lin[j][M2I]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][D2M]+= t.tr_lin[j][D2M]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2D]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      break;
// 	    case GD: // MM->GD
// 	      Q.tr_lin[i][M2I]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2M]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][M2I]+= t.tr_lin[j][M2I]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][D2M]+= t.tr_lin[j][D2M]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2D]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      break;
// 	    }
// 	  break;

// 	case MI: // if gap in template  
// 	  Weff_I = Weff(t.Neff_I[j]-1.0);
// 	  switch(states[step-1])
// 	    {
// 	    case MI:  // MI->MI
// 	      Q.tr_lin[i][M2M]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      break;
// 	    case MM:  // MI->MM
// 	      Q.tr_lin[i][M2M]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      break;
// 	    }
// 	  break;

// 	case DG:   
// 	  Weff_M = Weff(t.Neff_M[j]-1.0);
// 	  Weff_D = Weff(t.Neff_D[j]-1.0);
// 	  Weff_I = Weff(t.Neff_I[j]-1.0);
// 	  switch(states[step-1])
// 	    {
// 	    case DG:  // DG->DG
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2M]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2D]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][M2M]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][M2D]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      break;
// 	    case MM:  // DG->MM
// 	      Q.tr_lin[i][D2M]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][D2M]+= t.tr_lin[j][D2M]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2D]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][M2M]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      break;
// 	    }
// 	  break;
	  
// 	case IM: // if gap in query  
// 	  Weff_M = Weff(t.Neff_M[j]-1.0);
// 	  switch(states[step-1])
// 	    {
// 	    case IM:  // IM->IM
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][M2I]*wk[j]*Weff_M;
// 	      break;
// 	    case MM:  // IM->MM
// 	      Weff_D = Weff(t.Neff_D[j]-1.0);
// 	      Weff_I = Weff(t.Neff_I[j]-1.0);
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][M2I]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][D2M]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2D]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      break;
// 	    }
// 	  break;
	  
// 	case GD:   
// 	  Weff_M = Weff(t.Neff_M[j]-1.0);
// 	  switch(states[step-1])
// 	    {
// 	    case GD:  // GD->GD
// 	      Weff_I = Weff(t.Neff_I[j]-1.0);
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][M2I]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      break;
// 	    case MM:  // GD->MM
// 	      Weff_D = Weff(t.Neff_D[j]-1.0);
// 	      Weff_I = Weff(t.Neff_I[j]-1.0);
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][M2M]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][M2D]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][M2I]*wk[j]*Weff_M;
// 	      Q.tr_lin[i][D2D]+= t.tr_lin[j][D2D]*wk[j]*Weff_D;
// 	      Q.tr_lin[i][I2M]+= t.tr_lin[j][I2M]*wk[j]*Weff_I;
// 	      Q.tr_lin[i][I2I]+= t.tr_lin[j][I2I]*wk[j]*Weff_I;
// 	      break;
// 	    }
// 	  break;

// 	}
//     }
//   i = this->i[step];
//   j = this->j[step];
//   Weff_M = Weff(t.Neff_M[j]-1.0);
//   for (a=0; a<20; a++) Q.f[i][a] += t.f[j][a]*wk[j]*Weff_M;
// }


#ifdef CLUSTALO
/* @* Hit::ClobberGlobal (eg, hit)
 *
 */
void 
Hit::ClobberGlobal(void){

    if (i){
      //delete[] i; 
      i = NULL;
    }
    if (j){
      //delete[] j; 
      j = NULL;
    }
    if (states){
      //delete[] states; 
      states = NULL;
    }
    if (S){
      //delete[] S; 
      S = NULL;
    }
    if (S_ss){
      //delete[] S_ss; 
      S_ss = NULL;
    }
    if (P_posterior){
      //delete[] P_posterior; 
      P_posterior = NULL;
    }
    if (Xcons){
      //delete[] Xcons; 
      Xcons = NULL;
    }
    //  delete[] l;    l = NULL;
    i = j = NULL;
    states = NULL;
    S = S_ss = P_posterior = NULL;
    Xcons = NULL;
    if (irep==1) // if irep>1 then longname etc point to the same memory locations as the first repeat. 
      {          // but these have already been deleted.
	// 	printf("Delete name = %s\n",name);//////////////////////////
	//delete[] longname; 
	longname = NULL;
	//delete[] name; 
	name = NULL;
	//delete[] file; 
	file = NULL;
	//delete[] dbfile; 
	dbfile = NULL;
	/*for (int k=0; k<n_display; k++)
	  {
	  delete[] sname[k]; sname[k] = NULL;
	  delete[] seq[k]; seq[k] = NULL;
	  }*/
	//delete[] sname; 
	sname = NULL;
	//delete[] seq; 
	seq = NULL;
      }

    score = score_sort = score_aass = 0.0;
    Pval = Pvalt = Eval = Probab = 0;
    Pforward = sum_of_probs = 0.00;
    L = irep = nrep = n_display = nsteps = 0;
    i1 = i2 = j1 = j2 = matched_cols = min_overlap = 0;
}
#endif


/*
 * EOF hhhit-C.h
 */
