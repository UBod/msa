/**
 *
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <stdio.h>
#include "../alignment/Alignment.h"
#ifdef HAVE_MHASH_H
#include "mhash.h"
#endif

#if WIN32
	#include <time.h>
#endif

#include "Stats.h"


using namespace std;

namespace clustalw
{


Stats::Stats()
{
    enabled = false;
}


Stats::~Stats()
{
}




/* adopted from Sean Eddy'ssquid:aligneval.c */
/* Function: PairwiseIdentity()
 * 
 * Purpose:  Calculate the pairwise fractional identity between
 *           two aligned sequences s1 and s2. This is simply
 *           (idents / MIN(len1, len2)).
 *
 *           Note how many ways there are to calculate pairwise identity,
 *           because of the variety of choices for the denominator:
 *           idents/(idents+mismat) has the disadvantage that artifactual
 *             gappy alignments would have high "identities".
 *           idents/(AVG|MAX)(len1,len2) both have the disadvantage that 
 *             alignments of fragments to longer sequences would have
 *             artifactually low "identities".
 *           
 *           Case sensitive; also, watch out in nucleic acid alignments; 
 *           U/T RNA/DNA alignments will be counted as mismatches!
 */
/* float
PairwiseIdentity(char @s1, char @s2)
{
  int     idents;		/@ total identical positions  @/
  int     len1, len2;		/@ lengths of seqs            @/
  int     x;			/@ position in aligned seqs   @/

  idents = len1 = len2 = 0;
  for (x = 0; s1[x] != '\0' && s2[x] != '\0'; x++) 
  {
      if (!isgap(s1[x])) {
          len1++;
          if (s1[x] == s2[x]) idents++; 
      }
      if (!isgap(s2[x])) len2++;
  }
  if (len2 < len1) len1 = len2;
  return (len1 == 0 ? 0.0 : (float) idents / (float) len1);
}
*/


/* s1/s2 are unit-offset as usual
 *
 */
float
Stats::pairwiseIdentity(Alignment *alnObj, int s1, int s2)
{
    int idents;     /* total identical positions */
    int len1, len2; /* real lengths of seqs      */
    int x;          /* position in aligned seqs  */
    
    const vector<int>* seq1 = alnObj->getSequence(s1);
    const vector<int>* seq2 = alnObj->getSequence(s2);
    
    idents = len1 = len2 = 0;
    // cerr << "comparing " << alnObj->getName(s1).c_str() << ":" << alnObj->getName(s2).c_str() << " " << s1 << ":" << s2  << "\n";
    
    // sequence length should be identical, but be paranoid
    for (x = 1; x<=alnObj->getSeqLength(s1) && x<=alnObj->getSeqLength(s2); x++)
    {
        if (! alnObj->isGap(s1, x)) {
            len1++;
            //cerr << " pos " << x << ": " << (*seq1)[x] << ":" << (*seq2)[x] << "\n";
            if ((*seq1)[x] == (*seq2)[x])
                idents++; 
        }
        //DEBUG
        //else {
        //cerr << " gap at pos " << x << " (" << s1 << ")\n";
        //}
        if (! alnObj->isGap(s2, x))
            len2++;
        //DEBUG
        //else
        //    cerr << " gap at pos " << x << " (" << s2 << ")\n";
    }
    if (len2 < len1)
        len1 = len2;
    return (len1 == 0 ? 0.0 : (float) idents / (float) len1);
}
  


#ifdef HAVE_MHASH_H

string
Stats::ConcatInputHash(Alignment *alnObj)
{
    vector<string> rawSeqArray;
    string ret;
    char *hash;

    // collect all sequences and sort 
    const clustalw::SeqArray* seqArray = alnObj->getSeqArray();
    for(int s = 1; s <= alnObj->getNumSeqs(); s++)
    {
        string seq;
        for(int r = 1; r <= alnObj->getSeqLength(s); r++)
        {
            int val = (*seqArray)[s][r];
            seq.append(1, clustalw::userParameters->getAminoAcidCode(val));
        }
        rawSeqArray.push_back(seq);
    }
    std::sort(rawSeqArray.begin(), rawSeqArray.end());

    // concatenate sorted seqs
    string concatSeq;
    std::vector<string>::iterator iter;
    for(iter=rawSeqArray.begin(); iter != rawSeqArray.end(); ++iter)
        concatSeq.append(*iter);

    // build hash and return
    hash = Md5Hash(concatSeq.c_str());
    if (hash==NULL)
    {
        ret="HASHING_FAILURE";
    } else {
        for (int i=0; i<strlen(hash); i++)
        ret.append(1, hash[i]);
        free(hash);
    }
    return ret;
}


string
Stats::Md5ForSeq(Alignment *alnObj, int s)
{
    string seq;
    const clustalw::SeqArray* seqArray = alnObj->getSeqArray();
    char *hash;
    string ret;
    
    for(int l = 1; l <= alnObj->getSeqLength(s); l++)
    {
        int val = (*seqArray)[s][l];
        // continue if gap
        if((val < 0) || (val > userParameters->getMaxAA()))
            continue;
        seq.append(1, clustalw::userParameters->getAminoAcidCode(val));
    }
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

    hash = Md5Hash(seq.c_str());
    if (hash==NULL)
    {
        ret = "HASHING_FAILURE";
    } else {
        ret = hash;
    }
    return ret;
}

/* create md5 hash from input string
 * returns NULL on failure
 * user must free returned string
 */
char *
Stats::Md5Hash(const char *thread)
{
    MHASH td;
    char *tmpcstr;
    int i;
    unsigned char *hash;
    char *rethash;
    
    td = mhash_init(MHASH_MD5);
    if (td == MHASH_FAILED)
        return NULL;

    if (thread==NULL)
        return NULL;

    mhash(td, thread, strlen(thread));
    //mhash_deinit(td, hash);
    hash = (unsigned char*) mhash_end(td);
    
    rethash = (char*) calloc(mhash_get_block_size(MHASH_MD5)*2+1, sizeof(char));
    for (i = 0; i < mhash_get_block_size(MHASH_MD5); i++) {
        sprintf(&rethash[i*2], "%.2x", hash[i]);
    }

    return rethash;
}

#endif
// HAVE_MHASH_H




void
Stats::logCmdLine(int argc, char **argv)
{
    FILE *fp = fopen(logfilename.c_str(), "a");
    if (fp == NULL)
    {
        cerr << "couldn't open file " << logfilename << " for logging of stats\n";
        return;
    }

    if (argc > 1)
    {
        for (int i=1; i<argc; i++)
        {
            // remove non-interesting stuff
            if (strstr(argv[i], "-infile=")==NULL &&
                strstr(argv[i], "-outfile=")==NULL &&
                strstr(argv[i], "-stats=")==NULL &&
                strstr(argv[i], "-align")==NULL)
            {
                fprintf(fp, "cmdline non-default arg: %s\n", argv[i]);
            }
        }
    }
}


/* log some statistics for input sequences, i.e. alnObj should hold
 * the unaligned sequences
 *
 */
void
Stats::logInputSeqStats(Alignment *alnObj)
{
    int i;
    std::vector<double> lengths;
    time_t t = time(NULL);
    tm s = *localtime(&t);
    int shortest;
    string hash;
    
    FILE *fp = fopen(logfilename.c_str(), "a");
    if (fp == NULL)
    {
        cerr << "couldn't open file " << logfilename << " for logging of stats\n";
        return;
    }

    fprintf(fp, "logging job: %s on %s", userParameters->getSeqName().c_str(), asctime(&s));
    fprintf(fp, "clustal version: %s\n", userParameters->getRevisionLevel().c_str());

    fprintf(fp, "seq type: ");
    if (userParameters->getDNAFlag())
        fprintf(fp, "DNA");
    else
        fprintf(fp, "protein");
    fprintf(fp, "\n");
    
    fprintf(fp, "numseqs: %d\n", alnObj->getNumSeqs());

    // create a vector of seq lengths for later
    // and get shortest seq at the same time
    shortest=alnObj->getLengthLongestSequence();
    for (i = 1; i <= alnObj->getNumSeqs(); i++) {
        int l = alnObj->getSeqLength(i);
        lengths.push_back(l);
        if (l<shortest)
            shortest=l;
    }
    
    fprintf(fp, "seqlen longest: %d\n", alnObj->getLengthLongestSequence());
    fprintf(fp, "seqlen shortest: %d\n", shortest);

    fprintf(fp, "seqlen avg: %.2f\n", utilityObject->average(lengths));
    fprintf(fp, "seqlen std-dev: %.2f\n", utilityObject->stdDev(lengths));
    fprintf(fp, "seqlen median: %.2f\n", utilityObject->median(lengths));

    
#ifdef HAVE_MHASH_H
    //hash = concatInputHash(alnObj);
    //fprintf(fp, "seq hash: %s\n", hash.c_str());
    
    for (int s = 1; s <= alnObj->getNumSeqs(); s++) {
        string md5 = Md5ForSeq(alnObj, s);
        fprintf(fp, "md5 for seq %d: %s\n", s, md5.c_str());

    }
#else
    fprintf(fp, "md5: disabled\n");    
#endif

    fclose(fp);
}



/* log some statistics for aligned sequences, i.e. alnObj should hold
 * the already aligned sequences
 *
 */
void
Stats::logAlignedSeqStats(Alignment *alnObj)
{
    FILE *fp = fopen(logfilename.c_str(), "a");
    if (fp == NULL)
    {
        cerr << "couldn't open file " << logfilename << " for logging of stats\n";
        return;
    }
    
    // alignment length is the length of any sequence
    fprintf(fp, "aln len: %d\n", alnObj->getSeqLength(1));


    std::vector<double> pwIdents;
    double lowestPwId = 1.0;
    double hightestPwId = 0.0;

    // create vector of pairwise identities
    for(int s1 = 1; s1 <= alnObj->getNumSeqs(); s1++)
    {        
        for(int s2 = s1+1; s2 <= alnObj->getNumSeqs(); s2++)
        {
            double thisPwId = pairwiseIdentity(alnObj, s1, s2);
            pwIdents.push_back(thisPwId);
            if (thisPwId>hightestPwId)
                hightestPwId=thisPwId;
            if (thisPwId<lowestPwId)
                lowestPwId=thisPwId;
            //fprintf(fp, "ident %s:%s %d:%d=%f\n", alnObj->getName(s1).c_str(), alnObj->getName(s2).c_str(), s1, s2, PairwiseIdentity(alnObj, s1, s2));
        }
    }
    fprintf(fp, "aln pw-id highest: %.2f\n", hightestPwId);
    fprintf(fp, "aln pw-id lowest: %.2f\n", lowestPwId);
    fprintf(fp, "aln pw-id avg: %.2f\n",  utilityObject->average(pwIdents));
    fprintf(fp, "aln pw-id std-dev: %.2f\n",  utilityObject->stdDev(pwIdents));
    fprintf(fp, "aln pw-id median: %.2f\n",  utilityObject->median(pwIdents));

    fclose(fp);
}

}
    
