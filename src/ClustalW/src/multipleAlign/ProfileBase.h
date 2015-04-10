/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * The reason why we have 2 different types of profiles is because one of them 
 * (ProfileWithSub) has the substitution matrix information already in it. This
 * increases the performance when aligning 2 profile columns.
 */
#ifndef PROFILEBASE_H
#define PROFILEBASE_H

#include "../alignment/Alignment.h"

namespace clustalw
{

class ProfileBase
{
    public:
        /* Functions */
        ProfileBase(int prfLen, int firstS, int lastS);
        void calcGapCoeff(SeqArray* seqArray, vector<int>* gaps,  bool useStructPenalties,
                          vector<char>* gapPenaltyMask, int gapCoef, int lenCoef);
        const SeqArray* getProfilePtr(){return &profile;};
        void resetProfile(){for(int i = 0; i < (int)profile.size();i++)
                            {
                                profile[i].clear();
                            }
                            profile.clear();
                            };
        /* Attributes */

    protected:
        /* Functions */
        void calcVPenalties(SeqArray* aln, vector<int>* weight); 
        void calcResidueSpecificPen(SeqArray* aln, vector<int>* weight); 
        void calcHydrophilicPen(SeqArray* aln, vector<int>* weight); 
        int localPenalty(int penalty, int n, vector<int>* resWeight, vector<int>* hydWeight,
                         vector<int>* vWeight);  
        float percentId(vector<int>* s1, vector<int>* s2);

        /* Attributes */
        vector<vector<int> > profile;

        /* number of residues used for a window for the variable zone
           penalties */
        int vwindow;
        /* vll is the lower limit for the variable zone penalties
           (vll < pen < 1.0) */
        int vll;

        /* "Pascarella and Argos" residue specific gap modification
           factors. See Table 1 in the ClustalW 1994 NAR paper
           http://www.ncbi.nlm.nih.gov/pubmed/7984417
        */
        string pascarellaRes;
        vector<int> pascarellaProb;

        vector<vector<int> > vlut;
        static const int numLetters = 26;
        float reducedGap;
        bool nVarPen;
        bool nHydPen;
        bool nPrefPen;
        int gdist;
        int prfLength;
        int firstSeq, lastSeq;
    private:
        /* Functions */

        /* Attributes */

};

}
#endif
