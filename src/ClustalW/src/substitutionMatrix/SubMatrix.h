/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * This is the interface class to all the substitution matrices.
 * It provides the matrices in a form that the rest of the program can use.
 * It is also used to store the user defined matrix. This will be used mainly as an interface
 * to the matrices defined in matrices.h.
 * The way this class will work is the user can read in matrix series or a single matrix, 
 * or they can select one of the matrix series (e.g Blosum). This will then be used in the
 * alignment stages. There are separate matrices for amino acid pairwise and progressive, 
 * and for DNA alignments both pairwise and progressive.
 * It is possible to have a series of matrices that are user defined for amino acid
 * progressive ONLY!!
 * A single matrix is choosen for pairwise and for DNA alignments.
 * This class does 3 jobs. It reads in matrices from files/arrays, it provides
 * matrices in the formats that are used in the alignment stage, it allows users to select
 * which matrices they would like to use.
 */

#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include <vector>
#include <string>
#include "../general/clustalw.h"
#include "../general/userparams.h"
#include "../general/utils.h"
#include "../general/Array2D.h" 
#include <Rcpp.h>

namespace clustalw
{
using namespace std;

typedef vector<short> Xref;
typedef vector<short> Matrix;

class SubMatrix
{
    public:
        /* Functions */
        SubMatrix();
        ~SubMatrix();
        
        bool getUserMatFromR(Rcpp::NumericMatrix substitutionMatrix, int alignResidueType, int alignType);
        bool getAAScoreMatFromR(Rcpp::NumericMatrix substitutionMatrix);
        bool getDNAScoreMatFromR(Rcpp::NumericMatrix substitutionMatrix);
        bool getQTLowScoreMatFromR(Rcpp::NumericMatrix substitutionMatrix, bool dna);
        bool getUserMatSeriesFromR(Rcpp::NumericMatrix substitutionMatrix);

        bool getUserMatFromFile(char *str, int alignResidueType, int alignType);
		bool getAAScoreMatFromFile(char *str);
		bool getDNAScoreMatFromFile(char *str);
		bool getQTLowScoreMatFromFile(char *str, bool dna);
		bool getUserMatSeriesFromFile(char *str);

        void setCurrentNameAndNum(string _matrixName, int _matrixNum, int alignResidueType,
                                  int alignType);
        int getMatrixNumForMenu(int alignResidueType, int alignType);
        int getPairwiseMatrix(int matrix[NUMRES][NUMRES], PairScaleValues& scale, 
                              int& matAvg);
        int getProfileAlignMatrix(int matrix[NUMRES][NUMRES], double pcid, int minLen, 
                                  PrfScaleValues& scaleParam, int& matAvg);
        int getAlnScoreMatrix(int matrix[NUMRES][NUMRES]);
        // Access functions for the interactive menu. 
        int getMatrixNum();
        int getDNAMatrixNum();
        int getPWMatrixNum();
        int getPWDNAMatrixNum();
        void getQTMatrixForHistogram(int matrix[NUMRES][NUMRES]); 
                                    // NOTE Qt
        int getQTAAHistMatNum(){return QTAAHistMatNum;};
        int getQTDNAHistMatNum(){return QTDNAHistMatNum;};
        void setQTAAHistMatNum(int num){QTAAHistMatNum = num;};
        void setQTDNAHistMatNum(int num){QTDNAHistMatNum = num;};

        void getQTMatrixForLowScoreSeg(int matrix[NUMRES][NUMRES]);
        int getQTsegmentDNAMatNum(){return QTsegmentDNAMatNum;}
        void setQTsegmentDNAMatNum(int dnaMat){QTsegmentDNAMatNum = dnaMat;}
        int getQTsegmentAAMatNum(){return QTsegmentAAMatNum;}
        void setQTsegmentAAMatNum(int aaMat){QTsegmentAAMatNum = aaMat;}
        
        void tempInterface(int alignResidueType, int alignType);
        void setValuesToDefault();
        /* Attributes */

    private:
        
        /* Functions */
        int getMatrix(Matrix* matPtr, Xref* xref, int matrix[NUMRES][NUMRES],
                      bool negFlag, int scale, bool minimise = false); 
        int readMatrixSeriesFromR(const Rcpp::NumericMatrix substitutionmatrix, Matrix& userMat, Xref& xref);
        int readUserMatrixFromR(const Rcpp::NumericMatrix substitutionMatrix, Matrix& userMat, Xref& xref);
        int readMatrixSeries(const char *fileName, Matrix& userMat, Xref& xref);
        int readUserMatrix(const char *fileName, Matrix& userMat, Xref& xref);
        int getArgs(string line, char *args[], int max);
        void setUpCrossReferences();
        bool commentline(char* line);
        
        // The functions below are purely for testing purposes.
        void printGetMatrixResults(int mat[NUMRES][NUMRES]); 
        void compareMatrices(int mat1[NUMRES][NUMRES], int mat2[NUMRES][NUMRES]); 
        void printInFormat(vector<short>& temp, const char* name = "tempfile.out");
        void printVectorToFile(vector<short>& temp, const char* name = "tempfile.out");
        Matrix* getUserMatAddress(int alignResidueType, int alignType);
        Xref* getUserXrefAddress(int alignResidueType, int alignType);
        void checkResidueAndAlignType(int alignResidueType, int alignType);
                                     
        /* Attributes */
        bool userSeries;
        int matrixNum;
        int DNAMatrixNum;
        int pwMatrixNum;
        int pwDNAMatrixNum;
        
        string* matrixName;
        string* DNAMatrixName;
        string* pwMatrixName;
        string* pwDNAMatrixName;
                          
        // Matrix cross references.
        Xref defaultDNAXref;
        Xref defaultAAXref;
        Xref DNAXref; // User defined dna xref
        Xref AAXref;
        Xref pwAAXref; // pairwise
        Xref pwDNAXref;
        Xref QTscoreXref;
        Xref QTscoreDNAXref;
        Xref QTsegmentDNAXref;
        Xref QTsegmentAAXref;
        vector<Xref> AAXrefseries;
        
        vector<Matrix> userMatSeries;
        Matrix userMat;
        Matrix pwUserMat;
        Matrix userDNAMat;
        Matrix pwUserDNAMat;
        Matrix QTscoreUserMatrix;
        Matrix QTscoreUserDNAMatrix;
        Matrix QTsegmentDNAMatrix;
        Matrix QTsegmentAAMatrix;
        
        /* These are vectors to store the matrices defined in matrices.h */
        const int sizenAAMatrix; 
        const int sizeDNAMatrix; 
        Matrix* blosum30mtVec; 
        Matrix* blosum40mtVec;
        Matrix* blosum45mtVec;
        Matrix* blosum62mt2Vec;
        Matrix* blosum80mtVec;
        Matrix* pam20mtVec;
        Matrix* pam60mtVec;
        Matrix* pam120mtVec;
        Matrix* pam350mtVec;
        Matrix* idmatVec;
        Matrix* gon40mtVec;
        Matrix* gon80mtVec;
        Matrix* gon120mtVec;
        Matrix* gon160mtVec;
        Matrix* gon250mtVec;
        Matrix* gon350mtVec;
        Matrix* clustalvdnamtVec;
        Matrix* swgapdnamtVec;
        
        int matrixAvgScore; // NOTE Needed by other classes.
        UserMatrixSeries matSeries;
        string line2;
                
        int QTDNAHistMatNum;
        int QTAAHistMatNum;
        int QTsegmentDNAMatNum;
        int QTsegmentAAMatNum;
        
        // Temp, to hold current selection
        Matrix* mat;
        Xref* xref;
        Matrix* _matPtr;
        Xref* _matXref;

};
}
#endif

