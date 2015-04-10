/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 *
 * Changes:
 *
 * 13-02-07,Nigel Brown(EMBL): Increased maximum sequence identifier
 * width MAXNAMES from 30 to 150.
 * 20-12-07, Paul McGettigan: bug #53 change MAXNAMESTODISPLAY back to 10 from 30
 */

#ifndef CLUSTALW_H
#define CLUSTALW_H
/* Mark tidy up Nov 2005 */
/*********************CLUSTALW.H*********************************************/
/****************************************************************************/

/*
 *** AW NOT NEEDED ANYMORE since 2.0.9
 *** done via including config.h(clustalw) or DEFINE(clustalx)
 * Main header file for ClustalW.  Uncomment ONE of the following 4 lines
 * depending on which compiler you wish to use.
 */
/* NOT SUPPORTED #define VMS 1                 VAX or ALPHA VMS */
/* Think_C for Macintosh */
//#define MAC 1 */
/* Turbo C for PC's */
// #define WINDOWS 1
/* unix */
//#define UNIX 1
// d


#define DEBUGFULL 0
const bool DEBUGLOG = false;
/***************************************************************************/
/***************************************************************************/


#include "general.h"
#include "Array2D.h"
#include "SquareMat.h"
#include "SymMatrix.h"
#include <vector>
#include <string>

#include <Rcpp.h>
using namespace std;
namespace clustalw
{

typedef SymMatrix DistMatrix;
typedef std::vector<vector <int> > TreeGroups;

struct TreeNames
{
    string phylipName;
    string clustalName;
    string distName;
    string nexusName;
    string pimName;
};

struct AlignmentFileNames
{
    string treeFile;
    string profile2TreeFile;
    string clustalFile; 
    string nrbfFile;
    string gcgFile;
    string phylipFile;
    string gdeFile;
    string nexusFile;
    string fastaFile;
};

struct TreeNode
{
    // phylogenetic tree structure
    struct TreeNode *left;
    struct TreeNode *right;
    struct TreeNode *parent;
    float dist;
    int leaf;
    int order;
    string name;
};

struct PhyloTree
{
    TreeGroups treeDesc;
    vector<double> leftBranch;
    vector<double> rightBranch;
};
struct SeqInfo
{
    int firstSeq;
    int lastSeq;
    int numSeqs;
};

struct LowScoreSegParams
{
    int firstSeq; 
    int nSeqs;
    int lastSeq;
    int nCols;
    vector<int>* seqWeight;
    Array2D<int>* lowScoreRes;
    bool seqWeightCalculated;
};
/* Global constants */
const int extraEndElemNum = 2;
const int ENDALN = 127;
const int OK = -200;
const int CANNOTOPENFILE = -300;
const int NOSEQUENCESINFILE = -400;
const int OTHERERROR = -500;
const int ALLNAMESNOTDIFFERENT = -600;
const int MUSTREADINPROFILE1FIRST = -700;
const int EMPTYSEQUENCE = -800;
const int SEQUENCETOOBIG = -900;
const int BADFORMAT = -1000;
     
const int AABLOSUM = 0;
const int AAPAM = 1;
const int AAGONNET = 2;
const int AAIDENTITY = 3;
const int AAUSERDEFINED = 4;

const int PWAABLOSUM = 0;
const int PWAAPAM = 1;
const int PWAAGONNET = 2;
const int PWAAIDENTITY = 3;
const int PWAAUSER = 4;

const int DNAIUB = 0;
const int DNACLUSTALW = 1;
const int DNAUSERDEFINED = 2;

const int AAHISTIDENTITY = 0;
const int AAHISTGONNETPAM80 = 1;
const int AAHISTGONNETPAM120 = 2;
const int AAHISTGONNETPAM250 = 3;
const int AAHISTGONNETPAM350 = 4;
const int AAHISTUSER = 5;

const int QTAASEGGONNETPAM80 = 0;
const int QTAASEGGONNETPAM120 = 1;
const int QTAASEGGONNETPAM250 = 2;
const int QTAASEGGONNETPAM350 = 3;
const int QTAASEGUSER = 4;

const int MAXHYDRESIDUES = 9; // Only allowing 9 hyd residue choices
const int Protein = 0;
const int DNA = 1;
const int Pairwise = 0;
const int MultipleAlign = 1;

const int OUTSECST = 0;
const int OUTGAP = 1;
const int OUTBOTH = 2;
const int OUTNONE = 3;

const int MAXNAMES = 150;    /* Max chars read for seq. names */ //nige, was 30
//const int MAXNAMESTODISPLAY = 30; // Used for printout. Mark 18-7-07
//const int MAXNAMESTODISPLAY = 10; // Bug #53. Paul 20-12-07
const int MAXNAMESTODISPLAY = 30; //Paul replicate 1.83 behavour 9-2-08
const int MINNAMESTODISPLAY = 10; //Paul replicate 1.83 behavour 9-2-08
const int MAXTITLES = 60;      /* Title length */
const int FILENAMELEN = 256;             /* Max. file name length */

const int  UNKNOWN  = 0;
const int  EMBLSWISS = 1;
const int  PIR      = 2;
const int  PEARSON  = 3;
const int  GDE = 4;
const int  CLUSTAL = 5;    /* DES */
const int  MSF = 6; /* DES */
const int  RSF = 7;    /* JULIE */
const int  USER = 8;    /* DES */
const int  PHYLIP = 9;    /* DES */
const int  NEXUS = 10; /* DES */
const int  FASTA = 11; /* Ramu */

const int  NONE = 0;
const int  SECST = 1;
const int  GMASK = 2;

const int  PROFILE = 0;
const int  SEQUENCE = 1;

const int  BS_NODE_LABELS = 2;
const int  BS_BRANCH_LABELS = 1;

const int  PAGE_LEN = 22;   /* Number of lines of help sent to screen */

const int  PAGEWIDTH = 80;  /* maximum characters on output file page */
const int  LINELENGTH = 60;  /* Output file line length */
const int  GCG_LINELENGTH = 50;

const int NJ = 1;
const int UPGMA = 2;

const int ALIGNMENT = 1;
const int TREE = 2;

const int MinIdentifier = 1;

const string VALID_COMMAND_SEP = "-/";

#ifdef OS_MAC
    const char default_commandsep = '-';
    const char DIRDELIM = '/';
    const int INT_SCALE_FACTOR = 100;  /* Scaling factor to convert float to integer
        for profile scores */

#elif OS_WINDOWS
    const char  default_commandsep = '/';
    const char DIRDELIM = '\\';
    const int INT_SCALE_FACTOR = 100;  /* Scaling factor to convert float to integer
        for profile scores */

#elif OS_UNIX
    const char default_commandsep = '-';
    const char DIRDELIM = '/';
    const int INT_SCALE_FACTOR = 1000; /* Scaling factor to convert float to integer
        for profile scores */
#endif

       
const int NUMRES = 32; /* max size of comparison matrix */
const int INPUT = 0;
const int ALIGNED = 1;

const int LEFT = 1;
const int RIGHT = 2;

const int NODE = 0;
const int LEAF = 1;

const int GAPCOL = 32;        /* position of gap open penalty in profile */
const int LENCOL = 33;        /* position of gap extension penalty in profile */

typedef struct
{
   /* Holds values for the pairwise scales */
   float gapOpenScale; 
   float gapExtendScale;
   int intScale; 
}PairScaleValues;

typedef struct
{
    float scale;
    float intScale;
}PrfScaleValues;
 
typedef struct node
{
     /* phylogenetic tree structure */
    struct node *left;
    struct node *right;
    struct node *parent;
    float dist;
    int leaf;
    int order;
    char name[64];
} stree,  *treeptr;

typedef struct
{
    char title[30];
    char string[30];
} MatMenuEntry;

typedef struct
{
    int noptions;
    MatMenuEntry opt[10];
} MatMenu;

const int MAXMAT = 10;

typedef struct
{
    int llimit;
    int ulimit;
    vector<short>* matptr;
    vector<short>* AAXref;
} SeriesMat;

/*
 * UserMatSeries holds the number of matrices in the series and 
 */
typedef struct
{
    int nmat;
    SeriesMat mat[MAXMAT];
} UserMatrixSeries;

}
#endif

