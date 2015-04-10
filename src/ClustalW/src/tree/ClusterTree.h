/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * This is the interface class for the clustering part of the code. The 3 public functions
 * are the 3 things that can be done. 1) Generate a tree from alignment (treeFromAlignment), 
 * 2) Generate a tree from a distance matrix (treeFromDistMatrix), or
 * 3) Bootstrap a tree (bootstrapTree).
 **/
#ifndef CLUSTALTREEBASE_H
#define CLUSTALTREEBASE_H

#include <fstream>
#include <memory>
#include <iostream>
#include <vector>
#include <exception>
#include "../alignment/Alignment.h"
#include "NJTree.h"
#include "ClusterTreeOutput.h"
#include "../general/OutputFile.h"
#include "ClusterTreeAlgorithm.h"
#include <string>

namespace clustalw
{
class OutputFile;

class ClusterTree
{
    public:
        /* Functions */
        ClusterTree(); 
        /* Attributes */

    protected: // This is because we will have a derived class that needs these. 
        /* Functions */
        void overspillMessage(int overspill,int total_dists);
        void treeGapDelete(clustalw::Alignment *alignPtr);   
        int dnaDistanceMatrix(ofstream* treeFile, clustalw::Alignment *alignPtr);
        int protDistanceMatrix(ofstream* treeFile, clustalw::Alignment *alignPtr);
        bool isAmbiguity(int c);
        void calcPercIdentity(ofstream* pfile, clustalw::Alignment *alignPtr);
        void compareTree(clustalw::PhyloTree* tree1, clustalw::PhyloTree* tree2, vector<int>* hits, int n);
        //string getOutputFileName(const string prompt, string path, 
        //                              const string fileExtension);
        bool transition(int base1, int base2);
        void distanceMatrixOutput(ofstream* outFile, clustalw::DistMatrix* matToPrint,
                                  clustalw::Alignment *alignPtr);
        bool openFilesForBootstrap(clustalw::OutputFile* clustalFile, clustalw::OutputFile* phylipFile,
                         clustalw::OutputFile* nexusFile, clustalw::TreeNames* treeNames, string* path);
        bool openFilesForTreeFromAlignment(clustalw::OutputFile* clustalFile, clustalw::OutputFile* phylipFile,
                            clustalw::OutputFile* distFile, clustalw::OutputFile* nexusFile, clustalw::OutputFile* pimFile, 
                            clustalw::TreeNames* treeNames, string* path);
        int calcQuickDistMatForAll(ofstream* clustalFile, ofstream* phylipFile, 
                                   ofstream* nexusFile, ofstream* pimFile, 
                                   ofstream* distFile, clustalw::Alignment* alignPtr);
                                   
        int calcQuickDistMatForSubSet(ofstream* clustalFile, ofstream* phylipFile, 
                                      ofstream* nexusFile, clustalw::Alignment* alignPtr, 
                                      bool inBootLoop = false);
        void printBootstrapHeaderToClustalFile(ofstream* clustalFile);
        void promptForBoolSeedAndNumTrials();
        void printErrorMessageForBootstrap(int totalOverspill, int totalDists, int nfails);
        bool checkIfConditionsMet(int numSeqs, int min);
        /* Attributes */
        ClusterTreeAlgorithm* clusAlgorithm;
        auto_ptr<clustalw::DistMatrix> quickDistMat;
        
        vector<int> bootTotals; 
        vector<int> bootPositions;
        bool verbose;
        vector<int> treeGaps; 
        int numSeqs;
        int firstSeq;
        int lastSeq;
        string bootstrapPrompt;
        string bootstrapFileTypeMsg;
};

}
#endif
