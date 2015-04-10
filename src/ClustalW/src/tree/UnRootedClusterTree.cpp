/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "UnRootedClusterTree.h"
#include "../general/utils.h"
#include "RandomGenerator.h"
#include <math.h>
#include <sstream>
#include "../general/OutputFile.h"

namespace clustalw
{
UnRootedClusterTree::UnRootedClusterTree()
{

}
/*
 * The function treeFromAlignment is called to generate a tree from an alignment
 * without a distance matrix. It calculates the distance matrix from the alignment
 * quickly.
 */
 
// Note this function was called phylogenetic_tree before
void UnRootedClusterTree::treeFromAlignment(TreeNames* treeNames, Alignment *alignPtr)
{
    try
    {
        OutputFile phylipPhyTreeFile;   
        OutputFile clustalPhyTreeFile;
        OutputFile distancesPhyTreeFile;
        OutputFile nexusPhyTreeFile;
        OutputFile pimFile;
    
        string path;
        int i, j;
        int overspill = 0;
        int totalDists;
        numSeqs = alignPtr->getNumSeqs(); // NOTE class variable
        
        /**
         * Check if numSeqs is ok
         */
        if(!checkIfConditionsMet(numSeqs, 2))
        {
            return;
        }
              
        firstSeq = 1;
        lastSeq = numSeqs;
        phyloTree = new PhyloTree;
        
        // The SeqInfo struct is passed to reduce the number of parameters passed!
        SeqInfo info;
        info.firstSeq = firstSeq;
        info.lastSeq = lastSeq;
        info.numSeqs = numSeqs;
    
        outputTree = new ClusterTreeOutput(&info, 0); // No bootstrap!
    
        phyloTree->treeDesc.resize(numSeqs + 1, vector<int>(numSeqs + 1));
    
        TreeGroups saveTree(numSeqs + 1, vector<int>(numSeqs + 1));
    
        // NOTE at the moment there is only one type of clustering algorithm, but there will
        // be more!
        clusAlgorithm = new NJTree();

        utilityObject->getPath(userParameters->getSeqName(), &path);

        /**
         * Open the required output files.
         */
        if(!openFilesForTreeFromAlignment(&clustalPhyTreeFile, &phylipPhyTreeFile, 
                        &distancesPhyTreeFile, &nexusPhyTreeFile, &pimFile, treeNames, &path))
        {
            return; // Problem opeing one of the files, cannot continue!
        } 
    
        int _lenFirstSeq = alignPtr->getSeqLength(firstSeq);
    
        bootPositions.resize(_lenFirstSeq + 2);

        for (j = 1; j <= _lenFirstSeq; ++j)
        {
            bootPositions[j] = j;
        }

        /**
         * Calculate quickDist and overspill
         */
        overspill = calcQuickDistMatForAll(clustalPhyTreeFile.getPtrToFile(),
                       phylipPhyTreeFile.getPtrToFile(), nexusPhyTreeFile.getPtrToFile(),
                       pimFile.getPtrToFile(), distancesPhyTreeFile.getPtrToFile(), alignPtr);

        // check if any distances overflowed the distance corrections 
        if (overspill > 0)
        {
            totalDists = (numSeqs *(numSeqs - 1)) / 2;
            overspillMessage(overspill, totalDists);
        }

        if (userParameters->getOutputTreeClustal())
        {
            verbose = true;
        }
        // Turn on file output 
    

        if (userParameters->getOutputTreeClustal() || userParameters->getOutputTreePhylip() 
            || userParameters->getOutputTreeNexus())
        {
            
            clusAlgorithm->setVerbose(true);
            clusAlgorithm->generateTree(phyloTree, quickDistMat.get(), &info,
                                        clustalPhyTreeFile.getPtrToFile());
            clusAlgorithm->setVerbose(false);                            
        }

        for (i = 1; i < numSeqs + 1; i++)
            for (j = 1; j < numSeqs + 1; j++)
            {
                saveTree[i][j] = phyloTree->treeDesc[i][j];
            }

        if (userParameters->getOutputTreePhylip())
        {
            outputTree->printPhylipTree(phyloTree, phylipPhyTreeFile.getPtrToFile(), alignPtr,
                                        quickDistMat.get(), &bootTotals);
        }

        for (i = 1; i < numSeqs + 1; i++)
            for (j = 1; j < numSeqs + 1; j++)
            {
                phyloTree->treeDesc[i][j] = saveTree[i][j];
            }

        if (userParameters->getOutputTreeNexus())
        {
            outputTree->printNexusTree(phyloTree, nexusPhyTreeFile.getPtrToFile(), alignPtr, 
                                       quickDistMat.get(), &bootTotals);
        }
    
        /** Free up resources!!!!! */
    
        treeGaps.clear();
        bootPositions.clear();
          
        delete clusAlgorithm;
        delete phyloTree;
        delete outputTree;
    }
    catch(const exception& ex)
    {
        cerr << ex.what() << endl;
        utilityObject->error("Terminating program. Cannot continue\n");
        throw 1;
    }
}

/*
 * Routine for producing unrooted NJ trees from seperately aligned
 * pairwise distances.  This produces the GUIDE DENDROGRAMS in
 * PHYLIP format.
 */
void UnRootedClusterTree::treeFromDistMatrix(DistMatrix* distMat, 
                                Alignment *alignPtr, int seq1, int nSeqs,
                                 string& phylipName)
{
    OutputFile phylipPhyTreeFile;
    
    try
    {
        // Test to see if the inputs are valid
        if(seq1 < 1 || nSeqs < 1)
        {
            cerr << "Invalid inputs into treeFromDistMatrix \n"
                 << "seq1 = " << seq1 << " nSeqs = " << nSeqs << "\n"
                 << "Need to end program!\n";
            throw 1;
            return;
        }
        PhyloTree phyloTree;
        float dist;
        string path;
        verbose = false;
        firstSeq = seq1;
        lastSeq = firstSeq + nSeqs - 1;
    
        SeqInfo info;
        info.firstSeq = firstSeq;
        info.lastSeq = lastSeq;
        info.numSeqs = nSeqs;
        // Open up the phylipPhyTreeFile now!

        // NOTE that this is not exactly correct. This may cause a problem when outputing
        // a tree for each of the profiles. But then we can pass it in, maybe.
        utilityObject->getPath(userParameters->getSeqName(), &path);
        
        if(nSeqs >= 2)
        {
            string name = phylipName;
            if(!phylipPhyTreeFile.openFile(&name, 
                             "\nEnter name for new GUIDE TREE           file  ", &path, "dnd",
                             "Guide tree"))
            {
                return;
            }
            phylipName = name;                    
        }
        else
        {
            return;
        }
                
        // Not sure about bootstrapping here!
        clusAlgorithm = new NJTree();
        outputTree = new ClusterTreeOutput(&info, 0);
        
        ofstream* ptrToFile = phylipPhyTreeFile.getPtrToFile();
        
        if (nSeqs == 2)
        {
            dist = (*distMat)(firstSeq, firstSeq + 1) / 2.0;
            if(ptrToFile->is_open())
            {
                (*ptrToFile) <<  "(" << alignPtr->getName(firstSeq) << ":" 
                             << setprecision(5)
                             << dist << "," << alignPtr->getName(firstSeq + 1) << ":" 
                             << setprecision(5) << dist <<");\n";
            }
        }
        else
        {
            int dimensions = lastSeq - firstSeq + 2;
            phyloTree.treeDesc.resize(dimensions, vector<int>(dimensions));
/* debugging, also used when -OUTPUTTREE is used */
#if 0
            ofstream debuglog("debug.treeFromDistMatrix.txt", ios::out);
            clusAlgorithm->setVerbose(true);
            clusAlgorithm->generateTree(&phyloTree, distMat, &info, &debuglog);
#else
            clusAlgorithm->generateTree(&phyloTree, distMat, &info);
#endif            
            
            outputTree->printPhylipTree(&phyloTree, ptrToFile, alignPtr,
                                        distMat, &bootTotals);
        }
        delete clusAlgorithm;
        delete outputTree;
    
        //phylipPhyTreeFile.close();
    }
    catch(const exception &ex)
    {
        cerr << "ERROR: Error has occured in treeFromDistMatrix. " 
             << "Need to terminate program.\n"
             << ex.what();
        throw 1;
    }
    catch(...)
    {
        cerr << "ERROR: Error has occured in treeFromDistMatrix. " 
             << "Need to terminate program.\n";      
        throw 1;
    }
}

/*
 * Note before I had this function was accepting a distance matrix. But I think it
 * just uses the quickly generated one. So it doesnt need it anymore. But be careful !!!!
 */
void UnRootedClusterTree::bootstrapTree(TreeNames* treeNames, Alignment *alignPtr)
{
    int i, j;
    int ranno;
    string path;
    
    OutputFile clustalPhyTreeFile;
    ofstream* ptrToClustalFile; 
    OutputFile phylipPhyTreeFile;
    OutputFile nexusPhyTreeFile;
    
    try
    {
        phyloTree = new PhyloTree; 
        PhyloTree sampleTree;
        PhyloTree standardTree;
        PhyloTree saveTree;
        int totalDists, overspill = 0, totalOverspill = 0;
        int nfails = 0;
        numSeqs = alignPtr->getNumSeqs(); 
        firstSeq = 1;
        lastSeq = numSeqs;
        clusAlgorithm = new NJTree(); 

        SeqInfo info;
        info.firstSeq = firstSeq;
        info.lastSeq = lastSeq;
        info.numSeqs = numSeqs;    
    
        /**
         * Check if numSeqs is ok
         */
        if(!checkIfConditionsMet(numSeqs, 4))
        {
            return;
        } 

        if (!userParameters->getOutputTreeClustal() && !userParameters->getOutputTreePhylip()
            && !userParameters->getOutputTreeNexus())
        {
            utilityObject->error("You must select either clustal or phylip or nexus tree output format");
            return;
        }
        
        utilityObject->getPath(userParameters->getSeqName(), &path);

        if(!openFilesForBootstrap(&clustalPhyTreeFile, &phylipPhyTreeFile, &nexusPhyTreeFile,
                                  treeNames, &path))
        {
            return; // There was a problem opening the output files.
        }       
    
        int _lenFirstSeq = alignPtr->getSeqLength(firstSeq);
        bootTotals.clear();
        bootTotals.resize(numSeqs + 1);
        bootPositions.clear();
        bootPositions.resize(_lenFirstSeq + 2);

        for (j = 1; j <= _lenFirstSeq; ++j)
        // First select all positions for
        {
            bootPositions[j] = j;
        }
        
        // the "standard" tree
        overspill = calcQuickDistMatForSubSet(clustalPhyTreeFile.getPtrToFile(),
                phylipPhyTreeFile.getPtrToFile(), nexusPhyTreeFile.getPtrToFile(), alignPtr);

        // check if any distances overflowed the distance corrections
        if (overspill > 0)
        {
            totalDists = (numSeqs *(numSeqs - 1)) / 2;
            overspillMessage(overspill, totalDists);
        }

        treeGaps.clear();
    
        if (userParameters->getOutputTreeClustal())
        {
            verbose = true;
        }
        // Turn on screen output
        phyloTree->treeDesc.resize(numSeqs + 1, vector<int>(numSeqs + 1));
    
        // compute the standard tree 

        if (userParameters->getOutputTreeClustal() || userParameters->getOutputTreePhylip() ||
            userParameters->getOutputTreeNexus())
        {
            clusAlgorithm->setVerbose(true); 
            clusAlgorithm->generateTree(phyloTree, quickDistMat.get(), &info,
                                        clustalPhyTreeFile.getPtrToFile());
        }

        ptrToClustalFile = clustalPhyTreeFile.getPtrToFile();

        promptForBoolSeedAndNumTrials();

        RandomGenerator randGenerator(userParameters->getBootRanSeed());
        
        /**
         * Print bootstrap information to top of clustal bootstrap file!
         */
        if (userParameters->getOutputTreeClustal())
        {
            printBootstrapHeaderToClustalFile(ptrToClustalFile);
        }

        verbose = false; // Turn OFF screen output
        clusAlgorithm->setVerbose(false);
        
        sampleTree.treeDesc.resize(numSeqs + 1, vector<int>(numSeqs + 1));
    
        if (userParameters->getMenuFlag())
        {
            cout <<  "\n\nEach dot represents 10 trials\n\n";
        }
    
        totalOverspill = 0;
        nfails = 0;
        int lenSeq1 = alignPtr->getSeqLength(1);
        for (i = 1; i <= userParameters->getBootNumTrials(); ++i)
        {
            for (j = 1; j <= alignPtr->getSeqLength(firstSeq); ++j)
            {
                // select alignment positions for            
                ranno = randGenerator.addRand((unsigned long)lenSeq1) + 1;
                bootPositions[j] = ranno; // bootstrap sample 
            }
            
            overspill = calcQuickDistMatForSubSet(clustalPhyTreeFile.getPtrToFile(),
                phylipPhyTreeFile.getPtrToFile(), nexusPhyTreeFile.getPtrToFile(), alignPtr,
                true);
            
            if (overspill > 0)
            {
                totalOverspill = totalOverspill + overspill;
                nfails++;
            }

            treeGaps.clear();
        
            if (userParameters->getOutputTreeClustal() ||
                userParameters->getOutputTreePhylip() || userParameters->getOutputTreeNexus())
            {
                // NOTE this is bad to pass in clustalPhyTreeFile
                clusAlgorithm->generateTree(&sampleTree, quickDistMat.get(), &info,
                                            clustalPhyTreeFile.getPtrToFile());
            }

            sampleTree.leftBranch.clear();
            sampleTree.rightBranch.clear();

            compareTree(phyloTree, &sampleTree, &bootTotals, lastSeq - firstSeq + 1);
        
            if (userParameters->getMenuFlag())
            {
                if (i % 10 == 0)
                {
                    cout <<  ".";
                }
                if (i % 100 == 0)
                {
                    cout << "\n";
                }
            }
        }

        // check if any distances overflowed the distance corrections 
        if (nfails > 0)
        {
            totalDists = (numSeqs *(numSeqs - 1)) / 2;
            printErrorMessageForBootstrap(totalOverspill, totalDists, nfails);
        }

        bootPositions.clear();

        outputTree = new ClusterTreeOutput(&info, userParameters->getBootstrapFormat());

        /**
         * Print ClustalTree with bootTotals
         */
        if (userParameters->getOutputTreeClustal())
        {
            outputTree->printTree(phyloTree, clustalPhyTreeFile.getPtrToFile(), &bootTotals);
        }
    
        /**
         * Print phylip tree with boottotals.
         */
        if (userParameters->getOutputTreePhylip())
        {
            // Save the old tree!
            saveTree.treeDesc.resize(numSeqs + 1, vector<int>(numSeqs + 1));
            saveTree.treeDesc.assign(phyloTree->treeDesc.begin(), phyloTree->treeDesc.end());
        
            outputTree->printPhylipTree(phyloTree, phylipPhyTreeFile.getPtrToFile(),
                                        alignPtr, quickDistMat.get(), &bootTotals);
            // reassign the old values!                            
            phyloTree->treeDesc.assign(saveTree.treeDesc.begin(), saveTree.treeDesc.end());
        }

        /**
         * print nexus tree with boottotals
         */
        if (userParameters->getOutputTreeNexus())
        {
            outputTree->printNexusTree(phyloTree, nexusPhyTreeFile.getPtrToFile(), alignPtr, 
                                       quickDistMat.get(), &bootTotals);
        }
    
        delete phyloTree;
        delete clusAlgorithm;
        delete outputTree;
    }
    catch(const exception &ex)
    {

    }
}

}
