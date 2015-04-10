/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "RootedClusterTree.h"
#include "../../general/OutputFile.h"
#include "UPGMAAlgorithm.h"
#include "RootedTreeOutput.h"
//#include "../RandomGenerator.h"
namespace clustalw
{

auto_ptr<AlignmentSteps> RootedClusterTree::treeFromDistMatrix(RootedGuideTree* phyloTree,DistMatrix* distMat, Alignment *alignPtr, 
                                           int seq1, int nSeqs, string& phylipName)
{
    OutputFile phylipPhyTreeFile;
    auto_ptr<AlignmentSteps> progSteps;
    try
    {
        // Test to see if the inputs are valid
        if(seq1 < 1 || nSeqs < 1)
        {
            cerr << "Invalid inputs into treeFromDistMatrix \n"
                 << "seq1 = " << seq1 << " nSeqs = " << nSeqs << "\n"
                 << "Need to end program!\n";
            throw 1;
            return progSteps;
        }

        float dist;
        string path;
        verbose = false;
        firstSeq = seq1;
        lastSeq = firstSeq + nSeqs - 1;
    
        SeqInfo info;
        info.firstSeq = firstSeq;
        info.lastSeq = lastSeq;
        info.numSeqs = nSeqs;

        utilityObject->getPath(userParameters->getSeqName(), &path);
        
        if(nSeqs >= 2)
        {
            string name = phylipName;
            if(!phylipPhyTreeFile.openFile(&name, 
                             "\nEnter name for new GUIDE TREE           file  ", &path, "dnd",
                             "Guide tree"))
            {
                return progSteps;
            }
            phylipName = name;                    
        }
        else
        {
            return progSteps;
        }
                
        RootedTreeOutput outputTree(&info);
        
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
            progSteps.reset(new AlignmentSteps);
            vector<int> groups;
            groups.resize(nSeqs + 1, 0);
            groups[1] = 1;
            groups[2] = 2;
        }
        else
        {
            UPGMAAlgorithm clusAlgorithm;
            progSteps = clusAlgorithm.generateTree(phyloTree, distMat, &info, false);
            outputTree.printPhylipTree(phyloTree, ptrToFile, alignPtr, distMat);
        }
        return progSteps;
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

void RootedClusterTree::treeFromAlignment(TreeNames* treeNames, Alignment *alignPtr)
{
    try
    {
        OutputFile phylipPhyTreeFile;   
        OutputFile clustalPhyTreeFile;
        OutputFile distancesPhyTreeFile;
        OutputFile nexusPhyTreeFile;
        OutputFile pimFile;
        
        RootedGuideTree phyloTree;
        
        string path;
        int j;
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
        
        // The SeqInfo struct is passed to reduce the number of parameters passed!
        SeqInfo info;
        info.firstSeq = firstSeq;
        info.lastSeq = lastSeq;
        info.numSeqs = numSeqs;
    
        RootedTreeOutput outputTree(&info); // No bootstrap!

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
    
        bootPositions.clear();
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
    

        if (userParameters->getOutputTreeClustal() ||
            userParameters->getOutputTreePhylip() 
            || userParameters->getOutputTreeNexus())
        {
            UPGMAAlgorithm clusAlgorithm;
            clusAlgorithm.setVerbose(true);
            clusAlgorithm.generateTree(&phyloTree, quickDistMat.get(), &info, false,
                                       clustalPhyTreeFile.getPtrToFile());
            clusAlgorithm.setVerbose(false);
        }

        if (userParameters->getOutputTreePhylip())
        {
            outputTree.printPhylipTree(&phyloTree, phylipPhyTreeFile.getPtrToFile(), alignPtr,
                                        quickDistMat.get());
        }

        if (userParameters->getOutputTreeNexus())
        {
            outputTree.printNexusTree(&phyloTree, nexusPhyTreeFile.getPtrToFile(), alignPtr, 
                                       quickDistMat.get());
        }
    
        /** Free up resources!!!!! */
    
        treeGaps.clear();
        bootPositions.clear();
          
    }
    catch(const exception& ex)
    {
        cerr << ex.what() << endl;
        utilityObject->error("Terminating program. Cannot continue\n");
        throw 1;
    }
}
                                
}
