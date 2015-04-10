/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "ClusterTreeOutput.h"

namespace clustalw
{

ClusterTreeOutput::ClusterTreeOutput(SeqInfo* seqInfo, int boot)
: bootstrap(boot)
{
    firstSeq = seqInfo->firstSeq;
    lastSeq = seqInfo->lastSeq;
    numSeqs = seqInfo->numSeqs;
}

void ClusterTreeOutput::printTreeDesc(PhyloTree* phyloTree)
{
    for(int i = 1; i <= numSeqs; i++)
    {
        for(int j = 1; j <= numSeqs; j++)
        {
            cout << " " << phyloTree->treeDesc[i][j];
        }
        cout << "\n";
    }
}


/**
 * The function printPhylipTree is used to print out the unrooted clustered tree in
 * phylip format.
 * @param phyloTree A pointer to the PhyloTree struct that contains the description of the
 *                  tree. 
 * @param tree The file to print the phylip tree to. Must be open!
 * @param alignPtr The alignment object. Needed for the names. 
 * @param distMat The distance matrix that was used for the generation of the cluster. 
 * @param bootTotals Holds the bootstrap values. Only used if the tree has been bootstrapped
 */
void ClusterTreeOutput::printPhylipTree(PhyloTree* phyloTree, ofstream* tree,
      Alignment *alignPtr, DistMatrix* distMat, vector<int>* bootTotals)
{
    int oldRow;
    if (lastSeq - firstSeq + 1 == 2)
    {
        (*tree) << "(" << alignPtr->getName(firstSeq) << ":" << fixed <<setprecision(5) 
                << (*distMat)(firstSeq, firstSeq + 1) << "," << alignPtr->getName(firstSeq + 1) 
                << ":" << fixed << setprecision(5)  << (*distMat)(firstSeq, firstSeq + 1)
                << ");";
        return ;
    }

    (*tree) << "(\n";

    oldRow = twoWaySplit(phyloTree, tree, lastSeq - firstSeq + 1 - 2,
        1, alignPtr, bootTotals);
        
    (*tree) << ":" << fixed << setprecision(5) 
            << phyloTree->leftBranch[lastSeq - firstSeq + 1 - 2];
    
    if ((bootstrap == BS_BRANCH_LABELS) && (oldRow > 0) &&
        ((*bootTotals)[oldRow] > 0))
    {
        (*tree) << "[" << (*bootTotals)[oldRow] << "]";
    }
    (*tree) << ",\n";

    oldRow = twoWaySplit(phyloTree, tree, lastSeq - firstSeq + 1 - 2,
                         2, alignPtr, bootTotals);
    
    (*tree) << ":" << fixed << setprecision(5) 
            << phyloTree->leftBranch[lastSeq - firstSeq + 1 - 1];
    
    if ((bootstrap == BS_BRANCH_LABELS) && (oldRow > 0) &&
        ((*bootTotals)[oldRow] > 0))
    {
        (*tree) << "[" << (*bootTotals)[oldRow] << "]";
    }
    (*tree) << ",\n";

    oldRow = twoWaySplit(phyloTree, tree, lastSeq - firstSeq + 1 - 2,
                         3, alignPtr, bootTotals);
        
    (*tree) << ":" << fixed << setprecision(5) 
            << phyloTree->leftBranch[lastSeq - firstSeq + 1];
    
    if ((bootstrap == BS_BRANCH_LABELS) && (oldRow > 0) &&
        ((*bootTotals)[oldRow] > 0))
    {
        (*tree) << "[" << (*bootTotals)[oldRow] << "]";
    }
    (*tree) << ")";
    
    if (bootstrap == BS_NODE_LABELS)
    {
        (*tree) << "TRICHOTOMY";
    }
    (*tree) << ";\n";
}

int ClusterTreeOutput::twoWaySplit(PhyloTree* phyloTree, ofstream* tree, 
        int startRow, int flag, Alignment *alignPtr, vector<int>* bootTotals)
{
    int row, newRow = 0, oldRow, col, testCol = 0;
    bool singleSeq;

    if (startRow != lastSeq - firstSeq + 1-2)
    {
        (*tree) <<  "(\n";
    }

    for (col = 1; col <= lastSeq - firstSeq + 1; col++)
    {
        if (phyloTree->treeDesc[startRow][col] == flag)
        {
            testCol = col;
            break;
        }
    }

    singleSeq = true;
    for (row = startRow - 1; row >= 1; row--)
    if (phyloTree->treeDesc[row][testCol] == 1)
    {
        singleSeq = false;
        newRow = row;
        break;
    }

    if (singleSeq)
    {
        phyloTree->treeDesc[startRow][testCol] = 0;
        (*tree) << alignPtr->getName(testCol + firstSeq - 1);
        if (startRow == lastSeq - firstSeq + 1 - 2)
        {
            return (0);
        }

        (*tree) << ":" << fixed << setprecision(5) << phyloTree->leftBranch[startRow] 
                << ",\n"; 
    }
    else
    {
        for (col = 1; col <= lastSeq - firstSeq + 1; col++)
        {
            if ((phyloTree->treeDesc[startRow][col] == 1) &&
                (phyloTree->treeDesc[newRow][col] == 1))
            {
                phyloTree->treeDesc[startRow][col] = 0;
            }
        }
        oldRow = twoWaySplit(phyloTree, tree, newRow, 1, alignPtr, bootTotals);
        
        if (startRow == lastSeq - firstSeq + 1-2)
        {
            return (newRow);
        }

        (*tree) << ":" << fixed << setprecision(5) << phyloTree->leftBranch[startRow];
        
        if ((bootstrap == BS_BRANCH_LABELS) && ((*bootTotals)[oldRow] > 0))
        {
            (*tree) << "[" << (*bootTotals)[oldRow] << "]";
        }

        (*tree) << ",\n";
    }


    for (col = 1; col <= lastSeq - firstSeq + 1; col++)
    if (phyloTree->treeDesc[startRow][col] == flag)
    {
        testCol = col;
        break;
    }

    singleSeq = true;
    newRow = 0;
    for (row = startRow - 1; row >= 1; row--)
    if (phyloTree->treeDesc[row][testCol] == 1)
    {
        singleSeq = false;
        newRow = row;
        break;
    }

    if (singleSeq)
    {
        phyloTree->treeDesc[startRow][testCol] = 0;
        (*tree) << alignPtr->getName(testCol + firstSeq - 1);
        (*tree) << ":" << fixed << setprecision(5)  << phyloTree->rightBranch[startRow]
                <<")\n";
    }
    else
    {
        for (col = 1; col <= lastSeq - firstSeq + 1; col++)
        {
            if ((phyloTree->treeDesc[startRow][col] == 1) &&
                (phyloTree->treeDesc[newRow][col] == 1))
            {
                phyloTree->treeDesc[startRow][col] = 0;
            }
        }
        oldRow = twoWaySplit(phyloTree, tree, newRow, 1, alignPtr, bootTotals);
        (*tree) << ":" << fixed << setprecision(5)  << phyloTree->rightBranch[startRow];
        
        if ((bootstrap == BS_BRANCH_LABELS) && ((*bootTotals)[oldRow] > 0))
        {
            (*tree) << "[" << (*bootTotals)[oldRow] << "]";
        }

        (*tree) << ")\n";
    }
    if ((bootstrap == BS_NODE_LABELS) && ((*bootTotals)[startRow] > 0))
    {
        (*tree) << (*bootTotals)[startRow];
    }

    return (startRow);
}

void ClusterTreeOutput::printNexusTree(PhyloTree* phyloTree, ofstream* tree,
                   Alignment *alignPtr, DistMatrix* distMat, vector<int>* bootTotals)
{
    int i;
    int oldRow;

    (*tree) << "#NEXUS\n\n";

    (*tree) << "BEGIN TREES;\n\n";
    (*tree) << "\tTRANSLATE\n";
    
    for(i = 1; i < numSeqs; i++) 
    {
        (*tree) << "\t\t" << i << "\t" << alignPtr->getName(i) <<",\n";
    }
    (*tree) << "\t\t" << numSeqs << "\t" << alignPtr->getName(numSeqs) <<"\n";
    (*tree) << "\t\t;\n";

    (*tree) << "\tUTREE PAUP_1= ";

    if(lastSeq - firstSeq + 1 == 2) 
    {
        (*tree) << "(" << firstSeq << ":" << fixed << setprecision(5) 
        << (*distMat)(firstSeq, firstSeq + 1) << "," << firstSeq + 1 << ":" 
        << fixed << setprecision(5) << (*distMat)(firstSeq, firstSeq + 1) << ")";
    }
    else 
    {

        (*tree) << "(";
 
        oldRow = twoWaySplitNexus(phyloTree, tree,
                       lastSeq - firstSeq + 1 - 2, 1, alignPtr, bootTotals);
        
        (*tree) << ":" << fixed << setprecision(5) 
                << phyloTree->leftBranch[lastSeq-firstSeq+1-2];

        if ((bootstrap == BS_BRANCH_LABELS) && (oldRow > 0) && ((*bootTotals)[oldRow]>0))
        {
            (*tree) << "[" << (*bootTotals)[oldRow] << "]";
        }

        (*tree) << ",";

        oldRow = twoWaySplitNexus(phyloTree, tree, lastSeq - firstSeq + 1 - 2, 2, 
                                  alignPtr, bootTotals);
        (*tree) << ":" << fixed << setprecision(5) 
                       << phyloTree->leftBranch[lastSeq-firstSeq+1-1];

        if ((bootstrap==BS_BRANCH_LABELS) && (oldRow>0) && ((*bootTotals)[oldRow]>0))
        {
            (*tree) << "[" << (*bootTotals)[oldRow] << "]";
        }
        
        (*tree) << ",";

        oldRow = twoWaySplitNexus(phyloTree, tree,
                 lastSeq-firstSeq+1-2, 3, alignPtr, bootTotals);
                 
        (*tree) << ":" << fixed << setprecision(5) 
                << phyloTree->leftBranch[lastSeq-firstSeq+1];
        
        if ((bootstrap==BS_BRANCH_LABELS) && (oldRow>0) && ((*bootTotals)[oldRow]>0))
        {
            (*tree) << "[" << (*bootTotals)[oldRow] << "]";
        }
        
        (*tree) << ")";
        
        if (bootstrap == BS_NODE_LABELS)
        { 
            (*tree) << "TRICHOTOMY";
        }
        (*tree) << ";";
    }
    (*tree) << "\nENDBLOCK;\n";
}

int ClusterTreeOutput::twoWaySplitNexus(PhyloTree* phyloTree, ofstream* tree, 
                int startRow, int flag, Alignment *alignPtr, vector<int>* bootTotals)
{
    int row, newRow = 0, oldRow, col, testCol = 0;
    bool singleSeq;

    if (startRow != lastSeq - firstSeq + 1 - 2)
    {
        (*tree) <<  "(";
    }

    for (col = 1; col <= lastSeq - firstSeq + 1; col++)
    {
        if (phyloTree->treeDesc[startRow][col] == flag)
        {
            testCol = col;
            break;
        }
    }

    singleSeq = true;
    for (row = startRow - 1; row >= 1; row--)
    if (phyloTree->treeDesc[row][testCol] == 1)
    {
        singleSeq = false;
        newRow = row;
        break;
    }

    if (singleSeq)
    {
        phyloTree->treeDesc[startRow][testCol] = 0;
        (*tree) << testCol + firstSeq - 1;;
        if (startRow == lastSeq - firstSeq + 1-2)
        {
            return (0);
        }

        (*tree) << ":" << fixed << setprecision(5) << phyloTree->leftBranch[startRow] << ",";
    }
    else
    {
        for (col = 1; col <= lastSeq - firstSeq + 1; col++)
        {
            if ((phyloTree->treeDesc[startRow][col] == 1) &&
                (phyloTree->treeDesc[newRow][col] == 1))
            {
                phyloTree->treeDesc[startRow][col] = 0;
            }
        }
        oldRow = twoWaySplitNexus(phyloTree, tree, newRow, 1, alignPtr, bootTotals);
        
        if (startRow == lastSeq - firstSeq + 1-2)
        {
            return (newRow);
        }

        (*tree) << ":" << fixed << setprecision(5) << phyloTree->leftBranch[startRow];
        if ((bootstrap == BS_BRANCH_LABELS) && ((*bootTotals)[oldRow] > 0))
        {
            (*tree) << "[" << (*bootTotals)[oldRow] << "]";
        }

        (*tree) << ",";
    }


    for (col = 1; col <= lastSeq - firstSeq + 1; col++)
    if (phyloTree->treeDesc[startRow][col] == flag)
    {
        testCol = col;
        break;
    }

    singleSeq = true;
    newRow = 0;
    for (row = startRow - 1; row >= 1; row--)
    if (phyloTree->treeDesc[row][testCol] == 1)
    {
        singleSeq = false;
        newRow = row;
        break;
    }

    if (singleSeq)
    {
        phyloTree->treeDesc[startRow][testCol] = 0;
        (*tree) << testCol + firstSeq - 1;
        (*tree) << ":" << fixed << setprecision(5) << phyloTree->rightBranch[startRow] << ")";
    }
    else
    {
        for (col = 1; col <= lastSeq - firstSeq + 1; col++)
        {
            if ((phyloTree->treeDesc[startRow][col] == 1) &&
                (phyloTree->treeDesc[newRow][col] == 1))
            {
                phyloTree->treeDesc[startRow][col] = 0;
            }
        }
        oldRow = twoWaySplitNexus(phyloTree, tree, newRow, 1, alignPtr, bootTotals);
        
        (*tree) << ":" << fixed << setprecision(5) << phyloTree->rightBranch[startRow];
        if ((bootstrap == BS_BRANCH_LABELS) && ((*bootTotals)[oldRow] > 0))
        {
            (*tree) << "[" << (*bootTotals)[oldRow] << "]";
        }

        (*tree) << ")";
    }
    if ((bootstrap == BS_NODE_LABELS) && ((*bootTotals)[startRow] > 0))
    {
        (*tree) << (*bootTotals)[startRow];
    }

    return (startRow);
}

void ClusterTreeOutput::printTree(PhyloTree* phyloTree, ofstream* tree,
                                  vector<int>* totals)
{
    int row, col;

    (*tree) << "\n";

    for (row = 1; row <= lastSeq - firstSeq + 1 - 3; row++)
    {
        (*tree) << " \n";
        for (col = 1; col <= lastSeq - firstSeq + 1; col++)
        {
            if (phyloTree->treeDesc[row][col] == 0)
            {
                (*tree) << "*";
            }
            else
            {
                (*tree) << ".";
            }
        }
        if ((*totals)[row] > 0)
        {
            (*tree) << setw(7) << (*totals)[row];
        }
    }
    (*tree) << " \n";
    for (col = 1; col <= lastSeq - firstSeq + 1; col++)
    {
        (*tree) << setw(1) << phyloTree->treeDesc[lastSeq - firstSeq + 1 -2][col];
    }
    (*tree) << "\n";
}

}
