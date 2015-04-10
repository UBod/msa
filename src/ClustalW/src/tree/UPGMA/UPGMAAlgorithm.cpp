#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "UPGMAAlgorithm.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>

#include "../../general/SymMatrix.h"
#include "../../general/debuglogObject.h"
#include "../../general/clustalw.h"
#include "upgmadata.h"

namespace clustalw
{

UPGMAAlgorithm::UPGMAAlgorithm()
: overwriteMatrix(false),
  numSeqs(0),
  verbose(false),
  orderNode1(0),
  orderNode2(0),
  orderNewNode(0)
{}

auto_ptr<AlignmentSteps> UPGMAAlgorithm::generateTree(RootedGuideTree* phyTree,
                                          DistMatrix* distMat,
                                          SeqInfo* seqInfo, bool overwrite,
                                          ofstream* tree)
{
    if (tree == 0 || !tree->is_open())
    {
        verbose = false;
    }  

    if (verbose)
    {
        (*tree) << "\n\n\t\t\tUPGMA Method\n"
                << "\n\n This is a ROOTED tree\n"
                << "\n Numbers in parentheses are branch lengths\n\n";
    }
          
    progSteps.reset(new AlignmentSteps);
    
    Node** clusters; 
    Node* root;
    numSeqs = seqInfo->numSeqs;
    const int sizeDistMat = ((numSeqs + 1) * (numSeqs + 2)) / 2;
        
    double* elements = overwrite ? 
                      distMat->getDistMatrix(seqInfo->firstSeq, seqInfo->numSeqs) : 
                      (double *)memcpy(new double[sizeDistMat], 
                      distMat->getDistMatrix(seqInfo->firstSeq, seqInfo->numSeqs),
                      sizeDistMat * sizeof(double));
       
    clusters = initialiseNodes(elements, seqInfo->firstSeq);
    root = doUPGMA(clusters, tree);
    
    phyTree->setRoot(root);
    delete [] clusters;

    if(!overwrite)
    {
        delete [] elements;
    }
    distMat->clearSubArray(); 
    
    return progSteps;  
}


Node** UPGMAAlgorithm::initialiseNodes(double* distanceMatrix, int fSeq)
{    
    int firstSeq = fSeq; 

    
    Node** nodes = new Node*[numSeqs];
    Node** nodeIter = nodes;

    *nodes = new Node(firstSeq, 0, 0);
    
    distanceMatrix++;
    
    // Move to first position in distanceMatrix.
    for(int elementIndex = 1, e = numSeqs; elementIndex < e; 
        distanceMatrix += ++elementIndex) 
    {
        Node* newcluster = new Node(elementIndex + firstSeq, 
                                    distanceMatrix, elementIndex);             
        (*nodeIter++)->next = newcluster;
        *nodeIter = newcluster;
    }

  return nodes;
}

void UPGMAAlgorithm::printAllNodes(Node** nodes)
{
    int numNodes = 0;
    for(Node* nodeIter = *nodes; nodeIter; nodeIter = nodeIter->next)
    {
        numNodes++;
        cout << "Node " << numNodes << "\n";
        nodeIter->printElements();
        cout << "\n\n";
    }
    cout << "There are " << numNodes << " nodes\n";  
}

Node* UPGMAAlgorithm::doUPGMA(Node** nodes, ofstream* tree)
{
    if (tree == 0 || !tree->is_open())
    {
        verbose = false;
    } 
        
    string type1, type2;
    int step = 0;
    
    while((*nodes)->next) // While there is more than 1 node.
    {        
        step++;
        if (verbose)
        {
            (*tree) <<  "\n Cycle" << setw(4) << step << "     = ";
        }
        Node** ptrNodeWithMin;
        
        /** 
         * STEP 1 find node with min distance
         */
        ptrNodeWithMin = getNodeWithMinDist(nodes); 
        Node* nodeToJoin2 = *ptrNodeWithMin; 
        const int indexToMinDist = nodeToJoin2->indexToMinDist; 
        Node* const nodeToJoin1 = nodes[indexToMinDist];

        orderNode1 = nodeToJoin1->size;
        orderNode2 = nodeToJoin2->size;
        orderNewNode = orderNode1 + orderNode2;

        /**
         * STEP 2 Recompute the dist Matrix row for nodeToJoin1
         * example: Join nodes 2 and 6
         * 
         * 1  0
         * 2  X 0
         * 3  0 0 0
         * 4  0 0 0 0
         * 5  0 0 0 0 0
         * 6  0 0 0 0 0 0
         * 7  0 0 0 0 0 0 0
         */
        double* nodeToJoin2DistIter = nodeToJoin2->ptrToDistMatRow;
                
        if (indexToMinDist != 0)
        {
            recomputeNodeToJoin1DistMatRow(nodeToJoin1, &nodeToJoin2DistIter);
        }
        
        /**
         * STEP 3 Recompute all distances in column 
         * example: Join nodes 2 and 6
         * 
         * 1  0
         * 2  C 0
         * 3  0 X 0
         * 4  0 X 0 0
         * 5  0 X 0 0 0
         * 6  0 0 0 0 0 0
         * 7  0 X 0 0 0 0 0
         */
                 
        computeAllOtherDistsToNewNode(nodeToJoin1, nodeToJoin2, &nodeToJoin2DistIter);

        /**
         * STEP 4 Add the step to the progSteps.
         * This creates the multiple alignment steps.
         */     
        addAlignmentStep(&nodeToJoin1->allElements, &nodeToJoin2->allElements);

        
        double minDistance = (*ptrNodeWithMin)->minDist;
        
        double height = 0.0;                                  
        height = minDistance / 2.0;
        
        if(verbose)
        {
            if(nodeToJoin1->allElements.size() > 1)
            {
                type1 = "NODE: ";
            }
            else
            {
                type1 = "SEQ: ";
            }
            if(nodeToJoin2->allElements.size() > 1)
            {
                type2 = "NODE: ";
            }
            else
            {
                type2 = "SEQ: ";
            }
            (*tree) << type1 << nodeToJoin1->allElements[0] << " (" << setw(9) 
                    << setprecision(5) << height << ") joins " << type2 
                    << setw(4) << nodeToJoin2->allElements[0] << " (" 
                    << setw(9) << setprecision(5) << height << ")";
        }
        /**
         * STEP 5 merge 2 nodes
         */              
        nodeToJoin1->merge(ptrNodeWithMin, height);
  
    }

    return *nodes; 
}

/**
 * This function returns a Node object that has the minimum distance of
 * all the nodes left.
 */
Node** UPGMAAlgorithm::getNodeWithMinDist(Node** nodes)
{
    Node** ptrMinNode = NULL;
    double minDistance = numeric_limits<double>::max();
    //Start at node 1, check see if it points to something, then move on the next one.
    Node** nodeIter;
    for(nodeIter = &((*nodes)->next); *nodeIter; 
        nodeIter = &(*nodeIter)->next)
    {
        if ((*nodeIter)->getMinDist() < minDistance) 
        {
            minDistance = (*nodeIter)->getMinDist();
            ptrMinNode = nodeIter; // ptrMinNode will hold a ptr to node with sm'st dist
        }
    }
    return ptrMinNode;    
}

/**
 * This function is used to recompute the nodeToJoin1 dist mat row. Each row in the distance
 * matrix has the distances to all nodes before it, not after. For example the dist mat row
 * for Node 3 would have the distances to nodes 1 and 2.
         * 1  0
         * 2  0 0
         * 3  X X 0      Dists to node 1 and 2.
         * 4  0 0 0 0
         * 5  0 0 0 0 0
         * 6  0 0 0 0 0 0
         * 7  0 0 0 0 0 0 0 
 */
void UPGMAAlgorithm::recomputeNodeToJoin1DistMatRow(Node* nodeToJoin1, 
                                                    double** nodeToJoin2DistIter)
{
    double* nodeToJoin1DistIter = nodeToJoin1->ptrToDistMatRow;
    // calculate new distance
    *nodeToJoin1DistIter = calcNewDist(*nodeToJoin1DistIter, **nodeToJoin2DistIter);
                             
    const double* minIndex1 = nodeToJoin1DistIter;
    nodeToJoin1DistIter++; 
    (*nodeToJoin2DistIter)++;
    int numDistToUpdate = nodeToJoin1->numDists - 1;
            
    // For each of the distances in nodeToJoin1           
    while(numDistToUpdate > 0)
    {
        if (*nodeToJoin1DistIter >= 0) 
        {
            // Calculate the average
            *nodeToJoin1DistIter = calcNewDist(*nodeToJoin1DistIter, **nodeToJoin2DistIter);
            
            if (*nodeToJoin1DistIter < *minIndex1)
            {
                minIndex1 = nodeToJoin1DistIter;
            }
        }
        nodeToJoin1DistIter++;
        (*nodeToJoin2DistIter)++;
        numDistToUpdate--;
    }

    // We have found the minimum distance            
    nodeToJoin1->minDist = *minIndex1;
    nodeToJoin1->indexToMinDist = minIndex1 - nodeToJoin1->ptrToDistMatRow;
}

/**
 * This function is used to recompute all the other distances. It does this by calling
 * two other functions. We first compute all the distances until we get to node 2. 
 * Then we call another function to do the rest. This is because nodes after node 2 may
 * have had EITHER node1 of node2 as their min distance.
 */
void UPGMAAlgorithm::computeAllOtherDistsToNewNode(Node* nodeToJoin1, Node* nodeToJoin2,
                                                   double** nodeToJoin2DistIter)
{
    computeDistsUpToNodeToJoin2(nodeToJoin1, nodeToJoin2, nodeToJoin2DistIter);
    computeDistsForNodesAfterNode2(nodeToJoin2);
}

/**
 * STEP 3A:
 * This function recomputes the distances in column until we get to nodeToJoin2 
 * example: Join nodes 2 and 6
 * 
 * 1  0
 * 2  C 0
 * 3  0 X 0
 * 4  0 X 0 0
 * 5  0 X 0 0 0
 * 6  0 0 0 0 0 0
 * 7  0 0 0 0 0 0 0
 *
 * For each node until we get to nodeToJoin2    
 *    If newdistance is less than the old min distance, set it to this.
 *    else if its greater than old minDist and indexTominDist is the same, recompute min
 *    else leave the min the same as it hasnt changed.
 */
void UPGMAAlgorithm::computeDistsUpToNodeToJoin2(Node* nodeToJoin1, Node* nodeToJoin2, double** nodeToJoin2DistIter)
{
    const int indexToMinDist = nodeToJoin2->indexToMinDist;
    
    movePtrPastUnusedDistances(nodeToJoin2DistIter);
    
    Node* nodeIter;
    // For each node until we get to the second node we are joining           
    for(nodeIter = nodeToJoin1->next; nodeIter != nodeToJoin2; nodeIter = nodeIter->next) 
    {            
        (*nodeToJoin2DistIter)++; // Skip the distance to the node we are joining with
        movePtrPastUnusedDistances(nodeToJoin2DistIter);
                       
        double distToNode = nodeIter->ptrToDistMatRow[indexToMinDist];
            
        double newDistToNode = calcNewDist(distToNode, **nodeToJoin2DistIter);
        nodeIter->ptrToDistMatRow[indexToMinDist] = newDistToNode;
            
                        
        if (newDistToNode < nodeIter->minDist)
        { /** new value is smaller than current min. */
            nodeIter->minDist = newDistToNode;
            nodeIter->indexToMinDist = indexToMinDist;
        }
        else if ((newDistToNode > nodeIter->minDist) && 
                 (nodeIter->indexToMinDist == indexToMinDist)) 
        { /** indexToMinDist was the min dist, but it is now a larger num, recompute */
            nodeIter->findMinDist();
        }
            
    }

}

/**
 * STEP 3B Recompute distance for nodes after nodeToJoin2
 * example: Join nodes 2 and 6
 * 
 * 1  0
 * 2  C 0
 * 3  0 C 0
 * 4  0 C 0 0
 * 5  0 C 0 0 0
 * 6  0 0 0 0 0 0
 * 7  0 X 0 0 0 0 0
 *
 * For each node until we get to the end.            
 *     if dist is less than minDist, set mindist to new distance, set index to index
 *     else if dist is greater than mindist and the index was either nodetojoin1 or
 *     nodetojoin2, recompute distance, set entry for nodetojoin2 to NOTUSED
 *     else set the distance to unused. 
 */        

void UPGMAAlgorithm::computeDistsForNodesAfterNode2(Node* nodeToJoin2)
{
    int indexToNode2 = nodeToJoin2->numDists;
    const int indexToMinDist = nodeToJoin2->indexToMinDist;
    
    Node* nodeIter;       

    for(nodeIter = nodeToJoin2->next; nodeIter; nodeIter = nodeIter->next) 
    {              
        double &distUpdate = nodeIter->ptrToDistMatRow[indexToMinDist];
        
        distUpdate = 
	  ((distUpdate * orderNode1) + 
	   (nodeIter->ptrToDistMatRow[indexToNode2] * orderNode2)) 
	  / orderNewNode;

	/* The comparison (distUpdate > nodeIter->minDist) is unsafe. 
	 * Specifically, we get different results for 32/64bit machines, 
	 * which leads to different branching in the if/else statement 
	 * and nasty behaviour down the line. 
	 * Using all of Pfam as a benchmark distUpdate can 'wobble' by 40 ULPs
         * (Unit of Least Precision) which is difficult to translate into 
	 * a maximum relative error, so we pick COMPARE_REL_EPSILON 
	 * phenomenologically: approx 2E-15 was the biggest we saw in Pfam, 
	 * but suggest 1E-14 (for good measure). 
	 * Using ULPs, eg, *(long long int*)&(distUpdate),
	 * would be better (more elegant?) but gives problems 
	 * with aliasing and/or performance reduction. 
	 * FS, 2009-04-30 
	 */
#define COMPARE_REL_EPSILON 1E-14
        if ( (distUpdate < nodeIter->minDist) &&
	     ((nodeIter->minDist-distUpdate) /
	      nodeIter->minDist > COMPARE_REL_EPSILON) )
	  { /** new distance is smaller */
            nodeIter->minDist = distUpdate;
            nodeIter->indexToMinDist = indexToMinDist;
	  }
        else if (((distUpdate > nodeIter->minDist) &&
		  ((distUpdate-nodeIter->minDist) / 
		   distUpdate > COMPARE_REL_EPSILON) &&
		  (nodeIter->indexToMinDist == indexToMinDist)) ||
		 (nodeIter->indexToMinDist == indexToNode2))
	  { /** if old min dist was to either nodeToJoin1 or nodeToJoin2 */
            nodeIter->ptrToDistMatRow[indexToNode2] = UPGMA_BLANKDIST;
            nodeIter->findMinDist();
	  }
        else
	  {
            nodeIter->ptrToDistMatRow[indexToNode2] = UPGMA_BLANKDIST;
	  }
    }        

}

void UPGMAAlgorithm::addAlignmentStep(vector<int>* group1, vector<int>* group2)
{
    int sizeGroup1 = group1->size();
    int sizeGroup2 = group2->size();
    
    vector<int> groups;
    groups.resize(numSeqs + 1, 0);
    int sizeGroup = groups.size();
    
    for(int i = 0; i < sizeGroup1 && (*group1)[i] < sizeGroup; i++)
    {
        groups[(*group1)[i]] = 1; // Set to be apart of group1
    }
    
    for(int i = 0; i < sizeGroup2 && (*group2)[i] < sizeGroup; i++)
    {
        groups[(*group2)[i]] = 2; // Set to be apart of group1
    }
    
    progSteps->saveSet(&groups); 
}

/**
 * At the moment we are only using average distance, UPGMA, but we can also use this 
 * function to have different distance measures, min, max etc.
 */
double UPGMAAlgorithm::calcNewDist(double dist1, double dist2)
{
    return ((dist1 * orderNode1) + (dist2 * orderNode2)) / orderNewNode;
}
  
}
