#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "Node.h"
#include "../../general/VectorUtility.h"
#include "../../general/debuglogObject.h"
#include "../../general/clustalw.h"
#include <iostream>
#include <sstream>
namespace clustalw

{

Node::Node(int _seqNum, double *dists, int numDist)
    : next(0),
      left(0),
      right(0),
      size(1),
      seqNum(_seqNum),
      height(0.0),
      ptrToDistMatRow(dists),
      minDist(numeric_limits<double>::max()),
      indexToMinDist(-1),
      numDists(numDist),
      order(0)
{
    allElements.resize(1);
    allElements[0] = seqNum;
                    
    if (ptrToDistMatRow)
    {
        findMinDist();
    }
}


void Node::merge(Node **rightNode, double _height)
{
    left = new Node(*this);
    right = *rightNode;

    left->ptrToDistMatRow = 0;
    size = left->size + right->size;
    seqNum = -1;
    height = _height;
    left->height = height;
    right->height = height;
    
    vectorutils::mergeVectors(&allElements, &(right->allElements));
    right->allElements.clear();      
    
    if (next == right)
    {
        next = right->next;
    }
    else
    {    
        *rightNode = right->next;
    }
}


void Node::findMinDist()
{    
    double *distIterator = ptrToDistMatRow;
    double *minDistSoFar = distIterator++;
    
    // We search from the end of our area of the array       
    for(int i = numDists; --i; distIterator++) // When --i gets to zero it will stop
    {    
        if ((*distIterator >= 0) && (*distIterator < *minDistSoFar))
        {
            minDistSoFar = distIterator;
        }
    }
    
    minDist = *minDistSoFar;
    indexToMinDist = minDistSoFar - ptrToDistMatRow;
}
 
void Node::printElements()
{
    for(int i = 0; i < (int)allElements.size(); i++)
    {
        cout << " " << allElements[i];
    }
    cout << "\n";               
}

string Node::elementsToString()
{
    ostringstream elems;
    for(int i = 0; i < (int)allElements.size(); i++)
    {
        elems << " " << allElements[i];
    }
    return elems.str();
}

void Node::makeEmpty()
{
    makeEmpty(this);
}

void Node::makeEmpty(Node* t)
{
    if(t != 0)
    {
        makeEmpty(t->left);
        makeEmpty(t->right);
        delete t;
    }
    t = 0;
}
                                           
}
