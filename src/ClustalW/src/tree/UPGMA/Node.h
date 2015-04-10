#ifndef _NODE_H
#define _NODE_H
#include <limits>
#include <vector>
#include <string>
#include "upgmadata.h"
namespace clustalw
{

using namespace std;

class Node 
{
    public:
    
        Node(int seqNum, double *aptrToDistMatRow, int numDists);
        double getDist(int index){return ptrToDistMatRow[index];}
        void setDist(int index, double dist){ptrToDistMatRow[index] = dist;}
        double* getPtrToDistMatRow(){return ptrToDistMatRow;}
        void setDistMatRowToNull(){ptrToDistMatRow = 0;}
        int getNumDists(){return numDists;}
        double getMinDist(){return minDist;}
        int getIndexToMinDist(){return indexToMinDist;}
        void setMinDist(int index, double d){indexToMinDist = index; minDist = d;}
        int getOrder(){return order;}
        void setOrder(int o){order = o;}
        int getSeqNum(){return seqNum;}
        void setSeqNum(int sNum){seqNum = sNum;}
        void printNodeInfo();
        void printDistMatRow();
        void printMinDist();
        string elementsToString();        
        Node* getLeft(){return left;}
        Node* getRight(){return right;}
        void setLeft(Node* l){left = l;}
        void setRight(Node* r){right = r;}
        double getHeight(){return height;}
        void setHeight(double h){height = h;}
        vector<int>* getPtrToElements(){return &allElements;}
        int getFirstElement(){return allElements[0];}
        void clearElements(){allElements.clear();}
        int getFirstElem(){return allElements[0];}
        bool isLeafNode()
        {    
            if(left == 0 && right == 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        
        void merge(Node **rightNode, double _height);
        void findMinDist();
        void printElements();      
        void makeEmpty();
        void makeEmpty(Node* t);
        /* Attributes */ 
        Node *next; // For linked list of nodes
        Node *left, *right;
        int size;
        int seqNum;
        double height;
        vector<int> allElements; // THis is for the groups!!!   
        double *ptrToDistMatRow; 
        double minDist; // minimal distance
        int indexToMinDist; // index of minimal distance
        int numDists;
        int order;
};

}

#endif
