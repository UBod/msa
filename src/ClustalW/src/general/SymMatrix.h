/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * The SymMatrix class is used to store the distance matrix. It stores it in a double  
 * array.
 * This class throws an out_of_range exception if we try to access an element outside 
 * of the array bounds. It will throw a bad_alloc if we cannot allocate enough memory for 
 * the distance matrix.
 */
#ifndef SYMMATRIX_H
#define SYMMATRIX_H

#include<iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <new>
#include <cstdlib>

namespace clustalw
{
using namespace std;

class SymMatrix
{
    public:
        SymMatrix() : elements(0), numSeqs(0), subElements(0), firstSeq(0), numSeqsInSub(0){;}
        SymMatrix(int size)
         : elements(0),
           numSeqs(0),
           subElements(0),
           firstSeq(0),
           numSeqsInSub(0) 
        {
            // Size is the numSeqs + 1
            numSeqs = size - 1;
            sizeElements = ((numSeqs + 1) * (numSeqs + 2)) >> 1;
            try
            {
                elements = new double[sizeElements];
                setAllArrayToVal(elements, sizeElements, 0.0);
            }
            catch(bad_alloc& e)
            {
                cout << "Could not allocate a distance matrix for " << numSeqs 
                     << " seqs.\n";
                throw e;
            }
            
        }
        ~SymMatrix()
        {
            if(elements)
            {
                delete [] elements;
            }
            if(subElements)
            {
                delete [] subElements;
            }
        }
        
        inline void SetAt(int nRow, int nCol, const double& value)
        {
            index = getIndex(nRow, nCol, numSeqs);
            elements[index] = value; 
        }
        
        inline double GetAt(int nRow, int nCol) 
        {
            index = getIndex(nRow, nCol, numSeqs);
            return elements[index];
        }
   
        inline void ResizeRect(int size, double val = 0.0)
        {
            numSeqs = size - 1;
            sizeElements = ((numSeqs + 1) * (numSeqs + 2)) >> 1;
            if(elements)
            {
                delete [] elements;
            }
            try
            {
                elements = new double[sizeElements];
                setAllArrayToVal(elements, sizeElements, val);
            }
            catch(bad_alloc& e)
            {
                cout << "Could not allocate a distance matrix for " << numSeqs 
                     << " seqs. Need to terminate program.\n";
                throw e;
            }            
        }
        
        inline void setAllArrayToVal(double* array, int size, double val)
        {
            for(int i = 0; i < size; i++)
            {
                array[i] = val;
            }
        } 
        
        inline int getSize()
        {
            return numSeqs + 1;
        }

        inline void clearArray()
        {
            sizeElements = 0;
            numSeqs = 0;
            if(elements)
            {
                delete [] elements;
                elements = 0;
            }
            
            numSeqsInSub = 0;
            sizeSubElements = 0;
            
            if(subElements)
            {
                delete [] subElements;
                subElements = 0;                
            }
        }
         
        void printArray()
        {
            printArray(elements, sizeElements);
        }
        
        void printSubArray()
        {
            printArray(subElements, sizeSubElements);
        }
        
        void printArray(double* array, int size)
        {
            int numThisTime = 1;
            int numSoFar = 0;
            for(int i = 0; i < size; i++)
            {
                
                numSoFar++;
                
                if(numSoFar == numThisTime)
                {
                    cout << " " << setw(4) << array[i] << "\n";
                    numSoFar = 0;
                    numThisTime++;
                }
                else
                {
                    cout << " " << setw(4) << array[i];
                }    
            }
            cout << "\n";
        }
        
        inline int getIndex(const int &i, const int &j, const int &nSeqs) const
        {
            if(i == 0 || j == 0)
            {
                return 0;
            }
            
            int _i = i - 1, _j = j - 1;
            if (_i == _j) 
            {
                if ((_i >= nSeqs) || (_i < 0))
                {
                    throw out_of_range("index out of range\n");
                }
                return (_i * (_i + 3)) >> 1;
            }

            if (_i > _j) 
            {
                if ((_i >= nSeqs) || (_j < 0))
                {
                    throw out_of_range("index out of range\n");
                }
                return ((_i * (_i + 1)) >> 1) + _j;
            }

            if ((_j >= nSeqs) || (_i < 0))
            {
                throw out_of_range("index out of range\n");
            }
            return  ((_j * (_j + 1)) >> 1) + _i;            
        }
        
        //inline
        inline double operator() (unsigned row, unsigned col) const
        {
            return elements[getIndex(row, col, numSeqs)];
        }
        
        inline double& operator() (unsigned row, unsigned col)
        {
            return elements[getIndex(row, col, numSeqs)];
        }
        
        inline double* getElements()
        {
            return elements;
        }
        
        inline double* getDistMatrix(int fSeq, int nSeqsInSub)
        {
            if(fSeq == 1 && numSeqs == nSeqsInSub)
            {
                return elements; // We want the whole array.
            }
            else
            {
                // Calculate the subMatrix and return it!!!
                try
                {
                    if(subElements)
                    {
                        delete [] subElements;
                    }
                    sizeSubElements = ((nSeqsInSub + 1) * (nSeqsInSub + 2)) >> 1;
                    numSeqsInSub = nSeqsInSub;                    
                    
                    subElements = new double[sizeSubElements];
                    setAllArrayToVal(subElements, sizeSubElements, 0.0);
                    
                    int currIndex = 0;
                    subElements[0] = 0.0;
                    int lSeq = fSeq + numSeqsInSub - 1; // NOTE this is wrong!!!!! Need to fix

                    for(int i = fSeq; i <= lSeq; i++)
                    {
                        for(int j = i + 1; j <= lSeq; j++)
                        {
                            currIndex = getIndex(i - fSeq + 1, j - fSeq + 1, numSeqsInSub);
                            subElements[currIndex] = elements[getIndex(i, j, numSeqs)];    
                        }
                    }

                    return subElements;
                }
                catch(bad_alloc& e)
                {
                    cout << "Out of Memory!\n";
                    throw e;
                }                
            }
        }
        
        void clearSubArray()
        {
            if(subElements)
            {
                delete [] subElements;
            }
            subElements = 0; 
            sizeSubElements = 0;
            numSeqsInSub = 0;       
        }
        
        void makeSimilarityMatrix()
        {
            double value;
            for (int i = 0; i < numSeqs; i++)
            {
                SetAt(i + 1, i + 1, 0.0);
                for (int j = 0; j < i; j++)
                {
                    value = 100.0 - (GetAt(i + 1, j + 1)) * 100.0;
                    SetAt(i + 1, j + 1, value);
                    //distMat->SetAt(j + 1, i + 1, value);
                }
            }
        }
        
    private:
        double* elements;
        int sizeElements;
        int numSeqs;
        int index;
        double* subElements; // To be used to return a sub matrix.
        int firstSeq, numSeqsInSub;
        int sizeSubElements;
};

}

#endif
