/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef ARRAY2D_H
#define ARRAY2D_H

#include<iostream>
#include <vector>
namespace clustalw
{

using namespace std;

template <class T>
class Array2D
{
public:
   Array2D():m_dimRow(0), m_dimCol(0){;}
   Array2D(int nRow, int nCol) 
   {
       m_dimRow = nRow;
       m_dimCol = nCol;
       for (int i=0; i < nRow; i++)
       {
           vector<T> x(nCol);
           int y = x.size();
           m_2DVector.push_back(x);
       }
   }
   
   int getRowSize()
   {
       return m_dimRow;
   }
   
   int getColSize()
   {
       return m_dimCol;
   }
   
   void SetAt(int nRow, int nCol, const T& value)
   {
       m_2DVector[nRow][nCol] = value;
   }
   T GetAt(int nRow, int nCol) 
   {
       return m_2DVector[nRow][nCol];
   }
   void GrowRow(int newSize) 
   {
       if (newSize <= (int)m_dimRow)
           return;
       m_dimRow = newSize;
       for(int i = 0 ; i < newSize - (int)m_dimCol; i++)   
       {
           vector<T> x(m_dimRow);
           m_2DVector.push_back(x);
       }
   }
   void GrowCol(int newSize) 
   {
       if(newSize <= (int)m_dimCol)
           return;
       m_dimCol = newSize;
       for (int i=0; i < (int)m_dimRow; i++)
           m_2DVector[i].resize(newSize);
   }
   
   void ResizeRect(int row, int col)
   {
       GrowRow(row);
       GrowCol(col);
   }

   /* This is to get everything initialised to a value */
   void GrowRow(int newSize, const T& value) 
   {
       if (newSize <= m_dimRow)
           return;
       m_dimRow = newSize;
       for(int i = 0 ; i < newSize - m_dimCol; i++)   
       {
           vector<T> x(m_dimRow, value);
           m_2DVector.push_back(x);
       }
   }
   void GrowCol(int newSize, const T& value) 
   {
       if(newSize <= m_dimCol)
           return;
       m_dimCol = newSize;
       for (int i=0; i <m_dimRow; i++)
           m_2DVector[i].resize(newSize, value);
   }
   
   void ResizeRect(int row, int col, const T& value)
   {
       GrowRow(row, value);
       GrowCol(col, value);
   }   
      
   void printArray()
   {
       for(int row=0; row < m_dimRow; row++)
       {
           for(int col=0; col < m_dimCol; col++)
           {
               cout <<" "<< m_2DVector[row][col];  
           }
           cout<<"\n";
       }
   }
   
   void clearArray()
   {
       int size = m_2DVector.size();
       for(int i = 0; i < size; i++)
       {
           m_2DVector[i].clear();
       }
       m_2DVector.clear();
       m_dimRow = 0;
       m_dimCol = 0;
       
   }
   
   vector<T>& operator[](int x)    
   {
       return m_2DVector[x];
   }
private:
   vector< vector <T> > m_2DVector;
   unsigned int m_dimRow;
   unsigned int m_dimCol;
};

}
#endif

