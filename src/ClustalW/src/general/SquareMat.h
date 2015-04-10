
/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#include<iostream>
#include <vector>
#include <iomanip>
namespace clustalw
{
using namespace std;

template <class T>
class SquareMat
{
public:
   SquareMat():m_dimRow(0), m_dimCol(0){;}
   SquareMat(int size) 
   {
       m_dimRow = size;
       m_dimCol = size;
       for (int i=0; i < size; i++)
       {
           vector<T> x(size);
           int y = x.size();
           m_2DVector.push_back(x);
       }
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
       if (newSize <= m_dimRow)
           return;
       m_dimRow = newSize;
       for(int i = 0 ; i < newSize - m_dimCol; i++)   
       {
           vector<T> x(m_dimRow);
           m_2DVector.push_back(x);
       }
   }
   void GrowCol(int newSize) 
   {
       if(newSize <= m_dimCol)
           return;
       m_dimCol = newSize;
       for (int i=0; i <m_dimRow; i++)
           m_2DVector[i].resize(newSize);
   }
   
   void ResizeRect(int size)
   {
       GrowRow(size);
       GrowCol(size);
   }

   int getSize()
   {
       return m_dimRow;
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
         
   void printArray()
   {
       for(int row=0; row < m_dimRow; row++)
       {
           for(int col=0; col < m_dimCol; col++)
           {
               cout << setprecision(20) << " "<< m_2DVector[row][col];  
           }
           cout<<"\n";
       }
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

