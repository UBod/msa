/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef VECTORUTILITY_H
#define VECTORUTILITY_H
#include <vector>

namespace vectorutils
{

/**
 * The function mergeVectors will add the contents of 
 * the other two vectors to the end of the first vector.
 */
template<class T>
void mergeVectors(std::vector<T>* vecToAddTo, std::vector<T>* vector1, std::vector<T>* vector2)
{
    typename std::vector<T>::iterator beginInsert = vecToAddTo->end();
    typename std::vector<T>::iterator beginVec1 = vector1->begin();
    typename std::vector<T>::iterator endVec1 = vector1->end();
    typename std::vector<T>::iterator beginVec2 = vector2->begin();
    typename std::vector<T>::iterator endVec2 = vector2->end();
    
    // Add the first vector
    vecToAddTo->insert(beginInsert, beginVec1, endVec1);
    
    beginInsert = vecToAddTo->end();
    // Add the second vector
    vecToAddTo->insert(beginInsert, beginVec2, endVec2);
}

/**
 * The function mergeVectors will add the contents of the second vector 
 * to the end of the first vector.
 */
template<class T>
void mergeVectors(std::vector<T>* vecToAddTo, std::vector<T>* vector1)
{
    typename std::vector<T>::iterator beginInsert = vecToAddTo->end();
    typename std::vector<T>::iterator beginVec1 = vector1->begin();
    typename std::vector<T>::iterator endVec1 = vector1->end();
    
    // Add the first vector
    vecToAddTo->insert(beginInsert, beginVec1, endVec1);
}

}

#endif
