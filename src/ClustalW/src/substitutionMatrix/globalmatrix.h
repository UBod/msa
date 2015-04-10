/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef GLOBALMATRIX_H
#define GLOBALMATRIX_H
#include "SubMatrix.h"
namespace clustalw
{
// Mark Dec 13 2005
/*
 * The reason for having a global sub matrix object is that we have
 * user defined matrices that can be read in at any time and then used later in the 
 * MSA stage.
 */
extern SubMatrix *subMatrix;
}
#endif


