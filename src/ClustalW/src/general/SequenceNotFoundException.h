/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef SEQUENCENOTFOUNDEXCEPTION_H
#define SEQUENCENOTFOUNDEXCEPTION_H

// standard exceptions
#include <iostream>
#include <exception>
using namespace std;

class SequenceNotFoundException: public exception
{
  virtual const char* what() const throw ()
  {
    return "Could not find sequence with given id\n";
  }
};

#endif
