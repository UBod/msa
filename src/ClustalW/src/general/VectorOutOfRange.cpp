/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "VectorOutOfRange.h"
#include <string>
#include <sstream>
#include <iostream>

namespace clustalw
{

VectorOutOfRange::~VectorOutOfRange() throw()
{
    // Dont need to do anything
}

const char* VectorOutOfRange::what() const throw()
{
    std::ostringstream message;
    message << "\nIn Vector "<< _name << ", vector index " << _index << " exceeds bounds 1-" 
            << _max << "\n";
    std::string outputMessage = message.str();
    return outputMessage.c_str();
}

const char* VectorOutOfRange::what()
{
    std::ostringstream message;
    message << "\nIn Vector "<< _name << ", vector index " << _index << " exceeds bounds 1-" 
            << _max << "\n";
    std::string outputMessage = message.str();
    return outputMessage.c_str();
}

}


