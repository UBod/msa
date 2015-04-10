/**
 * Author: Nigel Brown
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * InFileStream subclasses std::ifstream, adding a check for the end-of-line
 * character convention in the input file. This is then used by the getline()
 * member as the line delimiter, unless the caller supplies an explicit
 * delimiter.
 *
 * Created: 09-02-07,Nigel Brown(EMBL)
 ***************************************************************************/
#ifndef INFILESTREAM_H
#define INFILESTREAM_H

#include <string>
#include <fstream>
#include <iostream>
#include <memory>

class InFileStream : public std::ifstream
{
  public:
    InFileStream();
    InFileStream(const char *filename);
    //- InFileStream(const InFileStream &copy);

    void open(const char *filename);
    void close();

    //int get();
    //bool is_open();
    std::istream& getline(char *s, std::streamsize n);/*{return ifstream::getline(s, n, delim);}*/
    std::istream& getline(char *s, std::streamsize n, char delim);
    /*{
        return ifstream::getline(s, n, delim);
    }*/

  protected:
    char findDelimiter();

  private:
    //disable copy-constructor
    InFileStream(const InFileStream &copy);
    std::string filename;
    //auto_ptr<ifstream> inFile;
    char delim;
};

#endif //INFILESTREAM_H


