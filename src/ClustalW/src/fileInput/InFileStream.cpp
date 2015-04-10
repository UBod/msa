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
 * Note: This is an ugly workaround; at present various operations repeatedly
 * construct/destruct an instance and open/close a sequence file up to 12
 * times! A cleaner class will probably derive this class from something like
 * 'istream' aggregating a 'filebuf' under control of istream::seekg().
 *
 * Created: 09-02-07,Nigel Brown(EMBL)
 * 
 * Changes:
 * Mark Larkin 13-2-07: I removed the dynamic cast from the getline functions. 
 ***************************************************************************/
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <string>
#include <fstream>
#include <iostream>
#include "InFileStream.h"
using namespace std;

const char LF = 0x0a;  //linefeed
const char CR = 0x0d;  //carriage return

InFileStream::InFileStream() :
    ifstream()
{
    delim = '\n'; // default
    //cout << "InFileStream() constructor 1" << endl;
}

InFileStream::InFileStream(const char *filename) :
    ifstream(filename, ios::in), filename(filename)
{
    //cout << "InFileStream(f) constructor 2" << endl;
    delim = findDelimiter();
}

//- copy-constructor: can't copy superclass private members
//- InFileStream::InFileStream(const InFileStream &copy) :
//-     ifstream(static_cast<const ifstream&>(copy))
//- {
//-     cout << "InFileStream() constructor 3" << endl;
//-     delim = copy.delim;
//- }

void InFileStream::open(const char *filename) 
{

    this->filename = filename;
    ifstream::open(filename, ios::in);
    if  (ifstream::fail())
        return;
    delim = findDelimiter();
}

//not necessary, but for symmetry to open()
void InFileStream::close() 
{
    ifstream::close();   
}


//getline with stored delimiter
std::istream& InFileStream::getline(char *s, streamsize n) 
{
    return ifstream::getline(s, n, delim);
}

//getline with caller supplied delimiter
std::istream& InFileStream::getline(char *s, streamsize n, char delim) 
{
    return ifstream::getline(s, n, delim);
}


/**
 * Mark 24-1-2007. I added the function findDelimiter to determine if '\r' or
 * '\n' will be used as the line delimiter when parsing the file.
 *
 * 25-01-07,Nigel Brown(EMBL): changed body of loop to check successive chars
 * in case of DOS/Windows
 *
 * 09-02-07,Nigel Brown(EMBL): moved member into new InFileStream subclassed
 * from std::ifstream, so this is called automatically for any file reader
 * that uses InFileStream in place of std::ifstream. Replaced if/then/else
 * with switch.
 */
char InFileStream::findDelimiter()
{
    ifstream in;
    int type = 0;
    
    in.open(filename.c_str(), ios::in);
    if (in.fail())
        return delim;
    
    in.seekg(0, ios::beg);

    //look for CR or LF or CRLF (or LFCR)
    if (in.is_open()) {
        char c;
        while (in.get(c)) {
            if (c == CR)
                type |= 1;
            else if (c == LF)
                type |= 2;
            else if (type)
                break;
        }
    }
    in.close();

    switch (type) {
	case 1:
	    //cout << "file is Mac System 9" << endl;
	    delim = '\r';
	    break;
	case 2:
	    //cout << "file is UNIX" << endl;
	    delim = '\n';
	    break;
	case 3:
	    //cout << "file is DOS" << endl;
	    delim = '\n';
	    break;
	default: //short or empty file
	    //cout << "file is UNIX (default)" << endl;
	    delim = '\n';
    }
    return delim;
}


