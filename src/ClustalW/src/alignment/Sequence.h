/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * This class contains the sequence information. It should also contain
 * the sequence name, and the encoded version.
 * A vector of Sequences is passed to the Alignment object to set it up.
 * CHANGES:
 * Mark 22-1-2007: Added a unique sequence identifier to help with correct output
 * order of sequences.   
 */
 
#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <vector>
#include <string>
#include "../general/userparams.h"
#include "../general/utils.h"

namespace clustalw
{

class Sequence
{
    public:
        /* Functions */
        Sequence(std::string& seq, std::string& name, std::string& title);
        Sequence(std::string& seq, std::string& name, std::string& title, 
                 unsigned long id);
        Sequence(std::vector<int>* encodedSequence, std::string& name, std::string& title, 
                 unsigned long id);
        void encodeSequence();
        void printSequence();
        std::vector<int>* getSequence();
        bool isEmpty();
        std::string getName();
        std::string getTitle();
        bool checkDNAFlag();
        unsigned long getIdentifier(){return identifier;}
        /* Attributes */

    private:
        /* Functions */
        void checkIntegrity();
        void copyStringIntoVector(std::vector<char>* _vectorTo, std::string* _stringFrom);
        
        /* Attributes */
        std::vector<char> _sequence;
        std::vector<int> _encodedSequence;
        std::string _name;
        std::string _title;
        unsigned long identifier;
};

}
#endif
