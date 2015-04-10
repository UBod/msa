/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <algorithm>
#include "Sequence.h"

using namespace std;

namespace clustalw
{

/**
 * 
 * @param seq 
 * @param name 
 * @param title 
 * @return 
 */
Sequence::Sequence(string& seq, string& name, string& title)
{
    copyStringIntoVector(&_sequence, &seq);
    encodeSequence();
    _name = name;
    _title = title;
    identifier = utilityObject->getUniqueSequenceIdentifier();
}

Sequence::Sequence(std::string& seq, std::string& name, std::string& title, unsigned long id)
{
    copyStringIntoVector(&_sequence, &seq);
    encodeSequence();
    _name = name;
    _title = title;
    identifier = id;
}

/**
 * This is an overloaded contructor that is used to construct a seq object from an
 * encoded sequenced instead of a string.
 * @param encodedSequence 
 * @param name 
 * @param title 
 * @param id The unique identifier from the previous sequence!!!
 * @return 
 */
Sequence::Sequence(std::vector<int>* encodedSequence, std::string& name, std::string& title,
                   unsigned long id)
{
    _encodedSequence = *encodedSequence;
    _name = name;
    _title = title;
    identifier = id;        
}

/**
 * 
 */
void Sequence::encodeSequence()
{
     /* code seq as ints .. use gapPos2 for gap */
    std::vector<char>::iterator it;

    _encodedSequence.push_back(0);
    
    for(it = _sequence.begin(); it != _sequence.end(); ++it) 
    {
        if (*it == '-')
        {
            _encodedSequence.push_back(userParameters->getGapPos2());
        }
        else
        {
            _encodedSequence.push_back(userParameters->resIndex(
                                   userParameters->getAminoAcidCodes(), *it));
        }
    }
}

/**
 * 
 * @param _vectorTo 
 * @param _stringFrom 
 */
void Sequence::copyStringIntoVector(vector<char>* _vectorTo, string* _stringFrom)
{
    _vectorTo->clear();

    for(int i = 0; i < (int)_stringFrom->size(); i++)
    {
        _vectorTo->push_back(_stringFrom->at(i));
    }
    
    if(_vectorTo->size() != _stringFrom->size())
    {
        std::cerr << "Error: In function copyStringIntoVector. Strings different length!\n";
        throw 1;
    }
}

/**
 * 
 */
void Sequence::printSequence()
{
    std::cout << "This is the sequence and the encoded sequence " << _name << std::endl;
    
    std::vector<char>::iterator itChar;    
    for(itChar = _sequence.begin(); itChar != _sequence.end(); ++itChar)
    {
        cout << *itChar;
    } 
    cout << std::endl;
        
    std::vector<int>::iterator itInt;
    for(itInt = _encodedSequence.begin(); itInt != _encodedSequence.end(); ++itInt)
    {
        cout << "  " << *itInt;
    } 
    cout << std::endl;
}

/**
 * 
 */
void Sequence::checkIntegrity()
{
    // The sequences should be the same length.
    if(_sequence.size() != _encodedSequence.size())
    {
        std::cerr << "Error: _sequence is not same size as _encodedSequence\n";
        throw 1;
    }
}

/**
 * 
 * @return the encoded sequence, this is what is used in the pairwise!
 */
std::vector<int>* Sequence::getSequence()
{
    return &_encodedSequence;
}

/**
 * 
 * @return 
 */
std::string Sequence::getName()
{
    return _name;
}

/**
 * 
 * @return 
 */
std::string Sequence::getTitle()
{
    return _title;
}

/**
 * 
 * @return 
 */
bool Sequence::isEmpty()
{
    if(_sequence.size() == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
   
}

/**
 * 
 * @return 
 */
bool Sequence::checkDNAFlag()
// check if DNA or Protein
// The decision is based on counting all A,C,G,T,U or N.
// If >= 85% of all characters (except -) are as above => DNA 
{
    int c, numResidues, numBases;
    float ratio;
    string dna_codes = "ACGTUN";

    numResidues = numBases = 0;
    
    vector<char>::iterator seqIterator = _sequence.begin();
    
    while (seqIterator != _sequence.end())
    {
        if (*seqIterator != '-')
        {
            numResidues++;
            if (*seqIterator == 'N')
            {
                numBases++;
            }
            else
            {
                c = userParameters->resIndex(dna_codes, *seqIterator);
                if (c >= 0)
                {
                    numBases++;
                }
            }
        }
        seqIterator++;
    }
    
    if ((numBases == 0) || (numResidues == 0))
    {
        return false;
    }
    ratio = (float)numBases / (float)numResidues;

    if (ratio >= 0.85)
    {
        return true;
    }
    else
    {
        return false;
    }
}

}

