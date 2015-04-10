/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#include <stdexcept>
#include <exception>
#include <string>
namespace clustalw
{

class VectorOutOfRange : public std::exception
{
    public:
        VectorOutOfRange(std::string vectorName, int index, int max)
            : _name(vectorName), _index(index), _max(max)
        {}
        ~VectorOutOfRange() throw();
        int index(){return _index;}
        int max(){return _max;}
        const char* what() const throw();
        const char* what();
    private:
        std::string _name;
       int _index;
        int _max;
 };

}

