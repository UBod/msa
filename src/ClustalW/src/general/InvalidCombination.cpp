/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/**
 * This class will be used to report invalid combinations of residue type and 
 * aligntype. Correct values are either 0 or 1.
 * Note: It will not be possible for a user to cause this exception, only programmer.
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <ostream>
#include <iostream>

namespace clustalw
{
    class InvalidCombination
    {
        public:
            InvalidCombination(int alignResidueType, int alignType)
                                : _alignResidueType(alignResidueType),
                                  _alignType(alignType) {}
            void whatHappened(std::ostream &os = std::cerr)
            {
                os << "Incorrect Combination of alignResidueType and alignType.\n"
                   << "Values should be 0 or 1\n"
                   << "alignResidueType = " << _alignResidueType << "\n"
                   << "alignType = " << _alignType << "\n";
            }
        private:
            int _alignResidueType;
            int _alignType;
    };
}


