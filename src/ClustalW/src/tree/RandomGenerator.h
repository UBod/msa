/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
/*
 *
 *  RandomGenerator
 *
 *  -linear and additive congruential random number generators
 *  (see R. Sedgewick, Algorithms, Chapter 35)
 *
 *  Implementation: R. Fuchs, EMBL Data Library, 1991
 *
 */
#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

namespace clustalw
{

class RandomGenerator
{
    public:
        /* Functions */
        RandomGenerator(unsigned long s);
        unsigned long addRand(unsigned long r);

        /* Attributes */

    private:
        /* Functions */
        unsigned long mult(unsigned long p, unsigned long q);

        /* Attributes */
        unsigned long j;
        unsigned long a[55];
        unsigned long m;
        unsigned long m1;
};

}
#endif
