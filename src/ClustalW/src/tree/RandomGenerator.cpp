/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include "RandomGenerator.h"

namespace clustalw
{

/**
 * The contructor initialises the random algorithm.
 * @param s 
 * @return 
 */
RandomGenerator::RandomGenerator(unsigned long s)
 : m(100000000), m1(10000)
{
    a[0] = s;
    j = 0;
    do
    {
        ++j;
        a[j] = (mult(31, a[j - 1]) + 1) % m;
    }
    while (j < 54);
}

/**
 * additive congruential method.
 * @param r 
 * @return unsigned long random number in the range 0 to r-1
 */
unsigned long RandomGenerator::addRand(unsigned long r)
{
    int x, y;
    j = (j + 1) % 55;
    x = (j + 23) % 55;
    y = (j + 54) % 55;
    a[j] = (a[x] + a[y]) % m;

    return (((a[j] / m1) *r) / m1);
}

/**
 * 
 * @param p 
 * @param q 
 * @return 
 */
unsigned long RandomGenerator::mult(unsigned long p, unsigned long q)
{
    unsigned long p1, p0, q1, q0;

    p1 = p / m1;
    p0 = p % m1;
    q1 = q / m1;
    q0 = q % m1;
    return ((((p0 *q1 + p1 * q0) % m1) *m1 + p0 * q0) % m);
}

}
