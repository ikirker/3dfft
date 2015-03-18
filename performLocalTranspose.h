/*
 *  performLocalTranspose.h
 *  Prototype for performLocalTranspose
 *
 *  Created by Ian Kirker on 15/04/2008.
 *
 */

/* Performs a transposition of a number of slabs on a cube of *
 *  complex numbers organised as follows:                     *
 *
 *     3   7   11  15/
 *    2   6   10  14//
 *   1   5   9   13///
 *  0   4   8   12////
 *  16  20  24  28///
 *  32  36  40  44//
 *  48  52  56  60/
 */
#ifndef HEADER_PERFORMLOCALTRANSPOSE

#include "libDefs.h" /* For complexType and complexSwap definition */

void performLocalTranspose(complexType *data, int extent, int numberOfSlabs);

#define HEADER_PERFORMLOCALTRANSPOSE
#endif