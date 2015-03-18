/*
 *  dataOps.h
 *  Making and checking the data before and after operation, respectively.
 *
 *  Created by Ian Kirker on 16/05/2008.
 *
 */

#include "libDefs.h"
#ifndef HEADER_DATAOPS

void makeDataArrays( complexType *data[2], int extent, int domainSize[2] );
int printData( complexType *data[2], int extent, int domainSize[2], int decompDims[2], int cartCoords[2] );
int checkData( complexType *data[2], int extent, int domainSize[2], int cartCoords[2], double tolerance, MPI_Comm comm );
void makeData( complexType *data[2], int extent, int domainSize[2], int cartCoords[2] );
void makeTestData( complexType *data[2], int extent, int domainSize[2], int cartCoords[2] );
void cleanUpData(complexType *data[2]);

#define HEADER_DATAOPS
#endif