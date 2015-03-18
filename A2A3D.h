/*
 *  A2A3D.h
 *  Calls used to perform the all-to-all transposes.
 *
 *  Created by Ian Kirker on 13/05/2008.
 *
 */

#ifndef HEADER_A2A3D
#define HEADER_A2A3D

#include <mpi.h>
#include "libDefs.h"

/* Direction of all to all indicator - goes in rearrangeType */
#define ROWS 0
#define COLS 1

/* Encapsulated data for All-to-All information */
typedef struct { MPI_Comm comm; int rearrangeDirection; } ataInfo;

int performDistTranspose(complexType *data, complexType *dataBuffer, int domainSize[2], int extent,
                         ataInfo *thisATA);

void ataRowRearrange(complexType *dataIn, complexType *dataOut, int domainSize[2], int extent);
void ataColRearrange(complexType *dataIn, complexType *dataOut, int domainSize[2], int extent);
void ataRowUnpack(complexType *dataIn, complexType *dataOut, int domainSize[2], int extent);
void ataColUnpack(complexType *dataIn, complexType *dataOut, int domainSize[2], int extent);

void freeATAcommsHandles(ataInfo *ataRow, ataInfo *ataCol);

#endif
