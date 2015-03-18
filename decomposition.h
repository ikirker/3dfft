/*
 *  decomposition.h
 *  Calls relating directly to the making of the decomposition
 *   of the data.
 *
 *  Created by Ian Kirker on 18/04/2008.
 *
 */

#ifndef HEADER_DECOMPOSITION

#include <mpi.h>
#include "A2A3D.h"

/* ataInfo struct defined in A2A3D.h */
void makeDecomposition(int decompDims[2], int domainSize[2], int extent, int decomp, 
					  int size, int cartCoords[2], ataInfo *rowInfo, ataInfo *colInfo, MPI_Comm *commAll);
					  				  					  
void divide2Ddomain(int dimensions[2], int processors);

#define HEADER_DECOMPOSITION
#endif