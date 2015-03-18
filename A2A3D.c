/*
 *  A2A3D.c
 *  Calls used to perform the all-to-all transposes.
 *
 *  Created by Ian Kirker on 13/05/2008.
 *
 */


#include "libDefs.h"
#include "A2A3D.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int performDistTranspose(complexType *data, complexType *dataBuffer, int domainSize[2], int extent,
                         ataInfo *thisATA)
{ /* Performs the whole tranpose, all to all, rearranging etc. Called from main.c */
	int elements;
	int err;
	
	if (thisATA->rearrangeDirection == ROWS)
	{
		ataRowRearrange(data, dataBuffer, domainSize, extent);
		elements = domainSize[0] * domainSize[0] * domainSize[1];
	} else if (thisATA->rearrangeDirection == COLS)
	{
		ataColRearrange(data, dataBuffer, domainSize, extent);
		elements = domainSize[0] * domainSize[1] * domainSize[1];
	}
	
	err = MPI_Alltoall(dataBuffer, elements * 2, MPI_DOUBLE, 
	             data, elements * 2, MPI_DOUBLE, thisATA->comm);
	
	if (thisATA->rearrangeDirection == ROWS)
	{
		ataRowUnpack(data, dataBuffer, domainSize, extent);
	} else if (thisATA->rearrangeDirection == COLS)
	{
		ataColUnpack(data, dataBuffer, domainSize, extent);
	}
	
	memcpy(data,dataBuffer,domainSize[0]*domainSize[1]*extent*sizeof(complexType));
	
	return 0;
}

void ataRowRearrange(complexType *dataIn, complexType *dataOut, int domainSize[2], int extent)
{ /* Rearranges the data in a domain such that all the data that needs to be *
   *  sent to one processor is contiguous and in the right order, for an     *
   *  all-to-all across rows of a 2D decomposition of a 3D array.            */
  /* The numbers in comments below refer to domainSize[] = {2,3}, extent=12  */ 
	int i;
	
	for(i=0;i<domainSize[0]*domainSize[1]*extent;i++)
	{
		complexAssign(&dataOut[
		                       (  ( i % domainSize[0] ) * domainSize[0] ) + // 0 3 6 every 1
							   (  ( ( i % extent ) / domainSize[0] ) * domainSize[0] * domainSize[0] * domainSize[1] ) + // 0 18 36 54 every 3
							   (  ( i / extent ) % domainSize[0] ) + // 0 1 2 every 12
							   (  ( i / ( domainSize[0] * extent ) ) * domainSize[0] * domainSize[0] ) // 0 9 every 36
							  ]
							  , dataIn[i]
							  );
	}
}

void ataColRearrange(complexType *dataIn, complexType *dataOut, int domainSize[2], int extent)
{
	int i;
	
	for(i=0;i<domainSize[0]*domainSize[1]*extent;i++)
	{
		
		complexAssign(&dataOut[
		                       ( i % extent ) * domainSize[0] * domainSize[1] +
							   ( ( i / extent ) % domainSize[0] ) * domainSize[1] +
							   ( i / ( domainSize[0] * extent ) )
							  ],dataIn[i]
							  );
	}
}


void ataRowUnpack(complexType *dataIn, complexType *dataOut, int domainSize[2], int extent)
{ /* Unpacks the data after the all-to-all. Performs the same operation as receiving *
   *  with a vector type would, but allows more flexibility, esp. in the case of the *
   *  row-wise. */
	int i;
	for(i=0;i<domainSize[0]*domainSize[1]*extent;i++)
	{
		complexAssign(&dataOut[
		                      (i%domainSize[0]) + 
							  ( ((i/domainSize[0]) % (domainSize[0] * domainSize[1])) * extent ) +
							  ( ( i / (domainSize[0] * domainSize[0] * domainSize[1] )) * domainSize[0] )
		                      ],dataIn[i]
							  );
	}
}

void ataColUnpack(complexType *dataIn, complexType *dataOut, int domainSize[2], int extent)
{
	int i;
	for(i=0;i<domainSize[0]*domainSize[1]*extent;i++)
	{
		complexAssign(&dataOut[
							  (i%domainSize[1]) + 
							  ( ((i/domainSize[1]) % (domainSize[0] * domainSize[1])) * extent ) +
							  ( ( i / (domainSize[0] * domainSize[1] * domainSize[1] )) * domainSize[1] )
		                      ],dataIn[i]
							  );
	}
}

void freeATAcommsHandles(ataInfo *ataRow, ataInfo *ataCol)
{
	MPI_Comm_free(&ataRow->comm);
	MPI_Comm_free(&ataCol->comm);
}

