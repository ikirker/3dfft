/*
 *  decomposition.c
 *  Calls relating directly to the making of the decomposition
 *   of the data.
 *
 *  Created by Ian Kirker on 18/04/2008.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h> /* For sqrt */
#include "comms.h"
#include "decomposition.h"
#include "A2A3D.h"


/* Set up domain sizes, processor arrangements, and column and row communicators. */
void makeDecomposition(int decompDims[2], int domainSize[2], int extent, int decomp, 
					  int size, int cartCoords[2], ataInfo *rowInfo, ataInfo *colInfo, MPI_Comm *commAll)
{	
	int cartRank;
	int periodicity[2] = {0,0};
	MPI_Comm tempComm;
	
	/* Work out size of domain */
	if ((decomp == 1)||(decomp==0))
	{ 
		decompDims[0] = 1;
		decompDims[1] = size;
	} else 
	if (decomp == 2)
	{
		divide2Ddomain(decompDims, size);
	};

	domainSize[0] = extent / decompDims[0];
	domainSize[1] = extent / decompDims[1];
	
	/* Check for a valid decomposition */
	if ( ( domainSize[0] * decompDims[0] != extent ) || 
	     ( domainSize[1] * decompDims[1] != extent ) )
	{ 
		if (amMaster(MPI_COMM_WORLD))
			fprintf(stderr, "Invalid decomposition obtained - check parameters.\n");
		MPI_Finalize();
		exit(6);
	}

	/* The creation of a cartesian communicator seems a little gratuitous  *
	 *  but it allows us generalisation. */
	MPI_Cart_create ( *commAll, 2, decompDims, periodicity, 1, &tempComm );
	*commAll = tempComm;

	/* Get this processor's position in the grid */
	MPI_Comm_rank(*commAll, &cartRank);
	MPI_Cart_coords(*commAll, cartRank, 2, cartCoords);

	/* Make the column and row communicators for the 2D case */
	/*	int MPI_Comm_split(MPI_Comm comm, int color, int key,
            MPI_Comm *newcomm) */
	MPI_Comm_split(*commAll, cartCoords[1], cartCoords[0], &(rowInfo->comm) );
	MPI_Comm_split(*commAll, cartCoords[0], cartCoords[1], &(colInfo->comm) );
	
	rowInfo->rearrangeDirection = ROWS;
	colInfo->rearrangeDirection = COLS;
	
	return;
}

void divide2Ddomain(int dimensions[2], int processors)
{ /* This divides up an processor count into two dimensions *
   *  and puts the result into dimensions. */
	int i;
	for(i = (int) sqrt( (double) processors);
		i>0;
		i--)
	{
		if ( (0 == i%2) && ( 0 == processors%i ) )
		{
			dimensions[0] = i;
			dimensions[1] = processors/i;
			return;
		}
	}
	fprintf(stderr, "Error in divide2Ddomain - processors = %d.\n"
					" No decomposition could be made.\n", processors);
	exit(6);
}
