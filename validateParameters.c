/*
 *  validateParameters.c
 *  A number of tests to make sure the options specified are valid.
 *
 *  Created by Ian Kirker on 14/04/2008.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "libDefs.h"
#include "comms.h"
#include "validateParameters.h"

void validateParameters(int size, int extent, int decomp)
{
	int temp;
	int failed = 0;
	
	/* Check for valid number of processors */
	/*  It must be a power of two. */
	temp = 1;
	while(temp<size)
	{
		temp*=2;
	}
	if ( temp != size )
	{
		if (amMaster(MPI_COMM_WORLD))
			fprintf(stderr, "Invalid number of processors specified - %d is not a power of 2.\n", size);
		failed = 1;
	}
	
	/* Check valid decomp - must be either 0, 1 or 2 */
	/* If 3 was passed, it has been altered to 1 in the options retrieval. */
	if ((decomp != 0) && (decomp != 1) && (decomp != 2))
	{
		if (amMaster(MPI_COMM_WORLD))
			fprintf(stderr, "Invalid decomposition specified - "
		                    "use 0 for automatic, 1 for slab, "
				            "2 for rod, 3 for slab with 2D FFT.\n");
		failed = 1;
	}
	
	if (decomp == 0)
	{
		if ( 0 == libraryHasAutomaticDecomposition() )
		{
			if (amMaster(MPI_COMM_WORLD))
				fprintf(stderr, "Invalid decomposition specified - "
			                    " this library does not support an"
				                " automatic decomposition.\n" );
			failed = 1;
		}
	}
	
	
	/* Check valid extent */
	
	/* In a slab decomposition, the extent must divide by the number of processors. */
	if (decomp == 1)
	{
		if ( size * (extent/size) != extent )
		{
			if (amMaster(MPI_COMM_WORLD))
				fprintf(stderr, "Invalid extent specified - extent must be a multiple of processor count.\n");
			failed = 1;
		}
	}
	
	/* In a rod decomposition, the extent^2 must divide by the number of processors, *
	 *  but in a way that makes sure we have rectangular domains across one face of  *
	 *  the data cube. Because of our size = 2^n restriction, this means that our    *
	 *  extent only has to be divisible by 2. */
	if (decomp == 2)
	{
		if ( 2 * (extent/2) != extent ) 
		{
			if (amMaster(MPI_COMM_WORLD))
				fprintf(stderr, "Invalid extent specified - extent must be divisible by 2.\n");
			failed = 1;
		}
	};
	
	if (failed == 1)
	{
		commsEnd();
		exit(2);
	}
	
	return;
};
