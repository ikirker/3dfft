/*
 *  dataOps.c
 *  Making and checking the data before and after operation, respectively.
 *  The complex functions used in this file should really be inlined using IPA.
 *
 *  Created by Ian Kirker on 16/05/2008.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libDefs.h"
#include "comms.h"
#include "dataOps.h"


void makeDataArrays( complexType *data[2], int extent, int domainSize[2] )
{
	/* Allocate storage space, checking for NULLs */
	/* Avoid this failing -- core dumps break IO handlers */
	if ( NULL == ( data[0] = malloc( extent * domainSize[0] * domainSize[1] * sizeof(complexType) ) ) )
	{
		fprintf(stderr, "Could not allocate primary data array.\n");
		commsEnd();
		exit(5);
	}
	
	if ( NULL == ( data[1] = malloc( extent * domainSize[0] * domainSize[1] * sizeof(complexType) ) ) )
	{
		fprintf(stderr, "Could not allocate secondary data array.\n");
		commsEnd();
		exit(5);
	}
}


void makeData( complexType *data[2], int extent, int domainSize[2], int cartCoords[2] )
{ /* Fills data array 0 with a trivariate multisine function. This should ideally give *
   *  a transform output that is easy to verify.  */
	int i,j,k;
	double ratio = 2.0 * 3.14159265358979323846 / ( (double) extent );
	
	/* Populate data field */
	for(i=0;i<domainSize[1];i++)
	{
		for(j=0;j<domainSize[0];j++)
		{
			for(k=0;k<extent;k++)
			{
				complexSet( &data[0][ i*domainSize[0]*extent + j*extent + k ],
							sin(
							ratio * (i + ( cartCoords[1] * domainSize[1] ) ) +
							ratio * (j + ( cartCoords[0] * domainSize[0] ) ) + 
							ratio * k
							)
							,
							0);
			}
		}
	}
	
	return;
	
}

void makeTestData( complexType *data[2], int extent, int domainSize[2], int cartCoords[2] )
{ /* Fills data array 0 such that the decomposed cube contains a simple counting up in the real *
   *  part, and the processor location in the imaginary part. For testing. */
	int i,j,k;
		
	/* Populate data field */
	for(i=0;i<domainSize[1];i++)
	{
		for(j=0;j<domainSize[0];j++)
		{
			for(k=0;k<extent;k++)
			{
				complexSet( &data[0][ i*domainSize[0]*extent + j*extent + k ],
							( (i + ( cartCoords[1] * domainSize[1] ) ) * extent * extent )  +
							( (j + ( cartCoords[0] * domainSize[0] ) ) * extent )  + k,
							cartCoords[0] * 100 + cartCoords[1]);
			}
		}
	}
	
	return;
	
}


int printData( complexType *data[2], int extent, int domainSize[2], int decompDims[2], int cartCoords[2] )
{ /* Simply prints out the data array in an understandable (though not always clean) format. */
  /* Best for extent<8 */
	_Complex double z;
	int i,j,k,m,n;
	for (m=0;m<decompDims[0];m++)
	{
		for (n=0;n<decompDims[1];n++)
		{
			if ( ( cartCoords[0] == m ) && ( cartCoords[1] == n ) )
			{
				printf("(%d,%d)\n",cartCoords[0],cartCoords[1]);
	
				for(i=0;i<domainSize[1];i++)
				{
					for(j=0;j<domainSize[0];j++)
					{
						for(k=0;k<extent;k++)
						{
							z = complexNative( data[0][ i*domainSize[0]*extent + j*extent + k ] );
							printf("%g,%g ", (abs(creal(z))>0.00001)?creal(z):0.0,(abs(cimag(z))>0.00001)?cimag(z):0.0);
						}
						printf("\n");
					}
					printf("\n\n");
				};
	
				fflush(stdout);
				MPI_Barrier(MPI_COMM_WORLD);
			} else { MPI_Barrier(MPI_COMM_WORLD); }
		}
	}
	return 1;
}

int checkData( complexType *data[2], int extent, int domainSize[2], int cartCoords[2], double tolerance, MPI_Comm comm )
{ /* Verifies that two peaks are in far corner and one off top near corner of array, *
   *  and that all other values are equal to zero.                                   */
	int i,j,k;
	double residue=0;
	double peaksize;
	
	/* Generate comparison data, first filling comparison array (data[1]) with zeroes... */
	for(i=0;i<domainSize[1];i++)
	{
		for(j=0;j<domainSize[0];j++)
		{
			for(k=0;k<extent;k++)
			{
				complexSet(data[1]+ i*domainSize[0]*extent + j*extent + k, 0, 0);
			}
		}
	}
	
	/* And then setting the peaks. First we have to find out where element 1,1,1 is, though. */
	/* The near peak is set to -i * 0.5 * extent^3, the far to i*0.5*extent^3                */
	/* NB: Cast these all to doubles so that we never need to worry about integer overflow   *
	 *  mid-multiply.                                                                        */
	peaksize = 0.5 * (double)extent * (double) extent * (double) extent;
	if ( domainSize[0] > 1 ) 
	{
		if ( cartCoords[0] == 0 )		
		{
			if ( domainSize[1] > 1 ) 
			{
				if ( cartCoords[1] == 0 )
				{   /* If it's on processor (0,0), set element 1,1,1. */
					complexSet(data[1] + domainSize[0]*extent + extent + 1, 0, -1 * peaksize);
				}
			} else {
				if ( cartCoords[1] == 1 )
				{  /* If it's on processor (0,1), set element 0,1,1 */
					complexSet(data[1] + extent + 1, 0, -1 * peaksize);
				}
			}
		}
	} else  if ( cartCoords[0] == 1 ) 
	{
		if ( domainSize[1] > 1 ) 
		{
			if ( cartCoords[1] == 0 )
			{  /* If it's on processor (1,0), set element 1,0,1 */
				complexSet(data[1] + domainSize[0]*extent + 1, 0, -1 * peaksize);
			}
		} else {
			if ( cartCoords[1] == 1 )
			{ /* If it's on processor (1,1), set element 0,0,1 */
				complexSet(data[1] + 1, 0, -1 * peaksize);
			}
		}
	}
	
	/* And then the far peak */
	if ( ( cartCoords[0] == extent/domainSize[0] - 1 ) && ( cartCoords[1] == extent / domainSize[1] - 1 ) )
	{ /* If I'm the bottom right processor */
		complexSet(data[1] + domainSize[1]*domainSize[0]*extent - 1 , 0, peaksize);
	}
	
	/* Now generate the sum of the absolute differences between the two... */
	for(i=0;i<domainSize[1];i++)
	{
		for(j=0;j<domainSize[0];j++)
		{
			for(k=0;k<extent;k++)
			{
				residue += complexAbsNorm(*(data[1]+ i*domainSize[0]*extent + j*extent + k), 
				                          *(data[0]+ i*domainSize[0]*extent + j*extent + k));
			}
		}
	}
	
	/* Reduce over processors. */
	doubleGlobalSum(&residue, comm);

	/* Normalise for matrix size. */
	/* Cast to doubles as before to eliminate problems with integer overflow mid-multiply. */
	residue/= (double) extent * (double) extent * (double) extent;

	if (amMaster(comm))
		fprintf(stderr, "Residue = %g\n", residue);
	
	/* And then if the residue is outwith acceptable limits, terminate without a result line. */
	if ( residue < tolerance ) 
	{
		return 1; 
	} else {
		commsEnd();
		exit(1);
	}
}

void cleanUpData(complexType *data[2])
{
	free(data[0]);
	free(data[1]);
}
