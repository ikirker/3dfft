/*
 *  performLocalTranspose.c
 *  Transposes numberOfSlabs slabs of extent*extent sized
 *   contiguous arrays of complexType data.
 *
 *  Created by Ian Kirker on 15/04/2008.
 *
 */

#include "libDefs.h"
#include "performLocalTranspose.h"

/* Transposes multiple extent*extent 2D complex arrays stored contiguously in memory. */
void performLocalTranspose(complexType *data, int extent, int numberOfSlabs)
{
	int i,j,k;
	for(i=0;i<numberOfSlabs;i++)
	{
		for(j=0;j<extent;j++)
		{
			for(k=0;k<j;k++)
			{
				complexSwap(&data[i*extent*extent + j*extent + k], &data[i*extent*extent + k*extent + j]);
			}
		}
	}
}
