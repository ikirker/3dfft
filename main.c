/**************************************
 * 3D FFT Benchmark                   *
 *                                    *
 *                                    *
 *  Ian Kirker ian.kirker@gmail.com   *
 **************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include <ctype.h>

#include "A2A3D.h"
#include "comms.h"
#include "dataOps.h"
#include "decomposition.h"
#include "libDefs.h"
#include "options.h"
#include "performLocalTranspose.h"
#include "validateParameters.h"

#define TOLERANCE 1e-10

int main (int argc, char ** argv) {

	/* Double buffer data - AlltoAll cannot be performed in-place */
	complexType *data[2];
	
	int extent;         /* Size of whole problem cube (extent*extent*extent)    */
	int domainSize[2];  /*  and per processor along each decomposable dimension */
	
	char decompName[5]; /* For output string */
	int decomp;         /* Decomposition type - 1 for slab, 2 for rod */
	int use2DFFT = 0;   /* 1 if we're using the library's 2D FFT, otherwise 0 */
	
	int skip    = 0;    /* Skip all work */
	int skipFFT = 0;    /* Skip FFTs, just ATA */
	int printOut= 0;    /* Print out the data instead of checking it at the end */
	
	double phaseTime[6]; /* Tracks time for each phase of FFT */
	
	/*** MPI Variables ***/
	int size;          /* Global number of tasks */
	int cartCoords[2]; /* Coordinates within the Cartesian communicator */
	int decompDims[2]; /* Number of processors along each dimension of the decomp */
	
	MPI_Comm commAll = MPI_COMM_WORLD;
	
	ataInfo ataRow, ataCol; /* Stored All-to-All information */
	
	/********* Preparation **********/
	
	
	/* Fire up the MPI handler and set standard variables. */
	commsInit(&argc, &argv);
	size = getSize(commAll);
	
	/* Get Command Line Options */
	getOptions(&argc, &argv, &extent, &decomp, &use2DFFT, &skip, &skipFFT, &printOut);
	
	/* Check all the parameters before going ahead */
	validateParameters(size,extent,decomp);

	/* Prepares a whole bunch of stuff -            */
	makeDecomposition(decompDims, domainSize, extent, decomp, 
					  size, cartCoords, &ataRow, &ataCol, &commAll);

	
	/* Create the data to use for the FFTs.         *
	 *  prepareFFTs may destroy some of the data    *
	 *  in the array, so we populate it afterwards. */           
	makeDataArrays(data, extent, domainSize);
	prepareFFTs(data[0], decomp, use2DFFT, extent, domainSize, ataCol.comm);
	
	if ( ( skipFFT==1 ) || ( skip==1 ) )
	{ /* If we're skipping bits, use the test data. */
		makeTestData(data, extent, domainSize, cartCoords);	
	} else {
		makeData(data, extent, domainSize, cartCoords);	
	}
	
	/* Print out the job parameters in a human understandable format and continue */
	if (amMaster(commAll))
	{
		if (decomp == 0) sprintf(decompName, "%s", "auto");
		if (decomp == 1) sprintf(decompName, "%s", "slab");
		if (decomp == 2) sprintf(decompName, "%s", "rod");
		fprintf(stderr, 
			"Running MPI 3D FFT Benchmark with %d processors.\n"
			" Problem size:  \t%dx%dx%d\n"
			" Decomposition: \t%s: %dx%d\n"
			" Each array:    \t%dx%dx%d\n"
			" Library:       \t%s\n"
			" Using 2D FFT call: \t%s\n",
			size,
			extent,extent,extent,
			decompName,
			decompDims[0],decompDims[1],
			domainSize[1],domainSize[0],extent,
			FFT_NAME,
			((use2DFFT==1)?"yes":"no")
			);
		if (skip == 1)
			fprintf(stderr, " Skip is set, calculation will be skipped.\n");
		if (skipFFT == 1)
			fprintf(stderr, " SkipFFT is set, 1D FFTs will be skipped.\n");
	}
	
	/* Barrier before we start */
	commSync(commAll);
	
	/********* Actual FFTs **********/
	
	/* If !0, skip skips the whole operation, skipFFT skips any transforms */
	
	phaseTime[0] = MPI_Wtime();
	if ( (decomp == 1) && (skip == 0) )
	{ /* Slab type decomp */
		if ( use2DFFT == 1 )
		{ /* With 2D FFT types in the slab dimensions */
			/* Note - data operated on this way may be transposed. */
			if (!skipFFT) perform2DFFT(data[0], data[1], extent, domainSize);
			phaseTime[1] = MPI_Wtime();
			phaseTime[2] = phaseTime[1];
			phaseTime[3] = phaseTime[1];
		} 
		else
		{ /* With 1D FFT types in the slab dimensions */
			if (!skipFFT) performFFTset(data[0], data[1], extent, domainSize);
			phaseTime[1] = MPI_Wtime();
			performLocalTranspose(data[0], extent, domainSize[1]);
			phaseTime[2] = MPI_Wtime();
			if (!skipFFT) performFFTset(data[0], data[1], extent, domainSize);
		}
		
		phaseTime[3] = MPI_Wtime();
		
		performDistTranspose(data[0], data[1], domainSize, extent, &ataCol);
		
		if (!skipFFT) performFFTset(data[0], data[1], extent, domainSize);
		
		phaseTime[4] = MPI_Wtime();
		
	} else if ( (decomp == 2) && (skip == 0) ) { 
		/* Rod decomp */
		if (!skipFFT) performFFTset(data[0], data[1], extent, domainSize);
		
		phaseTime[1] = MPI_Wtime();

		performDistTranspose(data[0], data[1], domainSize, extent, &ataRow);
				
		phaseTime[2] = MPI_Wtime();
				
		if (!skipFFT) performFFTset(data[0], data[1], extent, domainSize);
		
		phaseTime[3] = MPI_Wtime();
		
		performDistTranspose(data[0], data[1], domainSize, extent, &ataCol);
		
		phaseTime[4] = MPI_Wtime();
		
		if (!skipFFT) performFFTset(data[0], data[1], extent, domainSize);
	} else if ( (decomp == 0) && (skip == 0) && ( skipFFT == 0 ) ) { 
		/* Automatic Decomp */
		phaseTime[1] = phaseTime[0];
		phaseTime[2] = phaseTime[0];
		phaseTime[3] = phaseTime[0];
		phaseTime[4] = phaseTime[0];
		
		performAutomatic3DFFT(data[0], data[1], extent, domainSize);
	}
	phaseTime[5] = MPI_Wtime();


	/********* Output and finalisation **********/
	
    if ( printOut == 1 )
    { /* If we're skipping bits or requesting it, print the data instead. */
        printData( data, extent, domainSize, decompDims, cartCoords );
    } else if ( ( skipFFT==1 ) || ( skip==1 ) ) {
        if (amMaster(commAll)) {
            fprintf(stderr, "Skipping data checking because some steps have been skipped.\n");
        }
    } else {
        checkData( data, extent, domainSize, cartCoords, TOLERANCE, commAll );
    }
	
	/* Print out computer readable (CSV) job result string */
	if (amMaster(commAll))
	{
		printf("fft-results:%d,%d,%s,%s,%s,%g,%g,%g\n",
			size,
			extent,
			decompName,
			((use2DFFT==1)?"2DFFT":"1DFFT"),
			FFT_NAME,
			
			/* Communication/memory reorg time */ 
			(phaseTime[2] - phaseTime[1]) + 
			 (phaseTime[4] - phaseTime[3]),
			
			/* FFT time */
			(phaseTime[1] - phaseTime[0]) + 
			 (phaseTime[3] - phaseTime[2])+
			 (phaseTime[5] - phaseTime[4]),
			 
			 /* Total time */
			phaseTime[5] - phaseTime[0]
			);
	}
	
	/* Clean up all the parts */
	cleanUpData(data);
	cleanUpFFTs(decomp);
	freeATAcommsHandles(&ataRow, &ataCol);
	commsEnd();
	
	exit(0);
}

