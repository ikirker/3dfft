/*
 *  options.c
 *  Gets the options from the arguments.
 *
 *  Created by Ian Kirker on 15/05/2008.
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <unistd.h>

#include "libDefs.h"
#include "options.h"


int getOptions(int *argc, char ***argv, int *extent, int *decompType, int *use2DFFT, int *skip, int *skipFFT, int *printOut)
{   /* Get command-line options */
	
	int c;
	
	/* Defaults for testing. */
	*decompType = 1;
	*extent = 4;
	
	
	opterr = 0; /* Defined in unistd.h */
	while ((c = getopt (*argc, *argv, "x:d:lnhfp")) != -1)
	{
		switch (c)
		{
			/* argument -x sets the extent */
			case 'x':
			 *extent = atoi(optarg);
			 break;
			 
			/* -d sets type of decomposition */
			/* 0 indicates an automatic decomposition - *
			 * 1 indicates slab decomposition           * 
			 * 2 indicates rod decomposition            *
			 * 3 indicates slab using a 2D FFT call.    */
			case 'd':
			 *decompType = atoi(optarg);
			 if (*decompType == 3)
			 {
				*decompType = 1;
				*use2DFFT = 1;
			 }
			 break;
			 
			/* -l prints the FFT library used */
			case 'l':
			 printLib();
			 exit(0);
			 break;
			 
			/* -n makes the program skip all the actual work */
			case 'n':
			 *skip = 1;
			 break;
			 
			/* -f makes the program skip all the ffts */
			case 'f':
			 *skipFFT = 1;
			 break; 
			 
			/* -p prints data after operation rather than checking it */ 
			case 'p':
			 *printOut = 1;
			 break; 
			 
			/* Prints a list of options */ 
			case 'h':
			 printOptionList();
			 exit(0);
			 break;
			  
			/* Errant option handler */
			case '?':
			 if ((optopt == 'x')||(optopt == 'd'))
			 {
			  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			  exit(1);
			 }
			 else if (isprint (optopt))
			 {
			  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			  exit(1);
			 }
			 else
			 {
			  fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
			  exit(1);
			 }
			 
			default:
			 abort ();
		} /* end switch */
	} /* end while */
	return 0;
}

void printOptionList()
{
	printf("3D FFT benchmark options: \n"
	       "  -x<number>     Sets the size on one side of the global data cube.\n"
		   "  -d[0|1|2|3]    Sets the type of decomposition used:\n"
		   "                   0 - automatic (not universally available)\n"
		   "                   1 - slab\n"
		   "                   2 - rod\n"
		   "                   3 - slab with 2D FFTs used on each slab\n"
		   "  -f             Skips all FFT steps.\n"
		   "  -n             Skips all FFT and communication steps.\n"
		   "  -p             Prints data instead of checking.\n"
		   "  -h             Prints this message.\n"
		   );
}

