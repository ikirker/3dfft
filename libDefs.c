/*
 *  libDefs.c
 *  This file and its associated header contain any functions and statements 
 *   (within the actual program) that need to be modified when adding a 
 *   new library. This includes any work done on complex numbers, FFT
 *   planning and execution, and would also contain the definition of
 *   the appropriate MPI type to represent the complex number, except
 *   that I decided to assume they were all just two MPI_DOUBLEs.
 *   (Bit compatability can be verified.) 
 *
 *  Created by Ian Kirker on 14/04/2008.
 *
 */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "libDefs.h"

/* File scope plan variables - used in prepareFFTs and performFFTs */
#ifdef FFT_fftw2
	oneDplanType oneDplan;
	twoDplanType twoDplan;
	#ifdef HAS_AUTO
		parallelPlanType autoPlan;
	#endif
#else
	planType oneDplan;
	planType twoDplan;
	#ifdef HAS_AUTO
		parallelPlanType autoPlan;
	#endif
#endif

void prepareFFTs(complexType *data, int decomp, int use2DFFT, int extent, int domainSize[2], MPI_Comm commColumn)
{ /* Prepares plans for the FFTs */

	#ifdef FFT_fftw3
		int n0,n1,n2,alloc,local_n0,n_start;
		/* For the FFTW versions, we use single FFT plans and repeat them many times. We
		 *  could also use the fftw_plan_many_dft version.
		 */
		if (decomp != 0)
		{ /* We only *don't* need this when we're doing an automatic parallel call */
			/*	fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany, 
                                         fftw_complex *in, const int *inembed, 
										 int istride, int idist, 
										 fftw_complex *out, const int *onembed, 
                                         int ostride, int odist, 
                                         int sign, unsigned flags); */
			oneDplan = fftw_plan_many_dft( 1, &extent, domainSize[0]*domainSize[1], 
			                               data, NULL, 1, extent, data, NULL, 1, 
										   extent, FFTW_FORWARD, FFTW_MEASURE );
		}
		
		if (use2DFFT == 1)
		{ /* We only need this when we're doing a slab decomp with the 2d FFT. */
			/* Uses same plan on many sequences */
			twoDplan = fftw_plan_dft_2d( extent, extent, data, data, FFTW_FORWARD, FFTW_MEASURE );
		}
		
		#ifdef HAS_AUTO
		if (decomp == 0)
		{
			/* 3D MPI FFT only in alpha version. */
			fftw_mpi_init();
			/* extern fftw_plan fftw_mpi_plan_dft_3d (ptrdiff_t n0, 
			      ptrdiff_t n1, ptrdiff_t n2, fftw_complex *in, fftw_complex *out, 
				  MPI_Comm comm, int sign, unsigned flags);
			 */
			autoPlan = fftw_mpi_plan_dft_3d ( extent, extent,
			                                  extent, data, data, commColumn, 
							                  FFTW_FORWARD, FFTW_MEASURE );

		}
		#endif
	#endif

	#ifdef FFT_fftw2
		if (decomp != 0)
		{ /* We only *don't* need this when we're doing an automatic parallel call */
			/* Unlike FFTW3, you need to specifiy in-place here */
			oneDplan = fftw_create_plan( extent, FFTW_FORWARD, FFTW_MEASURE | FFTW_IN_PLACE);
		}
		
		if (use2DFFT == 1)
		{ /* We only need this when we're doing a slab decomp. */
			twoDplan = fftw2d_create_plan(extent, extent, FFTW_FORWARD, FFTW_MEASURE | FFTW_IN_PLACE);
		}
		
		#ifdef HAS_AUTO
			if (decomp == 0)
			{ /* We only need this when we're doing the automatic parallel FFT */
			  /* All MPI transforms are in-place. */
				autoPlan = fftw3d_mpi_create_plan(commColumn, extent, extent, extent,
                                              FFTW_FORWARD, FFTW_MEASURE);
			}
		#endif
	#endif

	#ifdef FFT_mkl
		long status;
		long twoDdims[2] = { extent,extent };
		long autoDims[3] = { extent,extent,extent };
		
		/* 1D */
		if (decomp != 0)
		{
			status = DftiCreateDescriptor( &oneDplan, DFTI_DOUBLE, DFTI_COMPLEX, 1, extent ); 
			status = DftiSetValue( oneDplan, DFTI_NUMBER_OF_TRANSFORMS, domainSize[0]*domainSize[1] ); 
			status = DftiSetValue( oneDplan, DFTI_INPUT_DISTANCE, extent ); 
			status = DftiSetValue( oneDplan, DFTI_OUTPUT_DISTANCE, extent ); 
			status = DftiCommitDescriptor( oneDplan );
		}
		
		/* 2D */
		if (use2DFFT == 1)
		{
			status = DftiCreateDescriptor( &twoDplan, DFTI_DOUBLE, DFTI_COMPLEX, 2, twoDdims );
			status = DftiSetValue( twoDplan, DFTI_TRANSPOSE, DFTI_ALLOW );
			status = DftiCommitDescriptor( twoDplan );
		}
		
		/* 3D Auto Parallel */
		#ifdef HAS_AUTO
		if (decomp == 0)
		{
			status = DftiCreateDescriptorDM(commColumn, &autoPlan, DFTI_DOUBLE, DFTI_COMPLEX, 3, autoDims );
			status = DftiCommitDescriptorDM( autoPlan );
		}
		#endif
	#endif

	#ifdef FFT_acml
		/* In ACML, the closest you can get to making a plan is to run
		 * a transform and have it save the data in the array it uses as a buffer.
		 */
		int err;
		
		oneDplan = malloc( ( (extent * 5) + 100 ) * sizeof(complexType) );
		
		/* The C prototype isn't in the manual, so I copy it here for helpfulness. */
		/* See also page 39 of the ACML User Guide */
		/* extern void zfft1mx(int mode, double scale, int inpl, 
								int nseq, int n, doublecomplex *x, 
								int incx1, int incx2, doublecomplex *y, 
								int incy1, int incy2, doublecomplex *comm, 
								int *info);*/
		zfft1mx( 100, (double)1.0, 1, domainSize[0]*domainSize[1],
		         extent, data, 1, extent, NULL, 1, extent, oneDplan, &err);
		
		if (use2DFFT == 1)
		{
			 twoDplan = malloc( ( (extent * extent) + (extent * 10) + 200 ) * sizeof(complexType) );
			
			/*
			 * extern void zfft2dx(int mode, double scale, int ltrans, 
			 *                     int inpl, int m, int n, doublecomplex *x, 
			 *                     int incx1, int incx2, doublecomplex *y, 
			 *                     int incy1, int incy2, doublecomplex *comm, int *info);
			 */

			zfft2dx( 100, 1.0, 0, 1, extent, extent, data, 1, extent, data, 1, extent, twoDplan, &err );
		}
	#endif

	#ifdef FFT_essl
		int colSize; /* Size of task column (and thus, total number of tasks)*/
		int colPos; /* Position in column - required for PESSL's 3D FFT prep.*/
		int tempBlacsMap; /* Temporary process mapping to map MPI procs to BLACS ones */
		int *procArray; /* New mapping map. Er. old BLACS proc -> new BLACS proc */
		int tempBlacsGridSize[2]; /* We don't actually care about this, but the
		                           * function returns it and I didn't want to NULL it. */
		int tempBlacsGridCoords[2]; /* This proc's coords in tempBlacsMap. */
		
		/* essl's dcft functions requires the function to be called once for init
		 *  before you actually use it to calculate the fft.
		 *  In the first call, it saves the data to an array (aux1 in the docs)
		 *   which has to be enormous. Seriously.
		 *  Then it -also- needs a working space, of equal magnitude, but this can
		 *   be dynamically allocated by the routine if you specify 0 for its size.
		 */
		   
		/* void   esvdcft   (int, void *, int, int, void *, int, int, int,
                   int, int, double, double *, int, double *, int) */
		
	
		dcft( 1,               /* Is this a planning call only? */
		      data,            /* Pointer to data in */
		      1,               /* Stride between elements */
			  extent,          /* Stride between sequences */
			  data,            /* Output target */
			  1,               /* Stride between elements */
			  extent,          /* Stride between sequences */
			  extent,          /* Length of sequences */
			  domainSize[0]*domainSize[1], /* Number of sequences */
			  +1,              /* Forward or backward? +/- */
			  (double)1.0,     /* Scale factor for output */
			  oneDplan,        /* Planning storage */
			  sizeof(oneDplan)/sizeof(double), /* Size of planning storage */
			  NULL,          /* Working area for calculation */
			  0					/* Size of working area */
			  );
			  
		if (use2DFFT == 1)
		{
			dcft2( 1,			/* Planning call */
				data,            /* Pointer to data in */
				1,               /* Stride between elements in first dimension */
				extent,          /* Stride between elements in second dimension */
				data,            /* Output target */
				1,               /* Stride between elements in first dimension */
				extent,          /* Stride between elements in second dimension */
				extent,          /* Length in the first dimension */
				extent,			 /* Length in the second dimension */
				+1,              /* Forward or backward? +/- */
				1.0,               /* Scale factor for output */
				twoDplan,        /* Planning storage */
				sizeof(twoDplan)/sizeof(double), /* Size of planning storage */
				NULL,          /* Working area for calculation */
				0				/* Size of working area */
			  );
		}
		
		#ifdef HAS_AUTO
		if (decomp == 0)
		{
			/* See page 821 of the PESSL reference. */
			/* Also http://www.netlib.org/blacs/BLACS/Examples.html#BLACS_GRIDMAP */
			/* For PESSL, you need to use a BLACS decomposition. 
			 *  This means you need to generate a BLACS decomp that
			 *  matches our 1D MPI_Cart decomp. (PESSL 3D FFT can
			 *  only manage a slab decomp.)
			 */
			MPI_Comm_size(commColumn, &colSize);
			MPI_Comm_rank(commColumn, &colPos);
			procArray = malloc(colSize * sizeof(int));
			if (procArray == NULL) 
			{ 
				fprintf(stderr,
					"Unable to alloc procArray in routine prepareFFTs (libDefs.c)\n"
					); 
				exit(1); 
			}
			
			/* Generate a temporary mapping to get process numbers */
			Cblacs_get(0,0,&tempBlacsMap);
			
			/* Our decomps assumes a column for the slab decomp, but PESSL
			 *  wants a row. */
			Cblacs_gridinit( &tempBlacsMap, "Col", colSize, 1 );
			/* void Cblacs_gridinfo(int context, int* np_row, int* np_col, 
			                        int* my_row, int* my_col); */
			Cblacs_gridinfo( tempBlacsMap, &(tempBlacsGridSize[0]), &(tempBlacsGridSize[1]),
			                  &(tempBlacsGridCoords[0]), &(tempBlacsGridCoords[1]) );

			memset(procArray, 0, colSize*sizeof(int));
			procArray[colPos] = tempBlacsGridCoords[0];
			MPI_Allreduce(procArray, procArray, colSize, MPI_INT, MPI_MAX, commColumn);
			
			Cblacs_get(0,0,&autoPlan);
			/* void Cblacs_gridmap(int* comm, int*  usrmap, int ldu, int nprow, int npcol); */
			Cblacs_gridmap(&autoPlan, procArray, 1, 1, colSize);
			
			free(procArray);
		}
		#endif /* HAS_AUTO */
	#endif
}


void performFFTset(complexType *data, complexType *buffer, int extent, int domainSize[2])
{ /* Performs a domain's worth of 1D FFTs */
	int i;

	#ifdef FFT_fftw3
		fftw_execute( oneDplan );
	#endif

	#ifdef FFT_fftw2
		/* If left to its own devices, FFTW2 will malloc a temporary 
		 *  array each call for working storage. This is less than 
		 *  efficient, so we use the other data buffer.
		 * Note that the FFTW_IN_PLACE flag was specified above, to 
		 *  make it use the normal *out argument like this.
		 */
		/*void fftw(fftw_plan plan, int howmany,
          fftw_complex *in, int istride, int idist,
          fftw_complex *out, int ostride, int odist);*/
		fftw( oneDplan, domainSize[0]*domainSize[1],
		      data, 1, extent, buffer, 1, extent );
	#endif


	#ifdef FFT_mkl
		DftiComputeForward( oneDplan, data );
	#endif


	#ifdef FFT_acml
		int err;
		/* Page 39 of the ACML User Guide */
		/* extern void zfft1mx(int mode, double scale, int inpl, 
								int nseq, int n, doublecomplex *x, 
								int incx1, int incx2, doublecomplex *y, 
								int incy1, int incy2, doublecomplex *comm, 
								int *info);*/
		zfft1mx( -1, (double)1.0, 1, domainSize[0]*domainSize[1],
		         extent, data, 1, extent, NULL, 1, extent, oneDplan, &err);
	#endif

	#ifdef FFT_essl
		int workingSize;
		
		if ( domainSize[0]*domainSize[1]*extent*2 > 20000 )
		{ /* If the already allocated buffer is big enough, use it for the FFT */
			workingSize = domainSize[0] * domainSize[1] * extent * 2;
		} else {
			/* Otherwise setting this to 0 means that the function
			 *  allocates its own memory. */
			workingSize = 0;
		}

		dcft( 0,               /* Is this a planning call only? */
		      data,            /* Pointer to data in */
			  1,               /* Stride between elements */
			  extent,          /* Stride between sequences */
			  data,            /* Output target */
			  1,               /* Stride between elements */
			  extent,          /* Stride between sequences */
			  extent,          /* Length of sequences */
			  domainSize[0]*domainSize[1], /* Number of sequences */
			  +1,              /* Forward or backward? +/- */
			  (double)1.0,               /* Scale factor for output */
			  oneDplan,        /* Planning storage */
			  sizeof(oneDplan)/sizeof(double), /* Size of planning storage */
			  NULL, /*(double*)(void*)buffer,*/          /* Working area for calculation */
			  0 /*workingSize*/		/* Size of working area */
			  );
	#endif
}

void perform2DFFT(complexType *data, complexType *buffer, int extent, int domainSize[2])
{ /* Performs a slab domain's worth of 2D FFTs */
	int i;
	int workingSize;
	
	#ifdef FFT_fftw3
		for(i=0;i<domainSize[1];i++)
		{
			/* From the 'Guru' interface - consider using fftw_plan_many_dft instead */
			/* void fftw_execute_dft( const fftw_plan p, fftw_complex *in, fftw_complex *out); */
			fftw_execute_dft( twoDplan, data + i*extent*extent, data + i*extent*extent );
		}
	#endif

	#ifdef FFT_fftw2
		for(i=0;i<domainSize[1];i++)
		{
			/*void fftwnd_one(fftwnd_plan p, fftw_complex *in, 
				fftw_complex *out); */
			fftwnd_one(twoDplan, data + i*extent*extent, data + i*extent*extent);
		}
	#endif

	#ifdef FFT_mkl
		for(i=0;i<domainSize[1];i++)
		{
			DftiComputeForward( twoDplan, data + i*extent*extent );
		}
	#endif

	#ifdef FFT_acml
		int err;
	
		for(i=0;i<domainSize[1];i++)
		{
			/* Initial argument (mode) = -1 means forward FFT w/ precalculated plan.
			 * See prepareFFTs for full prototype */
 			zfft2dx( -1, (double)1.0, 0, 1, extent, extent, 
			        data + i*extent*extent, 1, extent, 
					data + i*extent*extent, 1, extent, twoDplan, &err );
		}
	#endif

	#ifdef FFT_essl
		if ( extent > 251 )
		{
			if ( (domainSize[0] * domainSize[1] * extent * 2) > (20000 + (2*extent + 256) * (64+2.28)) )
			{ /* If the already allocated buffer is big enough, use it for the FFT */
				workingSize = domainSize[0] * domainSize[1] * extent * 2;
			} else {
				/* Otherwise setting this to 0 means that the function allocates its own memory. */
				workingSize = 0;
			}
		} else {
			if ( domainSize[0]*domainSize[1]*extent*2 > 20000 )
			{ 
				workingSize = domainSize[0] * domainSize[1] * extent * 2;
			} else {
				workingSize = 0;
			}
		}
		
		for(i=0;i<domainSize[1];i++)
		{		
			dcft2( 0,			/* Planning call */
				data + i*extent*extent,  /* Pointer to data in */
				1,               /* Stride between elements in first dimension */
				extent,          /* Stride between elements in second dimension */
				data + i*extent*extent,  /* Output target */
				1,               /* Stride between elements in first dimension */
				extent,          /* Stride between elements in second dimension */
				extent,          /* Length in the first dimension */
				extent,			 /* Length in the second dimension */
				+1,              /* Forward or backward? +/- */
				(double)1.0,     /* Scale factor for output */
				twoDplan,        /* Planning storage */
				sizeof(twoDplan)/sizeof(double), /* Size of planning storage */
				(double*)(void*)buffer,          /* Working area for calculation */
				workingSize		/* Size of working area */
			  );
		}
	#endif
}

void performAutomatic3DFFT(complexType *data, complexType *buffer, int extent, int domainSize[2])
{ /* Uses automatic routines from a given library to perform the whole FFT */
#ifdef HAS_AUTO
	#ifdef FFT_fftw3
		fftw_execute(autoPlan);
	#endif

	#ifdef FFT_mkl
		DftiComputeForwardDM(autoPlan, data);
	#endif

	#ifdef FFT_fftw2
		/* I'm mostly trusting that my data is ordered the same way
		 *  FFTW would do it - from what I can tell from the manual,
		 *  it should be okay. */
		/* All FFTW2 MPI transforms are in-place. */
		/*void fftwnd_mpi(fftwnd_mpi_plan p,
                int n_fields,
                fftw_complex *local_data, fftw_complex *work,
                fftwnd_mpi_output_order output_order);
		*/
		fftwnd_mpi(autoPlan, 1, data, buffer, FFTW_TRANSPOSED_ORDER);
		
	#endif

	#ifdef FFT_acml
	#endif

	#ifdef FFT_essl
		int ip[40];
		ip[0]=0;
		/* pdcft3 (x, y, n1, n2, n3, isign, scale, icontxt, ip); */
		pdcft3(data, buffer, extent, extent, extent, +1, 1.0, autoPlan, ip);
		memcpy(data, buffer, extent*domainSize[0]*domainSize[1]*2*sizeof(double));
	#endif
#endif /* endif HAS_AUTO*/
}

void cleanUpFFTs(int decomp)
{ /* If applicable, free memory associated with plans. */
  /* This may not actually be necessary, but "always free what you alloc". */

	#ifdef FFT_fftw3
		if (decomp != 0)
			fftw_destroy_plan(oneDplan);
		
		if (decomp == 1)
			fftw_destroy_plan(twoDplan);
			
		#ifdef FFT_fftw3_mpi
		    if (decomp == 0)
			{
				fftw_destroy_plan(autoPlan);
				fftw_mpi_cleanup();
			}
		#endif
	#endif

	#ifdef FFT_mkl
		long status;
		if (decomp != 0)
			status = DftiFreeDescriptor( &oneDplan );
		
		if (decomp == 1)
			status = DftiFreeDescriptor( &twoDplan );
		#ifdef HAS_AUTO
			if (decomp == 0)
				status = DftiFreeDescriptorDM( &autoPlan );
		#endif
	#endif

	#ifdef FFT_fftw2
		if (decomp != 0)
			fftw_destroy_plan(oneDplan);
			
		if (decomp == 1)
			fftwnd_destroy_plan(twoDplan);
		
		#ifdef HAS_AUTO	
			if (decomp == 0)
				fftwnd_mpi_destroy_plan(autoPlan);
		#endif
	#endif

	#ifdef FFT_acml
		free(oneDplan);
		if (decomp == 1)
			free(twoDplan);
	#endif

	#ifdef FFT_essl
	#endif
}

int libraryHasAutomaticDecomposition()
{ /* Used for validateParameters. */
	#ifdef FFT_fftw3
		#ifdef HAS_AUTO
			return 1;
		#else
			return 0;
		#endif
	#endif
	
	#ifdef FFT_fftw2
		#ifdef HAS_AUTO
			return 1;
		#else
			return 0;
		#endif
	#endif	

	#ifdef FFT_acml
		return 0;
	#endif

	#ifdef FFT_mkl
		#ifdef HAS_AUTO
			return 1;
		#else
			return 0;
		#endif
	#endif
	
	#ifdef FFT_essl
		#ifdef HAS_AUTO
			return 1;
		#else
			return 0;
		#endif
	#endif	
}


/*********************************
 * Complex Number Functions.     *
 *********************************/

void complexSwap(complexType *z, complexType *w)
{ /* Swaps the values of two complex numbers. */
	complexType swap;

	#ifdef FFT_fftw3
		swap = *z;
		*z = *w;
		*w = swap;
	#endif

	#ifdef FFT_mkl
		swap = *z;
		*z = *w;
		*w = swap;
		/*swap[0] = *z[0];
		swap[1] = *z[1];
		*z[0]   = *w[0];
		*z[1]   = *w[1];
		*w[0]   = swap[0];
		*w[1]   = swap[0];*/
	#endif

	#ifdef FFT_fftw2
		swap.re = z->re;
		swap.im = z->im;
		z->re  = w->re;
		z->im  = w->im;
		w->re  = swap.re;
		w->im  = swap.im;
	#endif

	#ifdef FFT_acml
		swap.real = z->real;
		swap.imag = z->imag;
		z->real  = w->real;
		z->imag  = w->imag;
		w->real  = swap.real;
		w->imag  = swap.imag;
	#endif

	#ifdef FFT_essl
		swap._data._re = z->_data._re;
		swap._data._im = z->_data._im;
		z->_data._re = w->_data._re;
		z->_data._im = w->_data._im;
		w->_data._re = swap._data._re; 
		w->_data._im = swap._data._im;
	#endif
}

void complexSet(complexType *z, double real, double imag)
{ /* Sets a complex number specified by reference *
   *  from real and imaginary values.             */

	#ifdef FFT_fftw3
		*z = real + imag * I;
	#endif

	#ifdef FFT_mkl
		*z = real + imag * I;
	#endif

	#ifdef FFT_fftw2
		z->re = real;
		z->im = imag;
	#endif

	#ifdef FFT_acml
		z->real = real;
		z->imag = imag;
	#endif

	#ifdef FFT_essl
		z->_data._re = real;
		z->_data._im = imag;
	#endif
}

void complexAssign( complexType *z, complexType w )
{ /* Sets a complex number specified by reference *
   *  from a passed in value.                     */
	#ifdef FFT_fftw3
		*z = w;
	#endif

	#ifdef FFT_mkl
		*z = w;
	#endif

	#ifdef FFT_fftw2
		z->re = w.re;
		z->im = w.im;
	#endif

	#ifdef FFT_acml
		z->real = w.real;
		z->imag = w.imag;
	#endif

	#ifdef FFT_essl
		z->_data._re = w._data._re;
		z->_data._im = w._data._im;
	#endif
}

double complexAbsNorm(complexType z, complexType w)
{
	/* Calculates the distance between two complex numbers.     *
	 * Used in the residue calculation.                         */
	#ifdef FFT_fftw3
		return cabs(z-w);
	#endif

	#ifdef FFT_mkl
		return cabs(z-w);
	#endif

	#ifdef FFT_fftw2
		return sqrt((z.re - w.re)*(z.re - w.re) + (z.im - w.im)*(z.im - w.im));
	#endif

	#ifdef FFT_acml
		return sqrt((z.real - w.real)*(z.real - w.real) + 
		            (z.imag - w.imag)*(z.imag - w.imag));
	#endif

	#ifdef FFT_essl
		return sqrt((z._data._re - w._data._re)*(z._data._re - w._data._re) +
		            (z._data._im - w._data._im)*(z._data._im - w._data._im));
	#endif
}

_Complex double complexNative( complexType z )
{
	/* In case this isn't clear, this returns a native ISO C99       *
	 *  complex number with value equal to that of the FFT library's * 
	 *  type of complex number passed in.                            *
	 * Used for testing.                                             *
	 * The prototype would be of type complex double, but acml       *
	 *  redefines complex. Helpful.                                  */

	#ifdef FFT_fftw3
		return z;
	#endif

	#ifdef FFT_mkl
		return z;
	#endif

	#ifdef FFT_fftw2
		return z.re + z.im * I;
	#endif

	#ifdef FFT_acml
		return z.real + z.imag * I;
	#endif

	#ifdef FFT_essl
		return z._data._re + z._data._im * I;
	#endif
}

void printLib()
{
	fprintf(stderr, "This executable uses the %s library.\n", FFT_NAME);
}