/*
 *  libDefs.h
 *  This file and its associated source file contain any functions and
 *   statements (within the actual program) that need to be modified 
 *   when adding a new library. This includes any work done on complex 
 *   numbers, FFT planning and execution, and would also contain the 
 *   definition of the appropriate MPI type to represent the complex 
 *   number, except that I decided to assume they were all just two 
 *   MPI_DOUBLEs. (Bit compatability can be verified.) 
 *
 * Defining HAS_AUTO in compilation is an instruction to use the parallel
 *  version of the library as well.
 *
 * We use the less eye-friendly version of the C complex number type
 *  to avoid namespace collisions -- complex is equivalent to _Complex.
 *
 *  Created by Ian Kirker on 14/04/2008.
 *
 */

#ifndef HEADER_LIBDEFS
#define HEADER_LIBDEFS

#include <complex.h>
#include <mpi.h>

#ifndef FFT_fftw2
#ifndef FFT_fftw3
#ifndef FFT_mkl
#ifndef FFT_essl
#ifndef FFT_acml
#define FFT_fftw3 /* Default to make test building easier. */
#endif
#endif
#endif
#endif
#endif

/* Include FFT library of choice */
#ifdef FFT_fftw2
	#ifdef FFTW_TYPE_SPECIFIED
		#include <dfftw.h>
		#ifdef HAS_AUTO
			#include <dfftw_mpi.h>
		#endif
	#else
		#include <fftw.h>
		#ifdef HAS_AUTO
			#include <fftw_mpi.h>
		#endif
	#endif
	
	#define FFT_NAME "fftw2"
	#define FFT_fftw2_LIBKEY 
	typedef fftw_complex complexType;
	typedef fftw_plan oneDplanType;
	typedef fftwnd_plan twoDplanType;
	#ifdef HAS_AUTO
		typedef fftwnd_mpi_plan parallelPlanType;
	#endif
#endif

#ifdef FFT_fftw3
	#include <fftw3.h>
	#define FFT_NAME "fftw3"
	#define FFT_fftw3_LIBKEY 1
	
	typedef fftw_complex complexType;
	typedef fftw_plan planType;
	#ifdef HAS_AUTO
		typedef fftw_plan parallelPlanType;
	#endif
#endif

#ifdef FFT_essl
	#include <essl.h>
	#define FFT_NAME "essl"
	#define FFT_essl_LIBKEY 2
	
	typedef dcmplx complexType;
	typedef double planType[20000];
	
	#ifdef HAS_AUTO
		#include <pessl.h>
		#include <Cblacs.h>
		/* This refers to the BLACS handle you get back from the grid init calls. */
		typedef int parallelPlanType; 
	#endif
#endif

#ifdef FFT_acml
	#include <acml.h>
	#define FFT_NAME "acml"
	#define FFT_acml_LIBKEY 3
	typedef doublecomplex complexType;
	
	/* In ACML, the closest you can get to making a plan is to run
	 * a transform and have it save the data in the array it uses as a buffer,
	 * so we alloc some space later but say that the planType is a pointer now.
	 */
	typedef complexType *planType;
#endif

#ifdef FFT_mkl
	#include <mkl_dfti.h>
	#define FFT_NAME "mkl"
	#define FFT_mkl_LIBKEY 4
	typedef _Complex double complexType;
	typedef DFTI_DESCRIPTOR_HANDLE planType;
	#ifdef HAS_AUTO
		typedef DFTI_DESCRIPTOR_DM_HANDLE parallelPlanType;
	#endif
#endif

void prepareFFTs(complexType *data, int decomp, int use2DFFT, int extent, int domainSize[2], MPI_Comm commColumn);
void performFFTset(complexType *data, complexType *buffer, int extent, int domainSize[2]);
void perform2DFFT(complexType *data, complexType *buffer, int extent, int domainSize[2]);
void performAutomatic3DFFT(complexType *data, complexType *buffer, int extent, int domainSize[2]);
void cleanUpFFTs(int decomp);
int libraryHasAutomaticDecomposition();

void printLib();
void complexSwap(complexType *, complexType *);
void complexSet(complexType *z, double real, double imag);
void complexAssign( complexType *z, complexType w );
double complexAbsNorm(complexType z, complexType w);
_Complex double complexNative( complexType z );

#endif