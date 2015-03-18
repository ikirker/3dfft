/*
 *  miscMPI.c
 * Contains simple MPI calls unrelated to the all-to-all transpositions.
 *
 *  Created by Ian Kirker on 15/05/2008.
 *
 */

#include "comms.h"
#include "libDefs.h"
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>


int amMaster(MPI_Comm comm)
{
	int rank;
	MPI_Comm_rank(comm, &rank);
	if (rank == 0) 
	{
		return 1;
	} else {
		return 0;
	}
}

int getSize(MPI_Comm comm)
{
	int size;
	MPI_Comm_size(comm, &size);
	return size;
}

int commsInit(int *argc, char ***argv)
{
	MPI_Init(argc, argv);
	return 0;
}

void commSync(MPI_Comm comm)
{
	MPI_Barrier(comm);
}

void doubleGlobalSum( double *amount, MPI_Comm comm)
{
        double amount_copy = *amount;
        MPI_Allreduce( &amount_copy, amount, 1, MPI_DOUBLE, MPI_SUM, comm );
}

void commsEnd()
{
	MPI_Finalize();
}
