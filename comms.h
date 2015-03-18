/*
 *  comms.h
 *  Communications functions not related to the decomposition.
 *
 *  Created by Ian Kirker on 15/05/2008.
 *
 */

#ifndef HEADER_COMMS 

#include <mpi.h>

int amMaster(MPI_Comm comm);
int commsInit(int *argc, char ***argv);
int getSize(MPI_Comm comm);
void commSync(MPI_Comm comm);
void doubleGlobalSum( double *amount, MPI_Comm comm);
void commsEnd();

#define HEADER_COMMS
#endif