/*
 *  options.h
 *  Prototypes for options.c
 *  Gets the options from the arguments.
 *  Created by Ian Kirker on 15/05/2008.
 * 
 *
 */

int getOptions(int *argc, char ***argv, int *extent, int *decompType, int *use2DFFT, int *skip, int *skipFFT, int *targetLoopCount, int *printOut);
void printOptionList();
