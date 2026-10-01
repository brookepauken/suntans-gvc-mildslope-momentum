/*
 * File: pressure.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for pressure.c
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _pressure_h
#define _pressure_h

// Public functions
void BiCGSolveQ(REAL **q, REAL **src, REAL **c, gridT *grid, physT *phys,
		       propT *prop, int myproc, int numprocs, MPI_Comm comm);
void CGSolveQ(REAL **q, REAL **src, REAL **c, gridT *grid, physT *phys, 
		     propT *prop, int myproc, int numprocs, MPI_Comm comm);
void Corrector(REAL **qc, gridT *grid, physT *phys, propT *prop, int myproc, 
		      int numprocs, MPI_Comm comm);
void ComputeQSource(REAL **src, gridT *grid, physT *phys, propT *prop, 
			   int myproc, int numprocs);

#endif
