/*
 * File: check.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for check.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _check_h
#define _check_h

#include "grid.h"
#include "phys.h"
#include "mympi.h"

int Check(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
int CheckDZ(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
int CheckDivergence(REAL **dzz, REAL **dzz_old, REAL **u_n_plus_1, REAL **u_n, REAL **u_n_minus_1,
		    REAL **w_im, gridT *grid, propT *prop, REAL tolerance, int myproc);
void Progress(propT *prop, int myproc, int numprocs);
void MemoryStats(gridT *grid, int myproc, int numprocs, MPI_Comm comm);
void CheckGrid(gridT *grid, int myproc, MPI_Comm comm);

#endif
