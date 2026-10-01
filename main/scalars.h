/*
 * File: scalars.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * ----------------------------------------
 * Header file for scalars.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _scalars_h
#define _scalars_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

void UpdateScalars(gridT *grid, physT *phys, propT *prop, REAL **w_im, REAL **SfH_tm1, REAL **SfH_tm2, REAL **SfHv_t, REAL **SfHv_tm1,  REAL **scal, 
		   REAL **scal_old, REAL **boundary_scal, int BCtop, int BCbot, REAL **Cn, 
		   REAL kappa, REAL kappaH, REAL **kappa_tv, REAL theta,
                   REAL **src1, REAL **src2, REAL *Ftop, REAL *Fbot, int alpha_top, int alpha_bot,
                   MPI_Comm comm, int myproc, int checkflag, int timestepping, int TVDscheme);

#endif
