/*
 * File: tvd.h
 * Author: Zhonghua Zhang and Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for tvd.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _tvd_h
#define _tvd_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"
#include "fileio.h"

// TVD Method.  0=No TVD (Uses method in original scheme, which is First-order upwind)
// Otherwise, TVD is implemented, with TVDMACRO=1: First-order upwind, 2: Lax-Wendroff, 3: Superbee, 4: Van Leer
// Note that this TVD implementation does not work with wetting and drying and is 
// not strictly tvd.
// This is now defined in defaults.h and can be set in suntans.dat for each of salt, temperature, and turbulence
//#define TVDMACRO 4

void HorizontalFaceScalars(gridT *grid, physT *phys, propT *prop, REAL **scal, REAL **U, REAL **boundary_scal, int TVD,
			   int dimen, MPI_Comm comm, int myproc);
void HorizontalFaceScalars_oldstep(gridT *grid, physT *phys, propT *prop, REAL **scal, REAL **faceval, REAL **boundary_scal, int TVD, 
				int timestep, int dimen, MPI_Comm comm, int myproc);
void GetApAm(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
	     REAL **w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot);
void GetApAmVert(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
	REAL *w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot);
void GetApAmNewAdv(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
	REAL *w, REAL **dzz, REAL **scal, int j, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot);
void GetApAm_linear(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
	     REAL **w, REAL *dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD);
void HorizontalFaceU(REAL **uc, gridT *grid, physT *phys, propT *prop, int TVDscheme,
			   MPI_Comm comm, int myproc); 
void GetPentDiag(REAL *pent_a, REAL *pent_b, REAL *pent_c, REAL *pent_d, REAL *pent_e, REAL *wp, REAL *wm,
		 REAL **w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot);
void GetPentDiagNewAdv(REAL *pent_a, REAL *pent_b, REAL *pent_c, REAL *pent_d, REAL *pent_e, REAL *wp, REAL *wm,
	     REAL *w, REAL **dzz, REAL **U, int j, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot);
void GetPentDiagVert(REAL *pent_a, REAL *pent_b, REAL *pent_c, REAL *pent_d, REAL *pent_e, REAL *wp, REAL *wm,
		 REAL *w, REAL **dzz, REAL **W, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot);
void GetScalarOnFace(REAL **SfHvf, REAL *wp, REAL *wm, REAL **w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot, REAL **sum_neighs);
void SumNeighborScalars(gridT *grid, physT *phys, REAL **scal, REAL **sum_neighs,
	MPI_Comm comm, int myproc);
#endif


