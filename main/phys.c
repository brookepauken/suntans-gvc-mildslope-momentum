/*
 * File: phys.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 *  --------------------------------
 * This file contains physically-based functions.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "util.h"
#include "initialization.h"
#include "memory.h"
#include "turbulence.h"
#include "boundaries.h"
#include "check.h"
#include "scalars.h"
#include "timer.h"
#include "profiles.h"
#include "state.h"
#include "diffusion.h"
#include "sources.h"
#include "mynetcdf.h"
#include "met.h"
#include "age.h"
#include "physio.h"
#include "merge.h"
#include "sediments.h"
#include "marsh.h"
#include "vertcoordinate.h"
#include "culvert.h"
#include "wave.h"
#include "subgrid.h"
#include "sendrecv.h"
#include "upredictor.h"
#include "pressure.h"

const int DEBUG=0;

//Cell-centered+interpolation scheme: 0
//Face-centered conservative scheme: 1
const int U_MOMENTUM_ADVECTION=1, W_MOMENTUM_ADVECTION=1;

/*
 * Private Function declarations.
 *
 */
static void GSSolve(gridT *grid, physT *phys, propT *prop, 
		    int myproc, int numprocs, MPI_Comm comm);
static void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop);
void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop);
static void EddyViscosity(gridT *grid, physT *phys, propT *prop, REAL **wnew, 
			  MPI_Comm comm, int myproc);
static void HorizontalSourceTerms_n1(gridT *grid, physT *phys, propT *prop,
          int myproc, int numprocs, MPI_Comm comm);
static void HorizontalSource(gridT *grid, physT *phys, propT *prop,
			     int myproc, int numprocs, MPI_Comm comm);
static void StoreVariables(gridT *grid, physT *phys, propT *prop);
static void NewCells(gridT *grid, physT *phys, propT *prop);
static void WPredictor(gridT *grid, physT *phys, propT *prop,
		       int myproc, int numprocs, MPI_Comm comm);
static void WSourceTerms(gridT *grid, physT *phys, propT *prop,
           int myproc, int numprocs, MPI_Comm comm);
void ComputeUC(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc,
	       interpolation interp, int kinterp, int subgridmodel);
void ComputeUCPerot(REAL **u, REAL **uc, REAL **vc, REAL *h, int kinterp,
			   int subgridmodel, gridT *grid);
static void ComputeUCLSQ(REAL **u, REAL **uc, REAL **vc, gridT *grid, physT *phys);
static void ComputeUCRT(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc);
static void ComputeNodalVelocity(physT *phys, gridT *grid, interpolation interp, int myproc);
static void  ComputeTangentialVelocity(physT *phys, gridT *grid, interpolation ninterp, interpolation tinterp,int myproc);
static void  ComputeQuadraticInterp(REAL x, REAL y, int ic, int ik, REAL **uc, 
    REAL **vc, physT *phys, gridT *grid, interpolation ninterp, 
    interpolation tinterp, int myproc);
static void ComputeRT0Velocity(REAL* tempu, REAL* tempv, REAL e1n1, REAL e1n2, 
    REAL e2n1, REAL e2n2, REAL Uj1, REAL Uj2);
static void BarycentricCoordsFromCartesian(gridT *grid, int cell, 
    REAL x, REAL y, REAL* lambda);
static void BarycentricCoordsFromCartesianEdge(gridT *grid, int cell, 
    REAL x, REAL y, REAL* lambda);
static REAL UFaceFlux(int j, int k, REAL **phi, REAL **u, gridT *grid, REAL dt, 
    int method);
static REAL HFaceFlux(int j, int k, REAL *phi, REAL **u, gridT *grid, REAL dt, 
    int method);
//static void SetFluxHeight(gridT *grid, physT *phys, propT *prop);
//void SetFluxHeight(gridT *grid, physT *phys, propT *prop, int dzfmeth, MPI_Comm comm, int myproc);
static void GetMomentumFaceValues(REAL **uface, REAL **ui, REAL **boundary_ui, REAL **U, gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc, int nonlinear, int TVD);
static void getTsurf(gridT *grid, physT *phys);
static void getchangeT(gridT *grid, physT *phys);
// added drag function Yun Zhang
static void InterpDrag(gridT *grid, physT *phys, propT *prop, int myproc);
static void OutputDrag(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static int UpdateBottomHeight(REAL *zB, REAL *zBold, REAL *zBold2, REAL *zBoffline, gridT *grid,
			      propT *prop, physT *phys, int myproc, MPI_Comm comm);
/*
 * Function: AllocatePhysicalVariables
 * Usage: AllocatePhysicalVariables(grid,phys,prop);
 * -------------------------------------------------
 * This function allocates space for the physical arrays but does not
 * allocate space for the grid as this has already been allocated.
 *
 */
void AllocatePhysicalVariables(gridT *grid, physT **phys, propT *prop)
{
  int flag=0, i, j, jptr, ib, Nc=grid->Nc, Ne=grid->Ne, Np=grid->Np, nf, k;

  // allocate physical structure
  *phys = (physT *)SunMalloc(sizeof(physT),"AllocatePhysicalVariables");

  // allocate  variables in plan
  (*phys)->u = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->uc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->vc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");

  // new variables for higher-order interpolation following Wang et al 2011
  (*phys)->nRT1u = (REAL ***)SunMalloc(Np*sizeof(REAL **),"AllocatePhysicalVariables");
  (*phys)->nRT1v = (REAL ***)SunMalloc(Np*sizeof(REAL **),"AllocatePhysicalVariables");
  (*phys)->nRT2u = (REAL **)SunMalloc(Np*sizeof(REAL*),"AllocatePhysicalVariables");
  (*phys)->nRT2v = (REAL **)SunMalloc(Np*sizeof(REAL*),"AllocatePhysicalVariables");
  (*phys)->tRT1 = (REAL **)SunMalloc(Ne*sizeof(REAL*),"AllocatePhysicalVariables");
  (*phys)->tRT2 = (REAL **)SunMalloc(Ne*sizeof(REAL*),"AllocatePhysicalVariables");

  // allocate rest of variables in plan
  (*phys)->uold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->vold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->uold2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->vold2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->D = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->utmp = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->u_old = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->u_old2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->u_ex = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");  
  (*phys)->ut = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_U = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_U2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables"); //AB3


  // for each variable in plan consider the number of layers it affects
  for(j=0;j<Ne;j++) {
    // the following line seems somewhat dubious...
    if(grid->Nkc[j] < grid->Nke[j]) {
      printf("Error!  Nkc(=%d)<Nke(=%d) at edge %d\n",grid->Nkc[j],grid->Nke[j],j);
      flag = 1;
    }
    // allocate memory for the max (cell-centered) quanity on the edge (from definition above)
    (*phys)->u[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->utmp[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->u_old[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->u_old2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->u_ex[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");    
    (*phys)->ut[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_U[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_U2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");//AB3
    /* new interpolation variables */
    // loop over the edges (Nkc vs Nke since for cells Nkc < ik < Nke there should be 0 velocity
    // on face to prevent mass from leaving the system)
    (*phys)->tRT1[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->tRT2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
  }
  // if we have an error quit MPI
  if(flag) {
    MPI_Finalize();
    exit(0);
  }

  // user defined variable
  (*phys)->user_def_nc = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->user_def_nc_nk = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  for(i=0;i<Nc;i++)
    (*phys)->user_def_nc_nk[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");

  // cell-centered physical variables in plan (no vertical direction)
  (*phys)->h = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->dhdt = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->zB = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->zBold = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->zBold2 = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");      
  (*phys)->zBoffline = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->Erosion = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");      
  (*phys)->Erosion_old = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");      
  (*phys)->Deposition = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");      
  (*phys)->hcorr = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->active = (unsigned char *)SunMalloc(Nc*sizeof(char),"AllocatePhysicalVariables");
  (*phys)->hold = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->h_old = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->htmp = (REAL *)SunMalloc(10*Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->htmp2 = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->htmp3 = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->hcoef = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->hfcoef = (REAL *)SunMalloc(grid->maxfaces*Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->Tsurf = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->dT = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocatePhysicalVariables");

  // cell-centered values that are also depth-varying
  (*phys)->w = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wnew = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->wtmp = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->w_old = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->w_old2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->w_ex = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");  
  (*phys)->w_im = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->w_s = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->w_sed = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->S_SfHv_t = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->S_SfHv_tm1 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->T_SfHv_t = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->T_SfHv_tm1 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->u_SfHv_t = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->u_SfHv_tm1 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->w_SfHv_t = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->w_SfHv_tm1 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->v_SfHv_t = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->v_SfHv_tm1 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_W = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_W2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables"); //AB3
  (*phys)->q = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->qc = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->qtmp = (REAL **)SunMalloc(grid->maxfaces*Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->s = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->T = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->s_old = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->T_old = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Ttmp = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->s0 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->rho = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_R = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Cn_T = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp3 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->stmp4 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->nu_tv = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->kappa_tv = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->SediKappa_tv = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->nu_lax = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  if(prop->turbmodel>=1) {
    (*phys)->qT = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->lT = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->qT_old = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->lT_old = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->Cn_q = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
    (*phys)->Cn_l = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  }
  (*phys)->tau_T = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->tau_B = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdT = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->CdB = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");

  //add phys->z0T and phys->z0B
  (*phys)->z0B = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->z0T = (REAL *)SunMalloc(Ne*sizeof(REAL),"AllocatePhysicalVariables");

  /* new interpolation variables */
  // loop over the nodes
  for(i=0; i < Np; i++) {
    // most complex one...
    (*phys)->nRT1u[i] = (REAL **)SunMalloc(grid->Nkp[i]*sizeof(REAL *),
        "AllocatePhysicalVariables");
    (*phys)->nRT1v[i] = (REAL **)SunMalloc(grid->Nkp[i]*sizeof(REAL *),
        "AllocatePhysicalVariables");
    // loop over all the cell over depth and allocate (note that for rapidly varying
    // bathymetery this will result in too much memory used in some places)
    for(k=0; k < grid->Nkp[i]; k++) {
      (*phys)->nRT1u[i][k] = (REAL *)SunMalloc(grid->numpcneighs[i]*sizeof(REAL),
          "AllocatePhysicalVariables");
      (*phys)->nRT1v[i][k] = (REAL *)SunMalloc(grid->numpcneighs[i]*sizeof(REAL),
          "AllocatePhysicalVariables");
    }
    // simpler one
    (*phys)->nRT2u[i] = (REAL *)SunMalloc(grid->Nkp[i]*sizeof(REAL),
        "AllocatePhysicalVariables");
    (*phys)->nRT2v[i] = (REAL *)SunMalloc(grid->Nkp[i]*sizeof(REAL),
        "AllocatePhysicalVariables");
  }
  // Netcdf write variables
  (*phys)->tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->tmpvarE = (REAL *)SunMalloc(grid->Ne*grid->Nkmax*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->tmpvarW = (REAL *)SunMalloc(grid->Nc*(grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
 
  // for each cell allocate memory for the number of layers at that location
  for(i=0;i<Nc;i++) {
    (*phys)->uc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->vc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->uold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->vold[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->uold2[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->vold2[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");    
    (*phys)->w[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wnew[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->wtmp[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w_old[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w_old2[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w_ex[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");    
    (*phys)->w_im[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w_s[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w_sed[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->S_SfHv_t[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->S_SfHv_tm1[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->T_SfHv_t[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->T_SfHv_tm1[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->u_SfHv_t[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->u_SfHv_tm1[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w_SfHv_t[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w_SfHv_tm1[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->v_SfHv_t[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->v_SfHv_tm1[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_W[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_W2[i] = (REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocatePhysicalVariables"); //AB3
    (*phys)->q[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->qc[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    for(nf=0;nf<grid->nfaces[i];nf++)
      (*phys)->qtmp[i*grid->maxfaces+nf] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->s[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->T[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->s_old[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->T_old[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Ttmp[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->s0[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->rho[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_R[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->Cn_T[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    if(prop->turbmodel>=1) {
      (*phys)->Cn_q[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
      (*phys)->Cn_l[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
      (*phys)->qT[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
      (*phys)->lT[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
      (*phys)->qT_old[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
      (*phys)->lT_old[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    }
    (*phys)->stmp[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->stmp2[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->stmp3[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->stmp4[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");    
    (*phys)->nu_tv[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->kappa_tv[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->SediKappa_tv[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->nu_lax[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
  }
 
  // allocate boundary value memory
  (*phys)->boundary_u = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_v = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_w = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_s = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_T = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_rho = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_tmp = (REAL **)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->boundary_h = (REAL *)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->boundary_flag = (REAL *)SunMalloc((grid->edgedist[5]-grid->edgedist[2])*sizeof(REAL),"AllocatePhysicalVariables");
  // allocate over vertical layers
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
    j=grid->edgep[jptr];

    (*phys)->boundary_u[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_v[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_w[jptr-grid->edgedist[2]] = (REAL *)SunMalloc((grid->Nke[j]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_s[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_T[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_tmp[jptr-grid->edgedist[2]] = (REAL *)SunMalloc((grid->Nke[j]+1)*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->boundary_rho[jptr-grid->edgedist[2]] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"AllocatePhysicalVariables");
    }

  // allocate coefficients
  (*phys)->ap = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->am = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->bp = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->bm = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->a = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->b = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->c = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->d = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->e = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");  

  //pent diagonal for now
  (*phys)->pent_a = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");  
  (*phys)->pent_b = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");  
  (*phys)->pent_c = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");  
  (*phys)->pent_d = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");  
  (*phys)->pent_e = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");  
  (*phys)->pent_source = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");  

  // Allocate for the face scalar
  (*phys)->SfHp = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->SfHm = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  
  //salinity and temp at previous timesteps
  (*phys)->S_SfH_tm1 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->S_SfH_tm2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->T_SfH_tm1 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->T_SfH_tm2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
 
  //same but horizontal and vertical velocity 
  (*phys)->u_SfH_tm1 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->u_SfH_tm2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->w_SfH_tm1 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->w_SfH_tm2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->v_SfH_tm1 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->v_SfH_tm2 = (REAL **)SunMalloc(Ne*sizeof(REAL *),"AllocatePhysicalVariables");

  for(j=0;j<Ne;j++) {
    (*phys)->SfHp[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->SfHm[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");

    (*phys)->S_SfH_tm1[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->S_SfH_tm2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->T_SfH_tm1[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->T_SfH_tm2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
  
    (*phys)->u_SfH_tm1[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->u_SfH_tm2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w_SfH_tm1[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->w_SfH_tm2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->v_SfH_tm1[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->v_SfH_tm2[j] = (REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocatePhysicalVariables");
  }

  // Allocate for TVD schemes
  (*phys)->Cp = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->Cm = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->rp = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->rm = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");

  (*phys)->wp = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");
  (*phys)->wm = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocatePhysicalVariables");

  (*phys)->gradSx = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->gradSy = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->sum_neighs = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocatePhysicalVariables");
  for(i=0;i<Nc;i++) {
    (*phys)->gradSx[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->gradSy[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
    (*phys)->sum_neighs[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocatePhysicalVariables");
  }


  // Allocate for least squares velocity fitting
  (*phys)->A = (REAL **)SunMalloc(grid->maxfaces*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->AT = (REAL **)SunMalloc(2*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->Apr = (REAL **)SunMalloc(2*sizeof(REAL *),"AllocatePhysicalVariables");
  (*phys)->bpr = (REAL *)SunMalloc(2*sizeof(REAL),"AllocatePhysicalVariables");
  for(i=0;i<grid->maxfaces;i++){
      (*phys)->A[i] = (REAL *)SunMalloc(2*sizeof(REAL),"AllocatePhysicalVariables");  
  }
  for(i=0;i<2;i++){
      (*phys)->AT[i] = (REAL *)SunMalloc(grid->maxfaces*sizeof(REAL),"AllocatePhysicalVariables");  
      (*phys)->Apr[i] = (REAL *)SunMalloc(2*sizeof(REAL),"AllocatePhysicalVariables");  
  }

}

/*
 * Function: FreePhysicalVariables
 * Usage: FreePhysicalVariables(grid,phys,prop);
 * ---------------------------------------------
 * This function frees all space allocated in AllocatePhysicalVariables
 *
 */
void FreePhysicalVariables(gridT *grid, physT *phys, propT *prop)
{
  int i, j, Nc=grid->Nc, Ne=grid->Ne, Np=grid->Np, nf;

  /* free variables for higher-order interpolation */
  // note that this isn't even currently called!
  // over each node
  for(i=0; i < Np; i++) {
    free(phys->nRT2u[i]);
    free(phys->nRT2v[i]);
    // over each layer
    for(j=0; j < grid->Nkp[i]; j++) {
      // free the memory
      free(phys->nRT1u[i][j]);
      free(phys->nRT1v[i][j]);
    }
  }
  // over each edge
  for(i=0; i < Ne; i++) {
    free(phys->tRT1[i]);
    free(phys->tRT2[i]);
  }

  // free all the arrays over depth for edge-oriented
  for(j=0;j<Ne;j++) {
    free(phys->u[j]);
    free(phys->utmp[j]);
    free(phys->u_old[j]);
    free(phys->u_old2[j]);
    free(phys->u_ex[j]);
    free(phys->ut[j]);
    free(phys->Cn_U[j]);
    free(phys->Cn_U2[j]); //AB3
  }

  // free all the arrays over depth for cell-oriented
  for(i=0;i<Nc;i++) {
    free(phys->uc[i]);
    free(phys->vc[i]);
    free(phys->wc[i]);
    free(phys->uold[i]);
    free(phys->vold[i]);
    free(phys->uold2[i]);
    free(phys->vold2[i]);
    free(phys->w[i]);
    free(phys->wnew[i]);
    free(phys->wtmp[i]);
    free(phys->w_old[i]);
    free(phys->w_old2[i]);
    free(phys->w_ex[i]);    
    free(phys->w_im[i]);
    free(phys->w_s[i]);
    free(phys->w_sed[i]);
    free(phys->S_SfHv_t[i]);
    free(phys->S_SfHv_tm1[i]);
    free(phys->T_SfHv_t[i]);
    free(phys->T_SfHv_tm1[i]);
    free(phys->u_SfHv_t[i]);
    free(phys->u_SfHv_tm1[i]);
    free(phys->w_SfHv_t[i]);
    free(phys->w_SfHv_tm1[i]);
    free(phys->v_SfHv_t[i]);
    free(phys->v_SfHv_tm1[i]);
    free(phys->Cn_W[i]);
    free(phys->Cn_W2[i]); //AB3
    free(phys->q[i]);
    free(phys->qc[i]);
    for(nf=0;nf<grid->maxfaces;nf++)
      free(phys->qtmp[i*grid->maxfaces+nf]);
    free(phys->s[i]);
    free(phys->T[i]);
    free(phys->s_old[i]);
    free(phys->T_old[i]);
    free(phys->s0[i]);
    free(phys->rho[i]);
    free(phys->Cn_R[i]);
    free(phys->Cn_T[i]);
    if(prop->turbmodel>=1) {
      free(phys->Cn_q[i]);
      free(phys->Cn_l[i]);
      free(phys->qT[i]);
      free(phys->lT[i]);
      free(phys->qT_old[i]);
      free(phys->lT_old[i]);
    }
    free(phys->stmp[i]);
    free(phys->stmp2[i]);
    free(phys->stmp3[i]);
    free(phys->stmp4[i]);    
    free(phys->nu_tv[i]);
    free(phys->kappa_tv[i]);
    free(phys->SediKappa_tv[i]);
    free(phys->nu_lax[i]);
    free(phys->user_def_nc_nk[i]);
  }
  free(phys->user_def_nc);
  free(phys->user_def_nc_nk);
  free(phys->h);
  free(phys->hcorr);
  free(phys->htmp);
  free(phys->htmp2);
  free(phys->htmp3);
  free(phys->h_old);
  free(phys->hold);
  free(phys->hcoef);
  free(phys->hfcoef);
  free(phys->zB);
  free(phys->zBold);
  free(phys->zBold2);
  free(phys->zBoffline);
  free(phys->uc);
  free(phys->vc);
  free(phys->wc);
  free(phys->w);
  free(phys->wnew);
  free(phys->wtmp);
  free(phys->w_old);
  free(phys->w_old2);
  free(phys->w_ex);  
  free(phys->w_im);
  free(phys->w_sed);
  free(phys->S_SfHv_t);
  free(phys->S_SfHv_tm1);
  free(phys->T_SfHv_t);
  free(phys->T_SfHv_tm1);
  free(phys->u_SfHv_t);
  free(phys->u_SfHv_tm1);
  free(phys->w_SfHv_t);
  free(phys->w_SfHv_tm1);
  free(phys->v_SfHv_t);
  free(phys->v_SfHv_tm1);
  free(phys->Cn_W);
  free(phys->Cn_W2); //AB3
  free(phys->q);
  free(phys->qc);
  free(phys->qtmp);
  free(phys->s);
  free(phys->T);
  free(phys->s_old);
  free(phys->T_old);
  free(phys->s0);
  free(phys->rho);
  free(phys->Cn_R);
  free(phys->Cn_T);
  if(prop->turbmodel>=1) {
    free(phys->Cn_q);
    free(phys->Cn_l);
    free(phys->qT);
    free(phys->lT);
    free(phys->qT_old);
    free(phys->lT_old);
  }  
  free(phys->stmp);
  free(phys->stmp2);
  free(phys->stmp3);
  free(phys->stmp4);  
  free(phys->nu_tv);
  free(phys->kappa_tv);
  free(phys->SediKappa_tv);
  free(phys->nu_lax);
  free(phys->tau_T);
  free(phys->tau_B);
  free(phys->CdT);
  free(phys->CdB);
  free(phys->z0T);
  free(phys->z0B);
  free(phys->u);
  free(phys->D);
  free(phys->utmp);
  free(phys->u_old);
  free(phys->u_old2);
  free(phys->u_ex);  
  free(phys->ut);
  free(phys->Cn_U);
  free(phys->Cn_U2);
  free(phys->ap);
  free(phys->am);
  free(phys->bp);
  free(phys->bm);
  free(phys->a);
  free(phys->b);
  free(phys->c);
  free(phys->d);

  free(phys->pent_a);
  free(phys->pent_b);
  free(phys->pent_c);
  free(phys->pent_d);
  free(phys->pent_e);
  free(phys->pent_source);

  // Free the horizontal facial scalar  
  for(j=0;j<Ne;j++) {
    free( phys->SfHp[j] );
    free( phys->SfHm[j] );

    free( phys->S_SfH_tm1[j] );
    free( phys->S_SfH_tm2[j] );
    free( phys->T_SfH_tm1[j] );
    free( phys->T_SfH_tm2[j] );

    free( phys->u_SfH_tm1[j] );
    free( phys->u_SfH_tm2[j] );
    free( phys->w_SfH_tm1[j] );
    free( phys->w_SfH_tm2[j] );
    free( phys->v_SfH_tm1[j] );
    free( phys->v_SfH_tm2[j] );
  }
  free(phys->SfHp);
  free(phys->SfHm);

  free( phys->S_SfH_tm1);
  free( phys->S_SfH_tm2);
  free( phys->T_SfH_tm1);
  free( phys->T_SfH_tm2);

  free( phys->u_SfH_tm1);
  free( phys->u_SfH_tm2);
  free( phys->w_SfH_tm1);
  free( phys->w_SfH_tm2);
  free( phys->v_SfH_tm1);
  free( phys->v_SfH_tm2);

  // Free the variables for TVD scheme
  free(phys->Cp);
  free(phys->Cm);
  free(phys->rp);
  free(phys->rm);
  free(phys->wp);
  free(phys->wm);

  free(phys->gradSx);
  free(phys->gradSy);

  free(phys->sum_neighs);

  free(phys);
}

/*
 * Function: InitializePhyiscalVariables
 * Usage: InitializePhyiscalVariables(grid,phys,prop,myproc,comm);
 * ---------------------------------------------------------------
 * This function initializes the physical variables by calling
 * the routines defined in the file initialize.c
 *
 */
void InitializePhysicalVariables(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm)
{
  int i, j, jptr, k, ktop, Nc=grid->Nc,nc1,nc2;
  REAL z, zu, *stmp, alpha_f, dmax, a, L0, L, zeta, znew;
  REAL *ncscratch;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];
  int Nci, Nki, T0;
  FILE *fid;

  prop->nstart=0;
  prop->n=prop->nstart;
  // Initialise the netcdf time
  prop->toffSet = getToffSet(prop->basetime,prop->starttime);
  prop->nctime = prop->toffSet*86400.0 + prop->nstart*prop->dt;

  if (prop->readinitialnc>0){
    ReadInitialNCcoord(prop,grid,&Nci,&Nki,&T0,myproc);

    printf("myproc: %d, Nci: %d, Nki: %d, T0: %d\n",myproc,Nci,Nki,T0);

    // Initialise a scratch variable for reading arrays
    ncscratch = (REAL *)SunMalloc(Nki*Nci*sizeof(REAL),"InitializePhysicalVariables");

  }

  if(DEBUG) printf("Updating dzz.\n");
  // Need to update the vertical grid and fix any cells in which
  // dzz is too small when h=0.
  if(prop->vertcoord==1 || prop->vertcoord==5)
   UpdateDZ(grid,phys,prop, -1);
  if(DEBUG) printf("Done with dzz\n");

  if(DEBUG) printf("Initializing free surface\n");
  // Initialize the free surface
  if (prop->readinitialnc){
     ReturnFreeSurfaceNC(prop,phys,grid,ncscratch,Nci,T0,myproc);
  }else{
    if((int)MPI_GetValue(DATAFILE,"init_fs_from_file","ReadGrid",myproc)) {
      MPI_GetFile(filename,DATAFILE,"fs_init_file","InitializePhysicalVariables",myproc);  
      sprintf(str,"%s.%d",filename,myproc);
      fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

      fread(phys->h,sizeof(REAL),grid->Nc,fid);
      
      fclose(fid);

      for(i=0;i<Nc;i++) {
        phys->dhdt[i]=0;
        phys->active[i]=1; // Need to make sure all cells are active by default                                                                                     
      }

    } else {        
      for(i=0;i<Nc;i++) {
	phys->dhdt[i]=0;
	phys->active[i]=1; // Need to make sure all cells are active by default
	phys->h[i]=ReturnFreeSurface(grid->xv[i],grid->yv[i],grid->dv[i]);
	if(phys->h[i]<-grid->dv[i] + DRYCELLHEIGHT){ 
	  phys->h[i]=-grid->dv[i] + DRYCELLHEIGHT;
	  phys->active[i]=0;  
	}
      }
    }
    if(DEBUG) printf("Done with h\n");
  }

  if(DEBUG) printf("Updating vert coords\n");
  // Need to update the vertical grid after updating the free surface.
  // The 1 indicates that this is the first call to UpdateDZ
  if(prop->vertcoord!=1){
   if(prop->vertcoord==5)
     UpdateDZ(grid,phys,prop, 1);

   if(DEBUG) printf("Vertcoordinate basic\n");
   VertCoordinateBasic(grid,prop,phys,myproc,comm);
   if(DEBUG) printf("Done vert coord basic\n");
  }
  else
   UpdateDZ(grid,phys,prop, 1);
  if(DEBUG) printf("Done updating vert coords\n");
  
  // initailize variables to 0 (except for filter "pressure")
  if(DEBUG) printf("Initializing to zero\n");
  for(i=0;i<Nc;i++) {
    phys->zB[i]=phys->zBold[i]=phys->zBold2[i]=phys->zBoffline[i]=0;    
    phys->Erosion[i]=phys->Erosion_old[i]=phys->Deposition[i]=0;    
    phys->w[i][grid->Nk[i]]=0;
    phys->w_s[i][grid->Nk[i]]=0;
    phys->w_sed[i][grid->Nk[i]]=0;
    for(k=0;k<grid->Nk[i];k++) {
      phys->w[i][k]=0;
      phys->w_s[i][k]=0;
      phys->w_sed[i][k]=0;
      phys->q[i][k]=0;
      phys->qc[i][k]=0;      
      phys->s[i][k]=0;
      phys->T[i][k]=0;
      phys->s_old[i][k]=0;
      phys->T_old[i][k]=0;
      phys->s0[i][k]=0;
    }
  }

  if(DEBUG) printf("Initializing u to 0\n");
  for(j=0;j<grid->Ne;j++){
    for(k=0;k<grid->Nke[j];k++){
      phys->u[j][k]=0;  
      phys->u_old[j][k]=0;
      phys->u_old2[j][k]=0;
    }
  }
  // Initialize the temperature, salinity, and background salinity
  // distributions.  Since z is not stored, need to use dz[k] to get
  // z[k].
  if(prop->readSalinity && prop->readinitialnc == 0) {
    stmp = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"InitializePhysicalVariables");
    if(fread(stmp,sizeof(REAL),grid->Nkmax,prop->InitSalinityFID) != grid->Nkmax)
      printf("Error reading stmp first\n");
    fclose(prop->InitSalinityFID);

    for(i=0;i<Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        phys->s[i][k]=stmp[k];
        phys->s0[i][k]=stmp[k];
      }
    SunFree(stmp,grid->Nkmax*sizeof(REAL),"InitializePhysicalVariables");
  } else if(prop->readinitialnc){
     ReturnSalinityNC(prop,phys,grid,ncscratch,Nci,Nki,T0,myproc);
  } else {
    if((int)MPI_GetValue(DATAFILE,"init_s_from_file","ReadGrid",myproc)) {
      MPI_GetFile(filename,DATAFILE,"s_init_file","InitializePhysicalVariables",myproc);  
      sprintf(str,"%s.%d",filename,myproc);  
      fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

      for(i=0;i<Nc;i++) {
	fread(phys->s[i],sizeof(REAL),grid->Nkmax,fid);
	for(k=0;k<grid->Nkmax;k++)
	  phys->s0[i][k]=phys->s[i][k];
      }
      fclose(fid);

      ISendRecvCellData3D(phys->s,grid,myproc,comm);
    } else {
      for(i=0;i<Nc;i++) {
	z = 0;
	zu = 0;
	for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	  z-=grid->dz[k]/2;

	  switch(prop->vertcoord) {
	  case 1:
	    phys->s[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],z);
	    phys->s0[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],z);                

	    break;
	  case 0: case 3:
	    phys->s[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],vert->zc[i][k]);
	    phys->s0[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],vert->zc[i][k]);           

	    break;
	  case 2: //case 4:
	    phys->s[i][k]=IsoReturnSalinity(grid->xv[i],grid->yv[i],z,zu,zu-grid->dz[k],i,k);
	    phys->s0[i][k]=IsoReturnSalinity(grid->xv[i],grid->yv[i],z,zu,zu-grid->dz[k],i,k);
	  case 4:

	    dmax = 500;

	    a = 200;
	    L0 = 500;2500;
	    L = 150e3;

	    zeta = -a*pow(cosh(grid->xv[i]/L0),-2.0);
	    znew = vert->zc[i][k]-zeta*sin(-PI*vert->zc[i][k]/dmax);

	    phys->s[i][k]=ReturnSalinity(grid->xv[i],grid->yv[i],znew); //sig version
	    phys->s0[i][k] = phys->s[i][k];
	    
	    break;
	  default:
	    if(myproc==0) {
	      printf("Error initializing salinity. Unknown vertcoord type %d.\n",prop->vertcoord);
	      MPI_Finalize();
	      exit(EXIT_FAILURE);
	    }

	    break;
	  }
	  z-=grid->dz[k]/2;
	  zu-=grid->dz[k];
	}
      }
    }
  }

  //printf("initializing q");

  //initialize q if desired 
  if((int)MPI_GetValue(DATAFILE,"init_q_from_file","ReadGrid",myproc)) {
    MPI_GetFile(filename,DATAFILE,"q_init_file","InitializePhysicalVariables",myproc);
    sprintf(str,"%s.%d",filename,myproc);
    fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

    for(i=0;i<Nc;i++) {
      fread(phys->q[i],sizeof(REAL),grid->Nkmax,fid);
    }
    fclose(fid);
  }

  //printf("initialized q");

  if(prop->readTemperature && prop->readinitialnc == 0) {
    stmp = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"InitializePhysicalVariables");
    if(fread(stmp,sizeof(REAL),grid->Nkmax,prop->InitTemperatureFID) != grid->Nkmax)
      printf("Error reading stmp second\n");
    fclose(prop->InitTemperatureFID);    

    for(i=0;i<Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
        phys->T[i][k]=stmp[k];

    SunFree(stmp,grid->Nkmax*sizeof(REAL),"InitializePhysicalVariables");
   } else if(prop->readinitialnc){
        ReturnTemperatureNC(prop,phys,grid,ncscratch,Nci,Nki,T0,myproc);
  } else {
    if((int)MPI_GetValue(DATAFILE,"init_T_from_file","ReadGrid",myproc)) {
      MPI_GetFile(filename,DATAFILE,"T_init_file","InitializePhysicalVariables",myproc);  
      sprintf(str,"%s.%d",filename,myproc);  
      fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

      for(i=0;i<Nc;i++) {
	fread(phys->T[i],sizeof(REAL),grid->Nkmax,fid);
      }
      fclose(fid);

      ISendRecvCellData3D(phys->T,grid,myproc,comm);
    } else {    
      if(DEBUG) printf("Init temperature\n");
      for(i=0;i<Nc;i++) {
	z = 0;

	for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	  z-=grid->dz[k]/2;

	  switch(prop->vertcoord) {
	  case 1:
	    phys->T[i][k]=ReturnTemperature(grid->xv[i],grid->yv[i],z,grid->dv[i]);

	    break;
	  case 0: case 3:
	    phys->T[i][k]=ReturnTemperature(grid->xv[i],grid->yv[i],vert->zc[i][k],grid->dv[i]);

	    break;
	  case 2: case 4:
	    phys->T[i][k]=IsoReturnTemperature(grid->xv[i],grid->yv[i],z,grid->dv[i],i,k);

	    break;
	  default:
	    if(myproc==0) {
	      printf("Error initializing temperature. Unknown vertcoord type %d.\n",prop->vertcoord);
	      MPI_Finalize();
	      exit(EXIT_FAILURE);
	    }

	    break;
	  }
	  z-=grid->dz[k]/2;	
	}
      }
    }
  }

  if(DEBUG) printf("Old T and S\n");
  // set the old values s^n-1 T^n-1=s^n T^n at prop->n==1
  if((int)MPI_GetValue(DATAFILE,"init_st_from_file","ReadGrid",myproc)){
    MPI_GetFile(filename,DATAFILE,"s_t1_init_file","InitializePhysicalVariables",myproc);  
    sprintf(str,"%s.%d",filename,myproc);  
    fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

    for(j=0;j<grid->Nc;j++) {
      fread(phys->s_old[j],sizeof(REAL),grid->Nkmax,fid);
    }
    fclose(fid);
    ISendRecvCellData3D(phys->s_old,grid,myproc,comm);

    //still use old T since we're not using T here
    for(i=0;i<Nc;i++)
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
      {
        phys->T_old[i][k]=phys->T[i][k];
        //phys->s_old2[i][k]=phys->s[i][k];
      }
  } else{
  for(i=0;i<Nc;i++)
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
    {
      phys->T_old[i][k]=phys->T[i][k];
      phys->s_old[i][k]=phys->s[i][k];      
        //phys->s_old2[i][k]=phys->s[i][k];      
      }
    }

  if(DEBUG) printf("Init velocity\n");
  // Initialize the velocity field
  if((int)MPI_GetValue(DATAFILE,"init_u_from_file","ReadGrid",myproc)) {
    MPI_GetFile(filename,DATAFILE,"u_init_file","InitializePhysicalVariables",myproc);  
    sprintf(str,"%s.%d",filename,myproc);  
    fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

    for(j=0;j<grid->Ne;j++) {
      fread(phys->u[j],sizeof(REAL),grid->Nkmax,fid);
    }
    fclose(fid);

    ISendRecvEdgeData3D(phys->u,grid,myproc,comm);
    //Set u old and u old 2 (not currently in use)
    if((int)MPI_GetValue(DATAFILE,"init_u_t1_from_file","ReadGrid",myproc)) { 

      //u_old 
      MPI_GetFile(filename,DATAFILE,"u0_tminus1_file","InitializePhysicalVariables",myproc);  
      sprintf(str,"%s.%d",filename,myproc);  
      fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

      for(j=0;j<grid->Ne;j++) {
        fread(phys->u_old[j],sizeof(REAL),grid->Nkmax,fid);
      }
      fclose(fid);
      if(DEBUG) printf("Initialized velocity t-1\n");

      ISendRecvEdgeData3D(phys->u_old,grid,myproc,comm);
      //u_old2
      MPI_GetFile(filename,DATAFILE,"u0_tminus2_file","InitializePhysicalVariables",myproc);  
      sprintf(str,"%s.%d",filename,myproc);  
      fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

      for(j=0;j<grid->Ne;j++) {
        fread(phys->u_old2[j],sizeof(REAL),grid->Nkmax,fid);
      }
      fclose(fid); 
      if(DEBUG) printf("Initialized velocity t-2\n");

      ISendRecvEdgeData3D(phys->u_old2,grid,myproc,comm);

      //Try turning off reading in w. 
      //w
      MPI_GetFile(filename,DATAFILE,"w_init_file","InitializePhysicalVariables",myproc);  
      sprintf(str,"%s.%d",filename,myproc);  
      fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

      //should this be Nkmax + 1
      for(j=0;j<grid->Nc;j++) {
        fread(phys->w[j],sizeof(REAL),grid->Nkmax+1,fid);
      }
      fclose(fid); 
      if(DEBUG) printf("Initialized w \n");
      // for(j=0;j<grid->Nc;j++){
      //   printf("w at top is %f and w at bottom is %f for i=%d \n", phys->w[j][0], phys->w[j][grid->Nk[j]], j);
      // }

      //w_old
      MPI_GetFile(filename,DATAFILE,"w0_tminus1_file","InitializePhysicalVariables",myproc);  
      sprintf(str,"%s.%d",filename,myproc);  
      fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

      for(j=0;j<grid->Nc;j++) {
        fread(phys->w_old[j],sizeof(REAL),grid->Nkmax+1,fid);
      }
      fclose(fid); 
      if(DEBUG) printf("Initialized w t-1\n");
      // for(j=0;j<grid->Nc;j++){
      //   printf("w_old at top is %f and w at bottom is %f for i=%d \n", phys->w_old[j][0], phys->w_old[j][grid->Nk[j]], j);
      // }

      //w_old2
      MPI_GetFile(filename,DATAFILE,"w0_tminus2_file","InitializePhysicalVariables",myproc);  
      sprintf(str,"%s.%d",filename,myproc);  
      fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

      for(j=0;j<grid->Nc;j++) {
        fread(phys->w_old2[j],sizeof(REAL),grid->Nkmax+1,fid);
      }
      fclose(fid); 

      // for(j=0;j<grid->Nc;j++){
      //   printf("w_old2 at top is %f and w at bottom is %f for i=%d \n", phys->w_old2[j][0], phys->w_old2[j][grid->Nk[j]], j);
      // }
      if(DEBUG) printf("Initialized w t-2\n");
    } else {
      for(j=0;j<grid->Ne;j++){
        for(k=0;k<grid->Nke[j];k++) {
          phys->u_old[j][k]=phys->u[j][k];
          phys->u_old2[j][k]=phys->u[j][k];
        }
      }

      for(i=0;i<grid->Nc;i++){
        for(k=0;k<grid->Nk[i];k++) {
          phys->w_old[i][k]=phys->w[i][k];
          phys->w_old2[i][k]=phys->w[i][k];
        }
      }

    }
  } else {
        for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      z = 0;
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      if(nc1==-1)
	nc1=nc2;
      if(nc2==-1)
	nc2=nc1;
      for(k=0;k<grid->Nke[j];k++) {
	z-=grid->dz[k]/2;
	if(prop->vertcoord==1)
	  phys->u[j][k]=ReturnHorizontalVelocity(
						 grid->xe[j],grid->ye[j],grid->n1[j],grid->n2[j],z);
	else  {
	  phys->u[j][k]=ReturnHorizontalVelocity(grid->xe[j],grid->ye[j],grid->n1[j],grid->n2[j],
						 0.5*(vert->zc[nc1][k]+vert->zc[nc2][k]));
	}
	z-=grid->dz[k]/2;
      }
    }
  }

  if(DEBUG) printf("Heat flux\n");
  // Initialise the heat flux arrays
  for(i=0;i<Nc;i++) {
    ktop = grid->ctop[i];
    phys->Tsurf[i] = phys->T[i][ktop];
    phys->dT[i] = 0.001; // Needs to be != 0
    
    for(k=0;k<grid->Nk[i];k++){
	phys->Ttmp[i][k]=phys->T[i][k];
    }
  }


  if(DEBUG) printf("UC\n");
  // Need to compute the velocity vectors at the cell centers based
  // on the initialized velocities at the faces.
  // since subgrid is not allocated yet, so just do the regular perot interpolation here
  ComputeUC(phys->uc, phys->vc, phys, grid, myproc, prop->interp,prop->kinterp,0);
  //added to calculate old us 
  ComputeUCPerot(phys->u_old, phys->uold, phys->vold, phys->h, prop->kinterp, 0, grid);
  //ComputeUC(phys->uold, phys->vold, phys, grid, myproc, prop->interp,prop->kinterp,0);

  if(DEBUG) printf("Send/recv\n");
  // send and receive interprocessor data
  ISendRecvCellData3D(phys->uc,grid,myproc,comm);
  ISendRecvCellData3D(phys->vc,grid,myproc,comm);
  ISendRecvCellData3D(phys->uold,grid,myproc,comm);
  ISendRecvCellData3D(phys->vold,grid,myproc,comm);
  ISendRecvWData(phys->w,grid,myproc,comm);
  ISendRecvWData(phys->w_old,grid,myproc,comm);
  ISendRecvWData(phys->w_old2,grid,myproc,comm);

  if(DEBUG) printf("smin and smax\n");
  // Determine minimum and maximum salinity
  phys->smin=phys->s[0][0];
  phys->smax=phys->s[0][0];

  // overall cells in plan
  for(i=0;i<grid->Nc;i++) {
    // and in the vertical
    for(k=0;k<grid->Nk[i];k++) {
      if(phys->s[i][k]<phys->smin) phys->smin=phys->s[i][k];      
      if(phys->s[i][k]>phys->smax) phys->smax=phys->s[i][k];      
    }
  }

  if(DEBUG) printf("Density\n");
  // Set the density from s and T using the equation of state 
  SetDensity(grid,phys,prop);

  if(DEBUG) printf("Eddy viscosity\n");
  // Initialize the eddy-viscosity and scalar diffusivity
  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i];k++) {
      phys->nu_tv[i][k]=0;
      phys->kappa_tv[i][k]=0;
      phys->SediKappa_tv[i][k]=0;
      phys->nu_lax[i][k]=0;
    }
  }

  if(prop->turbmodel>=1) {
    for(i=0;i<grid->Nc;i++) {
      for(k=0;k<grid->Nk[i];k++) {
        phys->qT[i][k]=0;
        phys->lT[i][k]=0;
        phys->qT_old[i][k]=0;
        phys->lT_old[i][k]=0;
	phys->Cn_l[i][k]=0;
        phys->Cn_q[i][k]=0;
      }
    }
  }

  // Free the scratch array
  if (prop->readinitialnc>0)
      SunFree(ncscratch,Nki*Nci*sizeof(REAL),"InitializePhyiscalVariables");
      //Close the initial condition netcdf file
      //MPI_NCClose(prop->initialNCfileID );

  if(DEBUG) printf("Done with init\n");
}

/*
 * Function: SetDragCoefficients
 * Usage: SetDragCoefficents(grid,phys,prop,myproc);
 * ------------------------------------------
 * Set the drag coefficients based on the log law as well as the applied shear stress.
 * to get CdT and CdB
 */
void SetDragCoefficients(gridT *grid, physT *phys, propT *prop) {
  int i, j, k,nc1,nc2;
  REAL z,dz,u_sum,A_sum,zfb;
  // if set Z0, it will calculate Cd from log law. if Cd is given directly, it will set value to 
  // phys->Cd
  // z0T
  if(prop->z0T==0) 
    for(j=0;j<grid->Ne;j++) 
      phys->CdT[j]=prop->CdT;
  else
    for(j=0;j<grid->Ne;j++) 
      phys->CdT[j]=pow(log(0.5*grid->dzf[j][grid->etop[j]]/phys->z0T[j])/KAPPA_VK,-2);

  // z0B
  if(prop->z0B==0) 
    for(j=0;j<grid->Ne;j++) 
      phys->CdB[j]=prop->CdB;
  else
    for(j=0;j<grid->Ne;j++) 
    {
      if(prop->vertcoord==1)
        if(!prop->subgrid)
          zfb=0.5*grid->dzf[j][grid->Nke[j]-1];
        else
          zfb=0.5*subgrid->dzboteff[j];
      else
        // need new modification for subgrid bathymetry
        zfb=0.5*vert->zfb[j];

      if(grid->Nke[j]>1 && grid->etop[j]!=(grid->Nke[j]-1))
        phys->CdB[j]=pow(log(zfb/phys->z0B[j])/KAPPA_VK,-2);
      else
        phys->CdB[j]=pow((log(2*zfb/phys->z0B[j])+phys->z0B[j]/2/zfb-1)/KAPPA_VK,-2);
    }

  if(prop->subgrid && grid->Nkmax==1)
    if(subgrid->dragpara)
      CalculateSubgridDragCoef(grid,phys,prop);
  
  for(j=0;j<grid->Ne;j++){
    if(prop->vertcoord==1)
      if(!prop->subgrid)
        zfb=0.5*grid->dzf[j][grid->Nke[j]-1];
      else
        zfb=0.5*subgrid->dzboteff[j];
    else
      // need new modification for subgrid bathymetry
      zfb=vert->zfb[j];
  
    if(prop->vertcoord==1)
    {  
      if(2*zfb<prop->BUFFERHEIGHT && grid->etop[j]==(grid->Nke[j]-1))
        phys->CdB[j]=100;  
    }else {
      // for the new vertical coordinate there is always constant layer numbers
      /*
      if(2*zfb<prop->BUFFERHEIGHT)
        phys->CdB[j]=100;
      */
    }
  }

  //ramp up bot stress, but leave alone if CdB=-1 or 100

  for(j=0;j<grid->Ne;j++){
    if(phys->CdB[j]!=-1 && phys->CdB[j]!=100){
      phys->CdB[j]*=erf((prop->rtime)/2830);
    }
  }
}

/*
* Function: InterpDrag(grid,phys,prop,myproc)
* usage: interpolate the value for z0T and z0B
* --------------------------------------------
* Author: Yun Zhang        
* Interpolate z0T and z0B according to z0t.dat and z0b.dat
* Intz0B==1, interpolate value
* Intz0B==2, read center point data directly
*/
static void InterpDrag(gridT *grid, physT *phys, propT *prop,int myproc)
{
   int n, Nb, Nt;
   REAL *xb, *yb, *z0b, *xt, *yt, *z0t;
   char str[BUFFERLENGTH];
   FILE *fid;
   // for z0B
   if(prop->Intz0B==1){
     Nb = MPI_GetSize(prop->INPUTZ0BFILE,"InterpDrag",myproc);
     xb = (REAL *)SunMalloc(Nb*sizeof(REAL),"InterpDrag");
     yb = (REAL *)SunMalloc(Nb*sizeof(REAL),"InterpDrag");
     z0b = (REAL *)SunMalloc(Nb*sizeof(REAL),"InterpDrag");
     fid = MPI_FOpen(prop->INPUTZ0BFILE,"r","InterpDrag",myproc);
     for(n=0;n<Nb;n++) {
       xb[n]=getfield(fid,str);
       yb[n]=getfield(fid,str);
       z0b[n]=getfield(fid,str); // but cd should be positive already
     }
     
     fclose(fid);
     // grid have xe ye
    /* for(n=0;n<grid->Ne;n++){
       // get the middle point of each edge
       xe[n]=0.5*(grid->xp[grid->edges[n*NUMEDGECOLUMNS]]+grid->xp[grid->edges[n*NUMEDGECOLUMNS+1]]);
       ye[n]=0.5*(grid->yp[grid->edges[n*NUMEDGECOLUMNS]]+grid->yp[grid->edges[n*NUMEDGECOLUMNS+1]]);
     }*/

     Interp(xb,yb,z0b,Nb,&(grid->xe[0]), &(grid->ye[0]), &(phys->z0B[0]), grid->Ne, grid->maxfaces);
     free(xb);
     free(yb);
     free(z0b);

   } else if(prop->Intz0B==2){
     sprintf(str,"%s-edge",prop->INPUTZ0BFILE);
     fid = MPI_FOpen(str,"r","InterpDrag",myproc);
     for(n=0;n<grid->Ne;n++) {
       getfield(fid,str);
       getfield(fid,str);
       phys->z0B[n]=getfield(fid,str);
     }
     fclose(fid);
   } else if(prop->Intz0B==0){
     for(n=0;n<grid->Ne;n++)
       phys->z0B[n]=prop->z0B;
   } else if(prop->Intz0B!=0){
     printf("Intz0B=%d, Intz0B can only be 0, 1 and 2\n",prop->Intz0B);
     MPI_Finalize();
     exit(EXIT_FAILURE);
   }

   // for z0T
   if(prop->Intz0T==1){

     Nt = MPI_GetSize(prop->INPUTZ0TFILE,"InterpDrag",myproc); // change to z0TFILE
     xt = (REAL *)SunMalloc(Nt*sizeof(REAL),"InterpDrag");
     yt = (REAL *)SunMalloc(Nt*sizeof(REAL),"InterpDrag");
     z0t = (REAL *)SunMalloc(Nt*sizeof(REAL),"InterpDrag");
     fid = MPI_FOpen(prop->INPUTZ0TFILE,"r","InterpDrag",myproc);
     for(n=0;n<Nt;n++) {
       xt[n]=getfield(fid,str);
       yt[n]=getfield(fid,str);
       z0t[n]=getfield(fid,str); // but cd should be positive already
     }
     fclose(fid);
     Interp(xt,yt,z0t,Nt,&(grid->xe[0]), &(grid->ye[0]),&(phys->z0T[0]),grid->Ne,grid->maxfaces);
     free(xt);
     free(yt);
     free(z0t);

   } else if(prop->Intz0T==2){
     sprintf(str,"%s-edge",prop->INPUTZ0TFILE);
     fid = MPI_FOpen(str,"r","InterpDrag",myproc);
     for(n=0;n<grid->Ne;n++) {
       getfield(fid,str);
       getfield(fid,str);
       phys->z0T[n]=getfield(fid,str);
     }
     fclose(fid);
   } else if(prop->Intz0T==0){
     
     for(n=0;n<grid->Ne;n++)
       phys->z0T[n]=prop->z0T;

   } else if(prop->Intz0T!=0){
     printf("Intz0T=%d, Intz0T can only be 0, 1 and 2\n",prop->Intz0T);
     MPI_Finalize();
     exit(EXIT_FAILURE);
   }
}

/*
* Function OutputDrag(grid,phys,prop,myproc)
* -------------------------------------------
* Author: Yun Zhang
* usage: output the phys->z0t and phys->z0b 
*         all the data is located 
*        at the center of cell
*/
static void OutputDrag(gridT *grid,physT *phys, propT *prop,int myproc, int numprocs, MPI_Comm comm)
{ 
  int n,nf,ne;
  REAL *z0b, *z0t;
  char str[BUFFERLENGTH], str1[BUFFERLENGTH], str2[BUFFERLENGTH];
  FILE *ofile;
  // for z0t
  if(prop->Intz0T==1 || prop->Intz0T==2){
    MPI_GetFile(str1,DATAFILE,"z0TFile","OutputDrag",myproc);
    z0t = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"OutputDrag");
    if(prop->mergeArrays)
      strcpy(str1,str);
    else
      sprintf(str1,"%s.%d",str,myproc);
    if(VERBOSE>2) printf("Outputting %s...\n",str); 
    ofile = MPI_FOpen(str,"w","OutputDrag",myproc);
    for(n=0;n<grid->Nc;n++) {
      z0t[n]=0;
      for(nf=0;nf<grid->nfaces[n];nf++) {
        ne = grid->face[n*grid->maxfaces+nf];
        z0t[n]+=phys->z0T[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
      }
      z0t[n]/=2*grid->Ac[n];
    }
    Write2DData(z0t,prop->mergeArrays,ofile,"Error outputting surface roughness data!\n",
		  grid,numprocs,myproc,comm);
    fclose(ofile);	
    free(z0t);
  }
  
  // for z0b
  if(prop->Intz0B==1 || prop->Intz0B==2){
    MPI_GetFile(str2,DATAFILE,"z0BFile","OutputDrag",myproc);
    z0b = (REAL *)SunMalloc(grid->Nc*sizeof(REAL),"OutputDrag");
    if(prop->mergeArrays)
      strcpy(str2,str);
    else
      sprintf(str2,"%s.%d",str,myproc);
    if(VERBOSE>2) printf("Outputting %s...\n",str); 
    ofile = MPI_FOpen(str,"w","OutputDrag",myproc);
    for(n=0;n<grid->Nc;n++) {
      z0b[n]=0;
      for(nf=0;nf<grid->nfaces[n];nf++) {
        ne = grid->face[n*grid->maxfaces+nf];
        z0b[n]+=phys->z0B[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
      }
      z0b[n]/=2*grid->Ac[n];
    }
    Write2DData(z0b,prop->mergeArrays,ofile,"Error outputting bottom roughness data!\n",
		  grid,numprocs,myproc,comm);
    fclose(ofile);
    free(z0b);
  }
}

/*
 * Function: InitializeVerticalGrid
 * Usage: InitializeVerticalGrid(grid);
 * ------------------------------------
 * Initialize the vertical grid by allocating space for grid->dzz and grid->dzzold
 * This just sets dzz and dzzold to dz since upon initialization dzz and dzzold
 * do not vary in the horizontal.
 *
 * This function is necessary so that gridding in veritcal can be redone each
 * simulation to account for changes in suntans.dat for Nkmax
 *
 */
void InitializeVerticalGrid(gridT **grid,int myproc)
{
  int i, j, k, Nc=(*grid)->Nc, Ne=(*grid)->Ne;

  // initialize in plan for the face (dzf dzfB) and 
  // cell centered (dzz dzzold) quantities
  (*grid)->stairstep = MPI_GetValue(DATAFILE,"stairstep","InitializeVerticalGrid",myproc);
  (*grid)->fixdzz = MPI_GetValue(DATAFILE,"fixdzz","InitializeVerticalGrid",myproc);
  (*grid)->dzsmall = (REAL)MPI_GetValue(DATAFILE,"dzsmall","InitializeVerticalGrid",myproc);
  (*grid)->smoothbot = (REAL)MPI_GetValue(DATAFILE,"smoothbot","InitializeVerticalGrid",myproc);  

  (*grid)->dzf = (REAL **)SunMalloc(Ne*sizeof(REAL *),"InitializeVerticalGrid");
  (*grid)->hf=(REAL *)SunMalloc(Ne*sizeof(REAL),"InitializeVerticalGrid");
  (*grid)->dzfB = (REAL *)SunMalloc(Ne*sizeof(REAL),"InitializeVerticalGrid");
  (*grid)->dzz = (REAL **)SunMalloc(Nc*sizeof(REAL *),"InitializeVerticalGrid");
  (*grid)->dzzold = (REAL **)SunMalloc(Nc*sizeof(REAL *),"InitializeVerticalGrid");
  //(*grid)->dzzold2 = (REAL **)SunMalloc(Nc*sizeof(REAL *),"InitializeVerticalGrid");
  (*grid)->dzbot = (REAL *)SunMalloc(Nc*sizeof(REAL),"InitializeVerticalGrid");


  // initialize over depth for edge-oriented quantities
  for(j=0;j<Ne;j++) {
    (*grid)->dzf[j]=(REAL *)SunMalloc(((*grid)->Nkc[j])*sizeof(REAL),"InitializeVerticalGrid");

    for(k=0;k<(*grid)->Nkmax;k++) {
      (*grid)->dzf[j][k]=(*grid)->dz[k];
    }
  }

  // initialize over depth for cell-centered quantities
  for(i=0;i<Nc;i++) {
    (*grid)->dzz[i]=(REAL *)SunMalloc(((*grid)->Nk[i])*sizeof(REAL),"InitializeVerticalGrid");
    (*grid)->dzzold[i]=(REAL *)SunMalloc(((*grid)->Nk[i])*sizeof(REAL),"InitializeVerticalGrid");
    
    for(k=0;k<(*grid)->Nk[i];k++) {
      (*grid)->dzz[i][k]=(*grid)->dz[k];  
      (*grid)->dzzold[i][k]=(*grid)->dz[k];  
    }
  }
}

/*
 * Function: UpdateDZ
 * Usage: UpdateDZ(grid,phys,0);
 * -----------------------------
 * This function updates the vertical grid spacings based on the free surface and
 * the bottom bathymetry.  That is, if the free surface cuts through cells and leaves
 * any cells dry (or wet), then this function will set the vertical grid spacing 
 * accordingly.
 *
 * If option==1, then it assumes this is the first call and sets dzzold to dzz at
 *   the end of this function.
 * Otherwise it sets dzzold to dzz at the beginning of the function and updates dzz
 *   thereafter.
 *
 */
void UpdateDZ(gridT *grid, physT *phys, propT *prop, int option)
{
  int i, j, k, ne1, ne2, Nc=grid->Nc, Ne=grid->Ne, flag, nc1, nc2;
  REAL z, dzz1, dzz2;

  // don't need to recompute for linearized FS
  if(prop->linearFS) {
    return;
  }

  // If this is not an initial call then set dzzold to store the old value of dzz
  // and also set the etopold and ctopold pointers to store the top indices of
  // the grid.
  if(!option) {
    for(j=0;j<Ne;j++)
      grid->etopold[j]=grid->etop[j];
    for(i=0;i<Nc;i++) {
      grid->ctopold[i]=grid->ctop[i];
      for(k=0;k<grid->ctop[i];k++)
        grid->dzzold[i][k]=0;
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }
  //fixdzz
  if(option==-1)
    for(i=0;i<Nc;i++)
      phys->h[i]=.0; 

  // First set the thickness of the bottom grid layer.  If this is a partial-step
  // grid then the dzz will vary over the horizontal at the bottom layer.  Otherwise,
  // the dzz at the bottom will be equal to dz at the bottom.
  for(i=0;i<Nc;i++) {
    z = 0;
    for(k=0;k<grid->Nk[i];k++)
      z-=grid->dz[k];
    grid->dzz[i][grid->Nk[i]-1]=grid->dz[grid->Nk[i]-1]+grid->dv[i]+z;
  }

  // Loop through and set the vertical grid thickness when the free surface cuts through 
  // a particular cell.
  if(grid->Nkmax>1) {
    for(i=0;i<Nc;i++) {
      z = 0;
      flag = 0;
      for(k=0;k<grid->Nk[i];k++) {
        z-=grid->dz[k];
        if(phys->h[i]>=z) 
          if(!flag) {
            if(k==grid->Nk[i]-1) {
              grid->dzz[i][k]=phys->h[i]+grid->dv[i];
              grid->ctop[i]=k;
            } else if(phys->h[i]==z) {
              grid->dzz[i][k]=0;
              grid->ctop[i]=k+1;
            } else {
              grid->dzz[i][k]=phys->h[i]-z;
              grid->ctop[i]=k;
            }
            flag=1;
          } else {
            if(k==grid->Nk[i]-1) 
              grid->dzz[i][k]=grid->dz[k]+grid->dv[i]+z;
            else 
              if(z<-grid->dv[i])
                grid->dzz[i][k]=0;
              else 
                grid->dzz[i][k]=grid->dz[k];
          } 
        else {
          // change 10/01/2014
          if(flag==0 && k==(grid->Nk[i]-1)){
            grid->dzz[i][k]=grid->dv[i]+phys->h[i];
            grid->ctop[i]=k;
          }  
          else
            grid->dzz[i][k]=0;
        }
      }
    }
  } else 
    for(i=0;i<Nc;i++) 
      grid->dzz[i][0]=grid->dv[i]+phys->h[i];

  // Now set grid->etop and ctop which store the index of the top cell  
  for(j=0;j<grid->Ne;j++) {
    ne1 = grid->grad[2*j];
    ne2 = grid->grad[2*j+1];
    if(ne1 == -1)
      grid->etop[j]=grid->ctop[ne2];
    else if(ne2 == -1)
      grid->etop[j]=grid->ctop[ne1];
    else if(grid->ctop[ne1]<grid->ctop[ne2])
      grid->etop[j]=grid->ctop[ne1];
    else
      grid->etop[j]=grid->ctop[ne2];
  }

  // If this is an initial call set the old values to the new values and
  // Determine the bottom-most flux-face height in the absence of h
  // Check the smallest dzz and set to minimum  
  if(option==-1) {
    for(i=0;i<Nc;i++){
      k=grid->Nk[i]-1;      
      grid->dzbot[i]=grid->dzz[i][k];
      if(!grid->stairstep && grid->fixdzz )   
        if(grid->dzz[i][k]<grid->dz[k]*grid->dzsmall) {
          grid->dv[i]+= (grid->dz[k]*grid->dzsmall-grid->dzz[i][k]);	
          grid->dzz[i][k]=grid->dz[k]*grid->dzsmall;
        }

    }
  }

  if(option==1) {    
    for(j=0;j<Ne;j++) 
      grid->etopold[j]=grid->etop[j];
    for(i=0;i<Nc;i++) {
      grid->ctopold[i]=grid->ctop[i];
      for(k=0;k<grid->Nk[i];k++)
        grid->dzzold[i][k]=grid->dzz[i][k];
    }

    for(j=0;j<Ne;j++) {
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      dzz1 = grid->dzz[nc1][grid->Nk[nc1]-1];
      dzz2 = grid->dzz[nc2][grid->Nk[nc2]-1];
      z=0;
      for(k=0;k<grid->Nke[j]-1;k++)
        z-=grid->dz[k];
      if(phys->h[nc1]<z) 
        dzz1 = dzz1-z+phys->h[nc1];
      if(phys->h[nc2]<z) 
        dzz2 = dzz2-z+phys->h[nc2];
      grid->dzfB[j] = Min(dzz1,dzz2);
    }
  }
}

/*
 * Function: DepthFromDZ
 * Usage: z = DepthFromDZ(grid,phys,i,k);
 * --------------------------------------
 * Return the depth beneath the free surface at location i, k.
 *
 */
REAL DepthFromDZ(gridT *grid, physT *phys, int i, int kind) {
  if(i==-1) {
    printf("!!Error with pointer => h[-1]!!\n");
    //return NAN;  // not consistent with all C compilers
    return -1;
  }
  else {
    int k;
    REAL z = phys->h[i]-0.5*grid->dzz[i][grid->ctop[i]];
    for(k=grid->ctop[i];k<kind;k++) {
      z-=0.5*grid->dzz[i][k-1];
      z-=0.5*grid->dzz[i][k];
    }
    //  printf("DepthfromDZ done\n");
    return z;
  }
}

/*
 * Function: Solve
 * Usage: Solve(grid,phys,prop,myproc,numprocs,comm);
 * --------------------------------------------------
 * This is the main solving routine and is called from suntans.c.
 *
 */
void Solve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{

  int i,k,j, n, blowup=0,ne,id,nc1,nc2,nf;
  REAL t0,sum1,v_average,flux,normal;
  metinT *metin;
  metT *met;
  averageT *average;
  int BCbot, BCtop;
  
  // Compute the initial quantities for comparison to determine conservative properties
  prop->n=0;
  // this make sure that we aren't loosing mass/energy
  if(DEBUG) printf("Debugging: ComputeConservatives...\n");
  ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);

  // Print out memory usage per processor and total memory if this is the first time step
  if(VERBOSE>2) MemoryStats(grid,myproc,numprocs,comm);
  // initialize theta0
  prop->theta0=prop->theta;

  // initialize the timers
  t_start=Timer();
  t_source=t_predictor=t_nonhydro=t_turb=t_transport=t_io=t_comm=t_check=0;

  // Set all boundary values at time t=nstart*dt;
  prop->n=prop->nstart;
  // initialize the time (often used for boundary/initial conditions)
  prop->rtime=prop->nstart*prop->dt;
 
  // Initialise the netcdf time (moved to InitializePhysicalVariables)
  // Get the toffSet property
  //printf("myproc: %d, starttime: %s\n",prop->starttime);
  //prop->toffSet = getToffSet(prop->basetime,prop->starttime);
  //prop->nctime = prop->toffSet*86400.0 + prop->nstart*prop->dt;
  //printf("myproc: %d, toffSet = %f (%s, %s)\n",myproc,prop->toffSet,&prop->basetime,&prop->starttime);

  // Initialise the boundary data from a netcdf file
  if(prop->netcdfBdy==1){
    AllocateBoundaryData(prop, grid, &bound, myproc,comm);
    InitBoundaryData(prop, grid, myproc,comm);
  }

  // get the boundary velocities (boundaries.c)
  if(DEBUG) printf("Debugging: BoundaryVelocities...\n");
  BoundaryVelocities(grid,phys,prop,myproc, comm); 
  // get the openboundary flux (boundaries.c)
  if(DEBUG) printf("Debugging: OpenBoundaryFluxes...\n");
  OpenBoundaryFluxes(NULL,phys->u,NULL,grid,phys,prop);
  // get the boundary scalars (boundaries.c)
  if(DEBUG) printf("Debugging: BoundaryScalars...\n");
  BoundaryScalars(grid,phys,prop,myproc,comm);
  // set the height of the face bewteen cells to compute the flux
  if(DEBUG) printf("Debugging: TVDFluxHeight...\n");
  
  //send receive dzz first                                                                                                                                          
  ISendRecvCellData3D(grid->dzz,grid,myproc,comm);
  ISendRecvCellData3D(grid->dzzold,grid,myproc,comm);

  if(prop->vertcoord!=1 && prop->vertcoord!=5)
    TvdFluxHeight(grid, phys, prop, vert->dzfmeth,comm, myproc);
    //uncommented on sherlock, but not locally...
  if(DEBUG) printf("Debugging: SetFluxHeight...\n");  

  //send receive dzz first
  //ISendRecvCellData3D(grid->dzz,grid,myproc,comm);
  //ISendRecvCellData3D(grid->dzzold,grid,myproc,comm);

  SetFluxHeight(grid,phys,prop,vert->dzfmeth,comm,myproc);
  ISendRecvEdgeData3D(grid->dzf,grid,myproc,comm);

  // find the bottom layer number which zfb>Bufferheight only works for new vertical coordinate
  if(DEBUG) printf("Debugging: FindBottomLayer...\n");
  if(prop->vertcoord!=1)
    FindBottomLayer(grid,prop,phys,myproc);

  // get the boundary velocities (boundaries.c)
  if(DEBUG) printf("Debugging: BoundaryVelocities...\n");  
  BoundaryVelocities(grid,phys,prop,myproc, comm); 
  // get the openboundary flux (boundaries.c)
  if(DEBUG) printf("Debugging: OpenBoundaryFluxes...\n");  
  OpenBoundaryFluxes(NULL,phys->u,NULL,grid,phys,prop);

  if(DEBUG) printf("Debugging: InitiaizeMerging...\n");  
  // Set up arrays to merge output
  if(prop->mergeArrays) {
    if(VERBOSE>2 && myproc==0) printf("Initializing arrays for merging...\n");
    InitializeMerging(grid,prop->outputNetcdf,numprocs,myproc,comm);
  }

  // culvert model
  if(prop->culvertmodel){
    if(DEBUG) printf("Debugging: SetupCulvertmodel...\n");    
    if(myproc==0)
      printf("\n\nculvert model has beed started\n\n");
    SetupCulvertmodel(grid,phys,prop,myproc);
  }

  // interp hmarsh and CdV and get hmarshleft and marshtop if consider marshmodel
  if(prop->marshmodel){
    if(DEBUG) printf("Debugging: SetupMarshmodel...\n");        
    SetupMarshmodel(grid,phys,prop,myproc,numprocs,comm);
    SetMarshTop(grid,phys,myproc);
  }

  // use subgrid method 
  if(prop->subgrid)
  {
    if(DEBUG) printf("Debugging: Setting up subgrid...\n");            
    SubgridBasic(grid,phys,prop,myproc,numprocs,comm);
    UpdateSubgridVeff(grid, phys, prop, myproc);
    UpdateSubgridFluxHeight(grid, phys, prop, myproc);
    UpdateSubgridAceff(grid, phys, prop, myproc);
    UpdateSubgridHeff(grid,phys,prop,myproc);
    UpdateSubgridVerticalAceff(grid, phys, prop, 0, myproc);
    if(prop->culvertmodel)
      SubgridCulverttopArea(grid, prop, myproc);
    SubgridFluxCheck(grid, phys, prop,myproc);
    //prop->thetaM=-1;
    //printf("Subgrid module is turned on, set thetaM=-1 which means vertical momentum advection calculation is explicit\n");
  }


  // interpolate z0B and z0T if Intz0t or Intz0B not zero
  if(DEBUG) printf("Debugging: InterpDrag...\n");              
  InterpDrag(grid,phys,prop,myproc);

  // set the drag coefficients for bottom friction
  if(DEBUG) printf("Debugging: SetDragCoefficients...\n");                
  SetDragCoefficients(grid,phys,prop);

  // add culvert part change drag coefficient for culvert part
  if(prop->culvertmodel) {
    if(DEBUG) printf("Debugging: SetCulvertDragCoefficient...\n");                  
    SetCulvertDragCoefficient(grid, phys, prop, myproc);
  }

  // output the drag coefficient for sunplot
  // including CdV z0B and z0T
  // if intz0B and intz0T is zero, z0B=prop->z0B and z0T=prop->z0T
  // output drag z0B zOT and CdV

  if(DEBUG) printf("Debugging: OutputDrag...\n");                    
  OutputDrag(grid,phys,prop,myproc,numprocs,comm);

  // for laxWendroff and central differencing- compute the numerical diffusion 
  // coefficients required for stability
  if(DEBUG) printf("Debugging: LaxWendroff...\n");                      
  if(prop->laxWendroff && prop->nonlinear==2) LaxWendroff(grid,phys,prop,myproc,comm);

  // Initialize the Sponge Layer
  if(DEBUG) printf("Debugging: InitSponge...\n");                        
  //printf("started init sponge \n");
  InitSponge(grid,myproc);
  //printf("finished init sponge \n");


  // Initialise the meteorological forcing input fields
  if(prop->metmodel>0)
  {
      if(DEBUG) printf("Debugging: Initializing met model...\n");                            
      if (prop->gamma==0.0){
        if(myproc==0) 
          printf("Warning gamma must be > 1 for heat flux model.\n");
      }else{
        if(myproc==0) 
          printf("Initial temperature = %f.\n",phys->T[0][0]);
      }
      AllocateMetIn(prop,grid,&metin,myproc);
      AllocateMet(prop,grid,&met,myproc);
      InitialiseMetFields(prop, grid, metin, met,myproc);
      
      // Initialise the heat flux variables
      updateMetData(prop, grid, metin, met, myproc, comm); 
     
    if(prop->metmodel>=2) {      
      updateAirSeaFluxes(prop, grid, phys, met, phys->T);

      //Communicate across processors
      ISendRecvCellData2D(met->Hs,grid,myproc,comm);
      ISendRecvCellData2D(met->Hl,grid,myproc,comm);
      ISendRecvCellData2D(met->Hsw,grid,myproc,comm);
      ISendRecvCellData2D(met->Hlw,grid,myproc,comm);
      ISendRecvCellData2D(met->tau_x,grid,myproc,comm);
      ISendRecvCellData2D(met->tau_y,grid,myproc,comm);
    }
  }

  // Initialise the output netcdf file metadata
  if(prop->outputNetcdf==1 &&prop->mergeArrays==0){
    if(DEBUG) printf("Debugging: InitialiseOutputNCugrid...\n");                                
    InitialiseOutputNCugrid(prop, grid, phys, met, myproc);
  }

  // Initialise the average arrays and netcdf file
  if(prop->calcaverage>0){
    if(DEBUG) printf("Debugging: Initializing averaging variables...\n");                                    
    AllocateAverageVariables(grid,&average,prop);
    ZeroAverageVariables(grid,average,prop);
    if(prop->mergeArrays==0)
      InitialiseAverageNCugrid(prop, grid, average, myproc);
  }
  // get the windstress (boundaries.c) - this needs to go after met data allocation -MR
  if(DEBUG) printf("Debugging: Initializing wind stress...\n");                                      
  WindStress(grid,phys,prop,met,myproc);

  // main time loop
  for(n=prop->nstart+1;n<=prop->nsteps+prop->nstart;n++) {
    prop->n = n;
    // compute the runtime 
    prop->rtime = n*prop->dt;
    
    // netcdf file time
    prop->nctime = prop->toffSet*86400.0 + n*prop->dt;
    //prop->nctime +=  n*prop->dt;
    //prop->nctime += prop->rtime;
    // Set nsteps<0 for debugging i/o without hydro/seds/etc...
    if(prop->nsteps>0) {

      // Ramp down theta from 1 to the value specified in suntans.dat over
      // the time thetaramptime specified in suntans.dat to damp out transient
      // oscillations
      if(prop->thetaramptime!=0)
        prop->theta=(1-exp(-prop->rtime/prop->thetaramptime))*prop->theta0+
          exp(-prop->rtime/prop->thetaramptime);

      if(prop->movingBed){
	UpdateBottomHeight(phys->zB,phys->zBold,phys->zBold2,phys->zBoffline,grid,prop,phys,myproc,comm);
	//printf("updated bed \n");
      }
      // Compute the horizontal source term phys->utmp which contains the explicit part
      // or the right hand side of the free-surface equation. 
      // begin the timer
      t0=Timer();
      /*
      // get the flux height (since free surface is changing) which is stored in dzf
      if(prop->vertcoord!=1 && prop->vertcoord!=5) {
	if(DEBUG) printf("Debugging - Main time loop: TvdFluxHeight...\n");	
        TvdFluxHeight(grid, phys, prop, vert->dzfmeth,comm, myproc);
      }
      if(DEBUG) printf("Debugging - Main time loop: SetFluxHeight...\n");	
      SetFluxHeight(grid,phys,prop,vert->dzfmeth,comm,myproc);
      */

      //need new function here to calculate source terms for HorizontalSource and Wpredictor BEFORE old variables are written over
      if(prop->n==1 && prop->readOldVelocity){
      //if(0){  
	//printf("reached \n");

	//share these amongst the cells before doing source terms, need dzz[nc1] and [nc2] here 
	//ISendRecvCellData3D(grid->dzz,grid,myproc,comm);
	//ISendRecvCellData3D(grid->dzzold,grid,myproc,comm);

        HorizontalSourceTerms_n1(grid,phys,prop,myproc,numprocs,comm);
        //printf("finished here \n");
        WSourceTerms(grid,phys,prop,myproc,numprocs,comm);
        //printf("and finished here \n");


	  ComputeOmegaOlds(grid, prop, phys,comm,myproc);
          ISendRecvWData(vert->omega,grid,myproc,comm);
          ISendRecvWData(vert->omega_old,grid,myproc,comm);
	  //printf("finished here \n"); 

        // for(i=0;i<grid->Nc;i++) {
        //   for(k=0;k<grid->Nk[i]+1;k++) 
        //   {
        //     phys->w_im[i][k]=prop->imfac2*phys->w_old[i][k]+prop->imfac3*phys->w_old2[i][k]+prop->imfac1*phys->w[i][k];
        //   }
        //   phys->w_im[i][k]=0;
        // }

        //printf("imfac1 = %f, imfac2 = %f, imfac3 = %f \n", prop->imfac1, prop->imfac2, prop->imfac3);


        // //also do omega here. 
        // printf("right before this gets called \n");
        // ComputeOmegaOlds(grid, prop, phys, myproc);
        // printf("right after this gets called \n");

      }

      
      // Store the old velocity and scalar fields
      // Store the old values of s, u, and w into stmp3, u_old, and w_old
      if(DEBUG) printf("Debugging - Main time loop: StoreVariables...\n");	
      StoreVariables(grid,phys,prop);
      
      // store old omega for the new vertical coordinate
      if(prop->vertcoord!=1) {
	if(DEBUG) printf("Debugging - Main time loop: StoreVertVariables...\n");		
        StoreVertVariables(grid,phys);
      }

      // find the bottom layer number which zfb>Bufferheight only works for new vertical coordinate
      if(prop->vertcoord!=1) {
	if(DEBUG) printf("Debugging - Main time loop: FindBottomLayer...\n");			
        FindBottomLayer(grid,prop,phys,myproc);
      }
      //printf("finished here \n");
     
      // compute CdB and CdT
      if(DEBUG) printf("Debugging - Main time loop: SetDragCoefficients...\n");			      
      SetDragCoefficients(grid,phys,prop);
      //printf("finished here \n");

      // use subgrid method
      if(prop->subgrid) {
	if(DEBUG) printf("Debugging - Main time loop: UpdateSubgridFluxHeight...\n");			      	
        UpdateSubgridFluxHeight(grid, phys, prop, myproc);
      }

      if(prop->culvertmodel)
      {
        // change phy->h back to pressure field
	if(DEBUG) printf("Debugging - Main time loop: StoreCulverPressure...\n");			      		
        StoreCulvertPressure(phys->h, grid->Nc, 0, myproc);
        // add culvert part change drag coefficient for culvert part
	if(DEBUG) printf("Debugging - Main time loop: SetCulvertDragCoefficient...\n");
	SetCulvertDragCoefficient(grid, phys, prop, myproc);  
      }

      // laxWendroff central differencing
      if(prop->laxWendroff && prop->nonlinear==2) {
	if(DEBUG) printf("Debugging - Main time loop: LaxWendroff...\n");	
        LaxWendroff(grid,phys,prop,myproc,comm);
      }
   
      // compute the horizontal source terms (like 
      /* 
       * 1) Old nonhydrostatic pressure gradient with theta m ethod
       * 2) Coriolis terms with AB2
       * 3) Baroclinic term with AB2
       * 4) Horizontal and vertical advection of horizontal momentum with AB2
       * 5) Horizontal laminar+turbulent diffusion of horizontal momentum
       */

      // calculate the preparation for HorizontalSource function due to the new vertical coordinate
      // time comsuming function!
      if(prop->vertcoord!=1) {
	if(DEBUG) printf("Debugging - Main time loop: VertCoordinateHorizontalSource...\n");
        VertCoordinateHorizontalSource(grid, phys, prop, myproc, numprocs, comm);
      }
      if(DEBUG) printf("Debugging - Main time loop: HorizontalSource...\n");      
      HorizontalSource(grid,phys,prop,myproc,numprocs,comm);
      //printf("and finished here \n");

      // add wave part 
      if(prop->wavemodel) {
	if(DEBUG) printf("Debugging - Main time loop: UpdateWave...\n");      	
        UpdateWave(grid, phys, prop, comm, blowup, myproc, numprocs);
      }

      // compute the time required for the source
      t_source+=Timer()-t0;

      // Use the explicit part created in HorizontalSource and solve for the 
      // free-surface
      // and hence compute the predicted or hydrostatic horizontal velocity field.  Then
      // send and receive the free surface interprocessor boundary data 
      // to the neighboring processors.
      // The predicted horizontal velocity is now in phys->u
      t0=Timer();

      // compute U^* and h^* (Eqn 40 and Eqn 31)
      if(DEBUG) printf("Debugging - Main time loop: Upredictor...\n");
      UPredictor(grid,phys,prop,myproc,numprocs,comm); //NOW we have new dzz, so dzzold = n, dzznew=n+1
      //printf("and finished here \n"); 

      ISendRecvCellData2D(phys->h_old,grid,myproc,comm);
      ISendRecvCellData2D(phys->h,grid,myproc,comm);

      //printf("through u pred \n");
      //now we have u*, dzf, and dzz

      t_predictor+=Timer()-t0;
      t0=Timer();
      if(DEBUG) printf("Debugging - Main time loop: CheckDZ...\n");      	            
      blowup = CheckDZ(grid,phys,prop,myproc,numprocs,comm);
      t_check+=Timer()-t0;

      // apply continuity via Eqn 82
      if(prop->vertcoord==1)
      {
        // w_im is calculated
	if(DEBUG) printf("Debugging - Main time loop: Continuity...\n");      	            	
        Continuity(phys->wnew,grid,phys,prop);	
        ISendRecvWData(phys->wnew,grid,myproc,comm);
      } else {
        // here only calculates the omega value which means the w* is still unknown.
        // the reason is that w is not used for hydrostatic calculation and scalar transport
        // after this calculation vert->omega=omega* vert->omega_old=omega^n vert->omega_old2=omega^n-1
        // omega_im is calculated
	if(DEBUG) printf("Debugging - Main time loop: LayerAveragedContinuity...\n");

        if(0){
          //also do omega here after upredictor. 
          printf("right before this gets called \n");
          // for(i=0;i<grid->Nc;i++){
          //   for(k=0;k<grid->Nk[i]+1;k++){
          //     printf("omega is %f, omega_old is %f, omega_old2 is %f, and omega_im is %f for i =%d, k=%d \n", vert->omega[i][k], 
          //       vert->omega_old[i][k], vert->omega_old2[i][k], vert->omega_im[i][k], i, k);
          // }
          // }
          ComputeOmegaOlds(grid, prop, phys, comm, myproc);
          printf("right after this gets called \n");
        //   for(i=0;i<grid->Nc;i++){
        //     for(k=0;k<grid->Nk[i]+1;k++){
        //       printf("after omega is %f, omega_old is %f, omega_old2 is %f, and omega_im is %f for i =%d, k=%d \n", vert->omega[i][k], 
        //         vert->omega_old[i][k], vert->omega_old2[i][k], vert->omega_im[i][k], i, k);
        //   }
        // }
          ISendRecvWData(vert->omega_old,grid,myproc,comm);
          ISendRecvWData(vert->omega_old2,grid,myproc,comm);

        }
	
        LayerAveragedContinuity(vert->omega,grid,prop,phys,comm,myproc);
        ISendRecvWData(vert->omega,grid,myproc,comm);

        //and now omega and omega_im
      }

      //after this point, have omega and h based on u predictor

      t0=Timer();
      // calculate flux in/out to each cell to ensure bounded scalar concentration under subgrid
      //if(prop->subgrid)
        //SubgridFluxCheck(grid, phys, prop,myproc);

      // Compute the eddy viscosity
      t0=Timer();

      if(DEBUG) printf("Debugging - Main time loop: EddyViscosity...\n");      	            		      
      if(prop->vertcoord==1)
        EddyViscosity(grid,phys,prop,phys->w_im,comm,myproc);
      else
        EddyViscosity(grid,phys,prop,vert->omega_im,comm,myproc);

      t_turb+=Timer()-t0;

      // Update the meteorological data
      if(prop->metmodel>0){
	if(DEBUG) printf("Debugging - Main time loop: UpdateMetData...\n");      	            		      	
        updateMetData(prop, grid, metin, met, myproc, comm);	
        //if(prop->metmodel==2){
        //    updateAirSeaFluxes(prop, grid, phys, met, phys->T);
        //}
      }
      
     // Update the age (passive) tracers
      if(prop->calcage>0){
	if(DEBUG) printf("Debugging - Main time loop: UpdateAge...\n");
	UpdateAge(grid,phys,prop,comm,myproc);
      }
     
      // Update the temperature only if gamma is nonzero in suntans.dat
      if(DEBUG && prop->gamma) printf("Debugging - Main time loop: Update Temperature...\n");	
      if(prop->gamma && prop->advectTemperature) {
        t0=Timer();
        getTsurf(grid,phys); // Find the surface temperature
        HeatSource(phys->wtmp,phys->stmp4,grid,phys,prop,met, myproc, comm);
        if(prop->vertcoord==1){
          UpdateScalars(grid,phys,prop,phys->w_im,phys->T_SfH_tm1,phys->T_SfH_tm2, phys->T_SfHv_t, phys->T_SfHv_tm1,phys->T,phys->T_old,phys->boundary_T,1,1,phys->Cn_R,
                prop->kappa_T,prop->kappa_TH,phys->kappa_tv,prop->theta,
			phys->stmp4,phys->wtmp,NULL,NULL,0,0,comm,myproc,0,0,prop->TVDtemp);
          getchangeT(grid,phys); // Get the change in surface temp
        }else{
	  //printf("pre advect T \n");
          UpdateScalars(grid,phys,prop,vert->omega_im,phys->T_SfH_tm1,phys->T_SfH_tm2, phys->T_SfHv_t, phys->T_SfHv_tm1,phys->T,phys->T_old,phys->boundary_T,1,1,phys->Cn_R, 
            prop->kappa_T,prop->kappa_TH,phys->kappa_tv,prop->theta,
			phys->stmp4,phys->wtmp,NULL,NULL,0,0,comm,myproc,0,0,prop->TVDtemp);
          getchangeT(grid,phys);
          //printf("post advect T \n");
        }
        ISendRecvCellData3D(phys->T,grid,myproc,comm);	
        ISendRecvCellData3D(phys->Ttmp,grid,myproc,comm);
        ISendRecvCellData2D(phys->dT,grid,myproc,comm);
        ISendRecvCellData2D(phys->Tsurf,grid,myproc,comm);

        t_transport+=Timer()-t0;
      }

      // Update the air-sea fluxes --> these are used for the previous time step source term and for the salt flux implicit term (salt tracer solver therefore needs to go next)
      if(prop->metmodel>=2){
	if(DEBUG) printf("Debugging - Main time loop: updateAirSeaFluxes...\n");		
        updateAirSeaFluxes(prop, grid, phys, met, phys->T);
        //Communicate across processors
        ISendRecvCellData2D(met->Hs,grid,myproc,comm);
        ISendRecvCellData2D(met->Hl,grid,myproc,comm);
        ISendRecvCellData2D(met->Hsw,grid,myproc,comm);
        ISendRecvCellData2D(met->Hlw,grid,myproc,comm);
        ISendRecvCellData2D(met->tau_x,grid,myproc,comm);
        ISendRecvCellData2D(met->tau_y,grid,myproc,comm);
      }

      // Update the salinity only if beta is nonzero in suntans.dat
      if(DEBUG && prop->beta) printf("Debugging - Main time loop: Update Salinity...\n");			  	
      if(prop->beta && prop->advectSalinity) {
        t0=Timer();
        if(prop->metmodel>0){
            SaltSource(phys->wtmp,phys->stmp4,grid,phys,prop,met);
            if(prop->vertcoord==1)
	    UpdateScalars(grid,phys,prop,phys->w_im,phys->S_SfH_tm1,phys->S_SfH_tm2, phys->S_SfHv_t, phys->S_SfHv_tm1,phys->s,phys->s_old,phys->boundary_s,1,1,phys->Cn_T,
                prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->theta,
			  phys->stmp4,phys->wtmp,NULL,NULL,0,0,comm,myproc,1,1,prop->TVDsalt);
            else
	    UpdateScalars(grid,phys,prop,vert->omega_im,phys->S_SfH_tm1,phys->S_SfH_tm2, phys->S_SfHv_t, phys->S_SfHv_tm1,phys->s,phys->s_old,phys->boundary_s,1,1,phys->Cn_T,
                prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->theta,
			  phys->stmp4,phys->wtmp,NULL,NULL,0,0,comm,myproc,1,1,prop->TVDsalt);        
        }else{ 
            if(prop->vertcoord==1)
	    UpdateScalars(grid,phys,prop,phys->w_im,phys->S_SfH_tm1,phys->S_SfH_tm2, phys->S_SfHv_t, phys->S_SfHv_tm1,phys->s,phys->s_old,phys->boundary_s,1,1,phys->Cn_T,
                prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->theta,
			  NULL,NULL,NULL,NULL,0,0,comm,myproc,1,1,prop->TVDsalt);
            else
	    UpdateScalars(grid,phys,prop,vert->omega_im,phys->S_SfH_tm1,phys->S_SfH_tm2, phys->S_SfHv_t, phys->S_SfHv_tm1,phys->s,phys->s_old,phys->boundary_s,1,1,phys->Cn_T,
                prop->kappa_s,prop->kappa_sH,phys->kappa_tv,prop->theta,
			  NULL,NULL,NULL,NULL,0,0,comm,myproc,0,0,prop->TVDsalt);
	}
	//printf("Debugging - Main time loop: through Update Salinity scalars...\n");
        ISendRecvCellData3D(phys->s,grid,myproc,comm);
	//MPI_Barrier(comm); //added just in case
	//printf("Debugging - Main time loop: sent and recieved salinity...\n");
      
        if(prop->metmodel>0){
          //Communicate across processors
          ISendRecvCellData2D(met->EP,grid,myproc,comm);
        }

        t_transport+=Timer()-t0;
      }

      // Compute sediment transport when prop->computeSediments=1
      if(prop->computeSediments){
	if(DEBUG) printf("Debugging - Main time loop: Sediment transport...\n");		
        t0=Timer(); 
        // update uc without non-hydrostatic pressure
        // ComputeUC(phys->uc, phys->vc, phys,grid, myproc, prop->interp,prop->kinterp,prop->subgrid);
        ComputeSediments(grid,phys,prop,myproc,numprocs,blowup,comm);
        t_transport+=Timer()-t0;
      }

      //now we've done scalar advection, advect momentum
      // U Advection.  Both blocks stand down when nonlinear!=0, in which case the stock
      // Eulerian advection in HorizontalSource() has already advected horizontal momentum.
      if(U_MOMENTUM_ADVECTION==0 && !prop->nonlinear) {
	int i, j, iptr, jptr, k, nc1, nc2;
	for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	  i = grid->cellp[iptr];
	  
	  for(k=0;k<grid->Nk[i];k++) 
	    phys->wtmp[i][k]=phys->stmp3[i][k]=phys->Cn_T[i][k]=0;
	  phys->wtmp[i][grid->Nk[i]]=0;
	}

  //rethink this. not really using the old vel here, no?
  //uncomment compute ucperot and switch the comments in for loop below to go back to old way

	// uc, vc, uold, vold store velocites at step ex
	//save uc from previous timestep for update scalar
  //uold, vold get set to uc, vc below, so need the previous timestep saved here 
	for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //phys->uold2[i][k]=phys->uc[i][k];
      //phys->vold2[i][k]=phys->vc[i][k];
      phys->uold2[i][k]=phys->uold[i][k];
      phys->vold2[i][k]=phys->vold[i][k];
	  }
  }

  //here do quick instead.

	
  ComputeUCPerot(phys->u_old,phys->uold2,phys->vold2,phys->h_old,prop->kinterp,prop->subgrid,grid);
	ComputeUCPerot(phys->u_ex,phys->uc,phys->vc,phys->h,prop->kinterp,prop->subgrid,grid);
	ISendRecvCellData3D(phys->uc,grid,myproc,comm);
	ISendRecvCellData3D(phys->vc,grid,myproc,comm);
      
	for(i=0;i<grid->Nc;i++) {
	  for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	    phys->uold[i][k]=phys->uc[i][k];
	    phys->vold[i][k]=phys->vc[i][k];
	  }
	}

	ISendRecvCellData3D(phys->uold2,grid,myproc,comm);
  ISendRecvCellData3D(phys->vold2,grid,myproc,comm);
        
	ISendRecvCellData3D(phys->uold,grid,myproc,comm);
	ISendRecvCellData3D(phys->vold,grid,myproc,comm);

	if(prop->CdB==-1)
	  BCbot = 2; //no slip
	else
	  BCbot = 1; //no flux

	BCtop = 1;

	if(prop->vertcoord==1) {
	  UpdateScalars(grid,phys,prop,phys->w_im,phys->u_SfH_tm1,phys->u_SfH_tm2, phys->u_SfHv_t, phys->u_SfHv_tm1,phys->uc,phys->uold2,phys->boundary_u,BCtop, BCbot, phys->Cn_T,
			0,0,NULL,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,1,1,prop->TVDmomentum);
	  UpdateScalars(grid,phys,prop,phys->w_im,phys->v_SfH_tm1,phys->v_SfH_tm2, phys->v_SfHv_t, phys->v_SfHv_tm1,phys->vc,phys->vold2,phys->boundary_v,BCtop,BCbot,phys->Cn_T,
			0,0,NULL,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,1,1,prop->TVDmomentum);
	} else {
	  UpdateScalars(grid,phys,prop,vert->omega_im,phys->u_SfH_tm1,phys->u_SfH_tm2, phys->u_SfHv_t, phys->u_SfHv_tm1,phys->uc,phys->uold2,phys->boundary_u,BCtop,BCbot,phys->Cn_T,
			0,0,NULL,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,1,1,prop->TVDmomentum);
	  UpdateScalars(grid,phys,prop,vert->omega_im,phys->v_SfH_tm1,phys->v_SfH_tm2, phys->v_SfHv_t, phys->v_SfHv_tm1,phys->vc,phys->vold2,phys->boundary_v,BCtop,BCbot,phys->Cn_T,
			0,0,NULL,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,1,1,prop->TVDmomentum);
	}
	ISendRecvCellData3D(phys->uc,grid,myproc,comm);
	ISendRecvCellData3D(phys->vc,grid,myproc,comm);

	REAL delta, delta_u, delta_v;
	for(j=0;j<grid->Ne;j++) {
	  nc1=grid->grad[2*j];
	  nc2=grid->grad[2*j+1];

	  if(grid->mark[j]==0 || grid->mark[j]==5) {
	    for(k=grid->etop[j];k<grid->Nke[j];k++) {
	      // This method does not take into account changes in h
	      /*
	      delta = 0.5*(phys->uc[nc1][k]+phys->uc[nc2][k]-phys->uold[nc1][k]-phys->uold[nc2][k])*grid->n1[j]
		+0.5*(phys->vc[nc1][k]+phys->vc[nc2][k]-phys->vold[nc1][k]-phys->vold[nc2][k])*grid->n2[j];
	      phys->u[j][k]+=delta;
	      */
	      // This method takes into account height changes for advection but not the old u
	      delta_u = (grid->dzz[nc1][k]*(phys->uc[nc1][k]-phys->uold[nc1][k])+
			  grid->dzz[nc2][k]*(phys->uc[nc2][k]-phys->uold[nc2][k]))/
		(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
	      delta_v = (grid->dzz[nc1][k]*(phys->vc[nc1][k]-phys->vold[nc1][k])+
			  grid->dzz[nc2][k]*(phys->vc[nc2][k]-phys->vold[nc2][k]))/
		(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
	      phys->u[j][k]+=delta_u*grid->n1[j]+delta_v*grid->n2[j];
	    }
	  }
	}
      }
      
      if(U_MOMENTUM_ADVECTION==1 && !prop->nonlinear) {
	int i, iptr, j, jptr, k, ne, nf;
	REAL u_im, fac1, fac2, fac3, u0, v0;
	REAL def1, def2;
	REAL div, divmax, err_max;
	REAL *J_eff_new = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"Debugging phys.c");
	REAL *J_eff_old = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"Debugging phys.c");
  
	REAL *a = phys->a;
	REAL *b = phys->b;
	REAL *c = phys->c;
	REAL *d = phys->d;
  REAL *e = phys->e;
  REAL *r = phys->pent_source;
  REAL *omega_eff = phys->pent_a;

	REAL *a0 = phys->pent_b;
	REAL *b0 = phys->bp;
	REAL *c0 = phys->bm;  		
  REAL *d0 = phys->pent_d;
  REAL *e0 = phys->pent_e; 
  REAL *ap = phys->ap;
  REAL *am = phys->am;
	
	fac1=prop->imfac1;
	fac2=prop->imfac2;
	fac3=prop->imfac3;	

	// Store u in utmp
	for(j=0;j<grid->Ne;j++) {
	  for(k=0;k<grid->Nke[j];k++) {
	    phys->utmp[j][k]=phys->u[j][k];
	    if(phys->u[j][k]!=phys->u[j][k]){
	       printf("u nan pre adv for j=%d, k=%d \n", j, k);
	    }
	    if(phys->u_old[j][k]!=phys->u_old[j][k]){
	      printf("u_old nan pre adv for j=%d, k=%d \n", j, k);
            }
	    if(phys->u_old2[j][k]!=phys->u_old2[j][k]){
	      printf("u_old2 nan pre adv for j=%d, k=%d \n", j, k);
            }
	    if(phys->u_ex[j][k]!=phys->u_ex[j][k]){
	      printf("u_ex nan pre adv for j=%d, k=%d \n", j, k);
            }
	     
	  }
	}
	  
	// Compute cell-centered horizontal advection:
	//    stmp = del dot (J U u) 
	//    stmp2 = del dot (J U v) 	
	for(i=0;i<grid->Nc;i++) {
	  for(k=0;k<grid->Nk[i];k++) {
	    phys->stmp[i][k]=phys->stmp2[i][k]=0;
	  }
	}

	//printf("Debugging - Main time loop: Advect u mom...\n");
	// Compute uc_ex and vc_ex
	ComputeUCPerot(phys->u_ex,phys->uold,phys->vold,phys->h,prop->kinterp,prop->subgrid,grid);
	ISendRecvCellData3D(phys->uold,grid,myproc,comm);
	ISendRecvCellData3D(phys->vold,grid,myproc,comm);

        //printf("Debugging - Main time loop: Advect u mom send and recv uold, vold from u->ex...\n");
	// Set U to a constant value to check CWC.
	if(prop->n==prop->nstart+prop->nsteps && CHECKCONSISTENCY) {
	  u0 = 0.01;
	  v0 = 0.0;
	  for(i=0;i<grid->Nc;i++) {
	    for(k=0;k<grid->Nk[i];k++) {
	      phys->uold[i][k]=u0;
	      phys->vold[i][k]=v0;	      
	    }
	  }
	  for(j=0;j<grid->Ne;j++) {
	    for(k=0;k<grid->Nke[j];k++) {
	      phys->utmp[j][k]=phys->u_ex[j][k]=u0*grid->n1[j]+v0*grid->n2[j];
	    }
	  }
	}
	
	/*
	 *
	 * Compute stmp = del dot (J U u) and stmp2 = del dot (J U v) 	 
	 *
	 */
	int nonlinear=5;
	int TVDmom = prop->TVDmomentum; //use for debugging
	int component, dimensions=2;
	REAL **stmp_pointer;
	// Loop through to compute horizontal advection of u and v to avoid repeated code
	for(component=0;component<dimensions;component++) {
	  // Interpolate the cell-centered velocity components onto the face and store
	  // in phys->ut using u^(im*)
	  if(component==0) {
	    GetMomentumFaceValues(phys->ut,phys->uold,phys->boundary_u,phys->u,
				  grid,phys,prop,comm,myproc,nonlinear,TVDmom);
	  } else {
	    GetMomentumFaceValues(phys->ut,phys->vold,phys->boundary_v,phys->u,
				  grid,phys,prop,comm,myproc,nonlinear,TVDmom);
	  }

	  if(prop->n==prop->nstart+prop->nsteps && CHECKCONSISTENCY) {
	    for(j=0;j<grid->Ne;j++) {
	      for(k=0;k<grid->Nke[j];k++) {
		if(component==0) {
		  phys->ut[j][k]=u0;
		} else {
		  phys->ut[j][k]=v0;
		}
	      }
	    }
	  }

	  if(component==0) {
	    stmp_pointer=phys->stmp;
	  } else {
	    stmp_pointer=phys->stmp2;
	  }
	  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	    i = grid->cellp[iptr];
	    
	    for(nf=0;nf<grid->nfaces[i];nf++) {
	      ne = grid->face[i*grid->maxfaces+nf];
	      
	      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
		
		u_im = fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k];

		stmp_pointer[i][k]+=1.0/grid->Ac[i]*
		  phys->ut[ne][k]*u_im*grid->dzf[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];
	      }
	    }
	  }
	  
	  ISendRecvCellData3D(stmp_pointer,grid,myproc,comm);
	  //printf("Debugging - Main time loop: Advect u mom send and recv stmp, for component %d...\n", component);
	}

	// Add explicit horizontal and vertical advection to utmp
	for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	  j = grid->edgep[jptr];
	  
	  nc1 = grid->grad[2*j];
	  nc2 = grid->grad[2*j+1];
	  Return_def(&def1,&def2,nc1,nc2,j,grid);
	  
	  for(k=grid->etop[j];k<grid->Nke[j];k++) {
	    J_eff_new[k]=(def1*grid->dzz[nc1][k]+def2*grid->dzz[nc2][k])/grid->dg[j];
	    J_eff_old[k]=(def1*grid->dzzold[nc1][k]+def2*grid->dzzold[nc2][k])/grid->dg[j];
	    if(J_eff_new[k]!=J_eff_new[k])
	      printf("J_eff_new[%d] nan for j %d, nc1=%d, nc2=%d \n", k, j, nc1, nc2);
	    if(J_eff_old[k]!=J_eff_old[k])
              printf("J_eff_old[%d] nan for j %d, nc1=%d, nc2=%d \n", k, j, nc1, nc2);
	  
	    if(J_eff_old[k]==0)
              printf("J_eff_old[%d] 0 for j %d, nc1=%d, nc2=%d \n", k, j, nc1, nc2);
	    if(J_eff_new[k]==0)
              printf("J_eff_new[%d] 0 for j %d, nc1=%d, nc2=%d \n", k, j, nc1, nc2);
	  }

	  //
	  // J_new u_new = J_old u_old - dt Adv(u) + J_new dt S
	  //
	  // S = (utmp - u_old)/dt
    // (u* - u_old)/dt???
	  //
	  // u_new = utmp + (J_old/J_new - 1)*u_old - dt Adv(u)
	  //
	  if(prop->n==prop->nstart+prop->nsteps && CHECKCONSISTENCY) {
	    // Only use this update to check for CWC since it is of the form
	    // J_new u_new = J_old/J_new u - dt Adv(u) which is not correct,
	    // although it does not require that we set u_old = U0
	    for(k=grid->etop[j];k<grid->Nke[j];k++) {
	      phys->utmp[j][k]=J_eff_old[k]/J_eff_new[k]*phys->utmp[j][k]-
		prop->dt/(J_eff_new[k]*grid->dg[j])*
		(grid->n1[j]*(def1*phys->stmp[nc1][k]+def2*phys->stmp[nc2][k])+    
		 grid->n2[j]*(def1*phys->stmp2[nc1][k]+def2*phys->stmp2[nc2][k]));
	    }
	  } else {
	    // This is the correct version although it will not be CWC because
	    // it requires that we set u_old = U0 to check for CWC which will violate
	    // the continuity equation (i.e. div \ne 0).
	    /* d/dt(J u_i) = J^{alpha} S_i */
	    REAL J_alpha;
	    for(k=grid->etop[j];k<grid->Nke[j];k++) {
	      //J_alpha = 0.5*(J_eff_old[k]+J_eff_new[k]); // Unstable with implicit or explicit
	      //J_alpha = J_eff_new[k];// Unstable with implicit or explicit
	      J_alpha = J_eff_old[k];// Stable with implicit or explicit

        //try at one half here too
	      phys->utmp[j][k]=J_alpha/J_eff_new[k]*phys->utmp[j][k]+
		(J_eff_old[k]-J_alpha)/J_eff_new[k]*phys->u_old[j][k]-
		prop->dt/(J_eff_new[k]*grid->dg[j])*
		(grid->n1[j]*(def1*phys->stmp[nc1][k]+def2*phys->stmp[nc2][k])+    
		 grid->n2[j]*(def1*phys->stmp2[nc1][k]+def2*phys->stmp2[nc2][k]));

	      if(phys->utmp[j][k]!=phys->utmp[j][k]){
		printf("utmp[%d][%d] is nan mid way through adv \n", j, k);
	      }
    	//   phys->utmp[j][k]=J_eff_new[k]/J_eff_new[k]*phys->utmp[j][k]-
      // prop->dt/(J_eff_new[k]*grid->dg[j])*
      // (grid->n1[j]*(def1*phys->stmp[nc1][k]+def2*phys->stmp[nc2][k])+    
      // grid->n2[j]*(def1*phys->stmp2[nc1][k]+def2*phys->stmp2[nc2][k]));
	    }
	  }

	  // omega interpolated to u edges
	  for(k=grid->etop[j];k<grid->Nke[j];k++) {
	    omega_eff[k] = (def1*vert->omega_im[nc1][k]+def2*vert->omega_im[nc2][k])/grid->dg[j];
	    if(omega_eff[k]!=omega_eff[k])
              printf("omega_eff[%d] nan for j %d, nc1=%d, nc2=%d \n", k, j, nc1, nc2);
	  }
    omega_eff[grid->Nke[j]]=0;

    if(nonlinear==2){
      // Tridiagonal coefficients for vertical advection with central differencing
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
        a0[k] = -0.5/J_eff_new[k]*omega_eff[k];
        b0[k] = -0.5/J_eff_new[k]*(omega_eff[k]-omega_eff[k+1]);
        c0[k] = 0.5/J_eff_new[k]*omega_eff[k+1];

      }

      // Top boundary is free-slip, i.e. u(k=ctop-1)=u(k=ctop)	  
      k=grid->etop[j];
      b0[k]+=a0[k];

      // Bottom boundary is no-flux, i.e. omega_eff[grid->Nke[j]]=0
      k=grid->Nke[j]-1;
      b0[k]=a0[k];
      //try BCBOT=2? was 1/
    } else if(nonlinear==5 && TVDmom>=5){
          //Pentadiagonal coefficients for QUICK/SHARP
      GetPentDiagNewAdv(a0, b0, c0, d0, e0, phys->wp,phys->wm,
        omega_eff,grid->dzz,phys->u_ex,j,grid->Nke[j],grid->etop[j],prop->dt,TVDmom,1,1);
    } else {
      //default to using TVD scheme 
      //need to make new version on edges!!!!
      GetApAmNewAdv(ap, am, phys->wp, phys->wm, phys->Cp,phys->Cm,phys->rp,phys->rm,
        omega_eff, grid->dzz,phys->u_ex, j,grid->Nke[j],grid->etop[j],prop->dt,TVDmom,1,1); 

      for(k=grid->etop[j];k<grid->Nke[j];k++) {
        //opposite negatives to previous code bc we add later
        a0[k] = -1/J_eff_new[k]*am[k];
        b0[k] = -1/J_eff_new[k]*(ap[k]-am[k+1]);
        c0[k] = 1/J_eff_new[k]*ap[k+1];
      }

      //not sure if we want these
      // Top boundary is free-slip, i.e. u(k=ctop-1)=u(k=ctop)	  
      k=grid->etop[j];
      b0[k]+=a0[k];

      // Bottom boundary is no-flux, i.e. omega_eff[grid->Nke[j]]=0
      //this was commented out..
      // k=grid->Nke[j]-1;
      // b0[k]=a0[k];
      // c0[k]=0;
    }

	  int EXPLICIT=0;
	  if(EXPLICIT) {
      if(nonlinear==2 || TVDmom<5){
        //central differencing
        for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++) {
          phys->utmp[j][k]+=prop->dt*(a0[k]*phys->u_ex[j][k-1]+
              b0[k]*phys->u_ex[j][k]+
              c0[k]*phys->u_ex[j][k+1]);
        }
        k=grid->etop[j];
        phys->utmp[j][k]+=prop->dt*(b0[k]*phys->u_ex[j][k]+
            c0[k]*phys->u_ex[j][k+1]);	  
        
        k=grid->Nke[j]-1;
        phys->utmp[j][k]+=prop->dt*(a0[k]*phys->u_ex[j][k-1]+
            b0[k]*phys->u_ex[j][k]);	  
      } else if(nonlinear==5 && TVDmom>=5){
        //QUICK/SHARP
        for(k=grid->etop[j]+2;k<grid->Nke[j]-2;k++) {
          phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(a0[k]*phys->u_ex[j][k-2]+
              b0[k]*phys->u_ex[j][k-1]+
              c0[k]*phys->u_ex[j][k] + d0[k]*phys->u_ex[j][k+1] + 
              e0[k]*phys->u_ex[j][k+2]);

        }

      //had commented all these out, lets test back in
      k=grid->etop[j]+1;
      phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(b0[k]*phys->u_ex[j][k-1]+
          c0[k]*phys->u_ex[j][k]+d0[k]*phys->u_ex[j][k+1]+
          e0[k]*phys->u_ex[j][k+2]);	 

      k=grid->etop[j];
      phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(c0[k]*phys->u_ex[j][k]+
          d0[k]*phys->u_ex[j][k+1]+e0[k]*phys->u_ex[j][k+2]);	  

      k=grid->Nke[j]-2;
      phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(a0[k]*phys->u_ex[j][k-2]+
          b0[k]*phys->u_ex[j][k-1]+c0[k]*phys->u_ex[j][k]+
          d0[k]*phys->u_ex[j][k+1]);	
      k=grid->Nke[j]-1;
      phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(a0[k]*phys->u_ex[j][k-2]+
          b0[k]*phys->u_ex[j][k-1]+c0[k]*phys->u_ex[j][k]);	
      }
	  } else {
	    REAL ifac1 = fac1, ifac2 = fac2, ifac3 = fac3;
	    // For debugging to test different implicit schemes
	    //  ifac1=1;
	    //	ifac2=0;
	    //	ifac3=0;

	    if(prop->n==prop->nstart+prop->nsteps && CHECKCONSISTENCY) {
	      if(myproc==0 && jptr==grid->edgedist[0]) {
		printf("********************************************************************\n");
		printf("* Warning!!! Setting u_old=u_old2=U0 --> Div check will not be zero!\n");
		printf("********************************************************************\n");
	      }
	      for(k=grid->etop[j];k<grid->Nke[j];k++) {
		phys->u_old[j][k]=phys->u_old2[j][k]=u0*grid->n1[j]+v0*grid->n2[j];
	      }
	    }
	    
	    // Right-hand side of implicit vertical advection
      if(nonlinear==2 || TVDmom<5){
        for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++) {
            phys->utmp[j][k]+=prop->dt*(a0[k]*(ifac2*phys->u_old[j][k-1]+ifac3*phys->u_old2[j][k-1])+
                 b0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k])+
               c0[k]*(ifac2*phys->u_old[j][k+1]+ifac3*phys->u_old2[j][k+1]));
          }

          k=grid->etop[j];
          phys->utmp[j][k]+=prop->dt*(b0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k])+
          c0[k]*(ifac2*phys->u_old[j][k+1]+ifac3*phys->u_old2[j][k+1]));
      
           k=grid->Nke[j]-1;
           phys->utmp[j][k]+=prop->dt*(a0[k]*(ifac2*phys->u_old[j][k-1]+ifac3*phys->u_old2[j][k-1])+
               b0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k]));
      } else{
	//try flipping sign here and for bcs???
          for(k=grid->etop[j]+2;k<grid->Nke[j]-2;k++) {
              phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(a0[k]*(ifac2*phys->u_old[j][k-2]+ifac3*phys->u_old2[j][k-2])+
                b0[k]*(ifac2*phys->u_old[j][k-1]+ifac3*phys->u_old2[j][k-1])+
                c0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k])+
                d0[k]*(ifac2*phys->u_old[j][k+1]+ifac3*phys->u_old2[j][k+1])+
	        e0[k]*(ifac2*phys->u_old[j][k+2]+ifac3*phys->u_old2[j][k+2]));
	      if(phys->utmp[j][k]!=phys->utmp[j][k]){
		printf("utmp[%d][%d] is nan pre bcs \n", j, k);
	      }
	  }

	    // k=grid->etop[j];
	    // phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(b0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k])+
			// 		c0[k]*(ifac2*phys->u_old[j][k+1]+ifac3*phys->u_old2[j][k+1]));
	    
	    // k=grid->Nke[j]-1;
	    // phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(a0[k]*(ifac2*phys->u_old[j][k-1]+ifac3*phys->u_old2[j][k-1])+
			// 		b0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k]));

      //again, had commented these all out so lets re test
      if(grid->Nke[j]>1){
        k=grid->etop[j]+1;
        phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(b0[k]*(ifac2*phys->u_old[j][k-1]+ifac3*phys->u_old2[j][k-1])+
            c0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k])+d0[k]*(ifac2*phys->u_old[j][k+1]+ifac3*phys->u_old2[j][k+1])+
            e0[k]*(ifac2*phys->u_old[j][k+2]+ifac3*phys->u_old2[j][k+2]));	 

	if(phys->utmp[j][k]!=phys->utmp[j][k]){
	  printf("utmp[%d][%d] is nan post bc \n", j, k);
	}

        k=grid->etop[j];
        phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(c0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k])+
          d0[k]*(ifac2*phys->u_old[j][k+1]+ifac3*phys->u_old2[j][k+1])+e0[k]*(ifac2*phys->u_old[j][k+2]+ifac3*phys->u_old2[j][k+2]));	  

	if(phys->utmp[j][k]!=phys->utmp[j][k]){
          printf("utmp[%d][%d] is nan post bc \n", j, k);
        }
	
        k=grid->Nke[j]-2;
        phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(a0[k]*(ifac2*phys->u_old[j][k-2]+ifac3*phys->u_old2[j][k-2])+
            b0[k]*(ifac2*phys->u_old[j][k-1]+ifac3*phys->u_old2[j][k-1])+c0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k])+
            d0[k]*(ifac2*phys->u_old[j][k+1]+ifac3*phys->u_old2[j][k+1]));	

	if(phys->utmp[j][k]!=phys->utmp[j][k]){
          printf("utmp[%d][%d] is nan post bc \n", j, k);
        }

        k=grid->Nke[j]-1;
        phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(a0[k]*(ifac2*phys->u_old[j][k-2]+ifac3*phys->u_old2[j][k-2])+
            b0[k]*(ifac2*phys->u_old[j][k-1]+ifac3*phys->u_old2[j][k-1])+c0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k]));	
	if(phys->utmp[j][k]!=phys->utmp[j][k]){
          printf("utmp[%d][%d] is nan post bc \n", j, k);
        }  
      }else{

        k=grid->etop[j];
        phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(c0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k])+
          d0[k]*(ifac2*phys->u_old[j][k+1]+ifac3*phys->u_old2[j][k+1])+e0[k]*(ifac2*phys->u_old[j][k+2]+ifac3*phys->u_old2[j][k+2]));	  


        k=grid->Nke[j]-1;
        phys->utmp[j][k]+=prop->dt/J_eff_new[k]*(a0[k]*(ifac2*phys->u_old[j][k-2]+ifac3*phys->u_old2[j][k-2])+
            b0[k]*(ifac2*phys->u_old[j][k-1]+ifac3*phys->u_old2[j][k-1])+c0[k]*(ifac2*phys->u_old[j][k]+ifac3*phys->u_old2[j][k]));	
        }
      }


	    // Create coefficients for the tridiagonal
      if(nonlinear==2 || TVDmom<5){ 
        for(k=grid->etop[j];k<grid->Nke[j];k++) {
          a[k]=-prop->dt*ifac1*a0[k];
          b[k]=1.0-prop->dt*ifac1*b0[k];
          c[k]=-prop->dt*ifac1*c0[k];
          d[k]=phys->utmp[j][k];
        }
      } else{
        for(k=grid->etop[j];k<grid->Nke[j];k++) {
          a[k]=-prop->dt/J_eff_new[k]*ifac1*a0[k];
          b[k]=-prop->dt/J_eff_new[k]*ifac1*b0[k];
          c[k]=1.0-prop->dt/J_eff_new[k]*ifac1*c0[k];
          d[k]=-prop->dt/J_eff_new[k]*ifac1*d0[k];
          e[k]=-prop->dt/J_eff_new[k]*ifac1*e0[k];
          r[k]=phys->utmp[j][k];

	  //printf("for j=%d, k=%d: a0=%f, b0=%f, c0=%f, d0=%f, e0=%f  on proc %d \n", j, k, a0[k], b0[k], c0[k], d0[k], e0[k], myproc);
          //printf("for j=%d, k=%d: a=%f, b=%f, c=%f, d=%f, e=%f, r=%f on proc %d \n", j, k, a[k], b[k], c[k], d[k], e[k], r[k], myproc);

	}
      }

      for(k=grid->etop[j];k<grid->Nke[j];k++){
	if(phys->utmp[j][k]!=phys->utmp[j][k]){
	  printf("utmp[%d][%d] is nan pre penta through adv \n", j, k);
	}
      }
      if(nonlinear==2 || TVDmom<5){
	    TriSolve(&(a[grid->etop[j]]),&(b[grid->etop[j]]),&(c[grid->etop[j]]),
		     &(d[grid->etop[j]]),&(phys->utmp[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);
      }else{
        //printf("u for j=%d, k=0 is %f before reduce penta \n", j, phys->utmp[j][0]);
        ReducePentadiag(&(a[grid->etop[j]]),&(b[grid->etop[j]]),&(c[grid->etop[j]]),
        &(d[grid->etop[j]]),&(e[grid->etop[j]]),&(r[grid->etop[j]]),&(phys->utmp[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);
        //printf("u for j=%d, k=0 is %f after reduce penta \n", j, phys->utmp[j][0]);
      }
      for(k=grid->etop[j];k<grid->Nke[j];k++){
	if(phys->utmp[j][k]!=phys->utmp[j][k]){
	  printf("utmp[%d][%d] is nan post penta through adv on proc %d. grid->xe[j]=%f, grid->ye[j]=%f. \n", j, k, myproc, grid->xe[j], grid->ye[j]);
	}
      }
	  }
	}
        //printf("Debugging - Main time loop: Advect u mom through matrix inversion...\n");	  

	if(prop->n==prop->nstart+prop->nsteps && CHECKCONSISTENCY) {
	  divmax = -INFTY;
	  err_max = -INFTY;
	  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	    i = grid->cellp[iptr];
	    
	    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	      div=vert->omega_im[i][k]-vert->omega_im[i][k+1];
	      for(nf=0;nf<grid->nfaces[i];nf++) {
		ne = grid->face[i*grid->maxfaces+nf];
		u_im = fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k];
		
		div+=1.0/grid->Ac[i]*
		  u_im*grid->dzf[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];
	      }
	      div=(grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt+div;
	      if(fabs(div)>divmax)
		divmax=fabs(div);

	      /*
	      if(myproc==0 && i==floor(grid->Nc/2)) {
		printf("k=%d, div=%.4e\n",k,div);
	      }
	      */
	    }
	  }

	  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	    j = grid->edgep[jptr];
	    
	    REAL U0 = u0*grid->n1[j]+v0*grid->n2[j];
	    REAL err_max_j = -INFTY;
	    for(k=grid->etop[j];k<grid->Nke[j];k++) {
	      if(fabs((phys->utmp[j][k]-U0)/U0)>err_max_j) {
		err_max_j = fabs((phys->utmp[j][k]-U0)/U0);
	      }
	    }

	    if(err_max_j>err_max) err_max = err_max_j;
	    /*
	    if(err_max_j>1e-6) {
	      k=0;
	      printf("k = %d, U = %.4e, U0 = %.2e, |(U-U0)/U0| = %.4e, mark=%d\n",
		     k,phys->utmp[j][k],U0,fabs((phys->utmp[j][k]-U0)/U0),grid->mark[j]);
	    }
	    */
	  }
	  printf("U Adv: proc: %d, divmax = %.4e, err_max = %.4e\n",myproc,divmax,err_max);
	}

	/*for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	  j = grid->edgep[jptr];	

	  for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->u[j][k]=phys->utmp[j][k];
    }
    } */

	//commented out above so this actuall works lol 11/9/25
  REAL tmp;
	for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
	  j = grid->edgep[jptr];	

	  for(k=grid->etop[j];k<grid->Nke[j];k++) {
	    tmp = phys->u[j][k];
	    phys->u[j][k]=phys->utmp[j][k];
	    phys->utmp[j][k] = tmp;
	    //undo
	    phys->utmp[j][k] = phys->u[j][k];
	  }
	}
        //printf("Debugging - Main time loop: set u to utmp then back to u...\n");

	//make sure boundaries are back to whatever u is 
	for(jptr=grid->edgedist[1];jptr<grid->edgedist[5];jptr++){
	  j = grid->edgep[jptr];
	   for(k=grid->etop[j];k<grid->Nke[j];k++) {
	     phys->utmp[j][k] = phys->u[j][k];
	   }
	}
        //printf("Debugging - Main time loop: Advect utmp on mom bounds back to u...\n");

  for(j=0;j<grid->Ne;j++){
    for(k=grid->etop[j];k<grid->Nke[j];k++){
      if(phys->u[j][k]!=phys->u[j][k]){
	printf("u nan for j=%d, k=%d \n", j, k);
      }
      if(phys->utmp[j][k]!=phys->utmp[j][k]){
        printf("utmp nan for j=%d, k=%d \n", j, k);
      }
    }
  }
  //printf("Debugging - Main time loop: Advect utmp nan check...\n");

	SunFree(J_eff_new,grid->Nkmax*sizeof(REAL),"U_MOMENTUM_ADVECTION==1");
	SunFree(J_eff_old,grid->Nkmax*sizeof(REAL),"U_MOMENTUM_ADVECTION==1");

	//printf("Debugging - Main time loop: freed Js...\n");

      }
      //or, try set flux height here THEN layer av cont
      
      ISendRecvEdgeData3D(phys->u,grid,myproc,comm); //send new u after adv

      if(DEBUG) printf("Debugging - Main time loop: Update omega after u adv...\n");
      // Need to update continuity after adding momentum advection to u
      if(prop->vertcoord==1) {
	Continuity(phys->wnew,grid,phys,prop);	
	ISendRecvWData(phys->wnew,grid,myproc,comm);
      } else {
	//SetFluxHeight(grid,phys,prop,vert->dzfmeth,comm,myproc); //added set flux height
	LayerAveragedContinuity(vert->omega,grid,prop,phys,comm,myproc);
	ISendRecvWData(vert->omega,grid,myproc,comm);
	if(DEBUG) printf("Debugging - Main time loop: send and recieved omega...\n");
      }
     //move after Wpredictor
      
      // update subgrid->Acceff and subgrid->Acveff for non hydrostatic
      if(prop->subgrid) {
	if(DEBUG) printf("Debugging - Main time loop: UpdateSubgridVerticalAceff...\n");				
        UpdateSubgridVerticalAceff(grid, phys, prop, 1, myproc);
      }
      // Compute vertical momentum and the nonhydrostatic pressure
      if(prop->nonhydrostatic && !blowup) {
        // Predicted vertical velocity field is in phys->w
	if(DEBUG) printf("Debugging - Main time loop: WPredictor...\n");
  //here, set u back to utmp 
	WPredictor(grid,phys,prop,myproc,numprocs,comm);
        ISendRecvWData(phys->w,grid,myproc,comm);

  //here, 

  // Need to update continuity after adding momentum advection to u
  // 	SetFluxHeight(grid,phys,prop,vert->dzfmeth,comm,myproc);

  //     if(prop->vertcoord==1) {
	// Continuity(phys->wnew,grid,phys,prop);	
	// ISendRecvWData(phys->wnew,grid,myproc,comm);
  //     } else {
	// LayerAveragedContinuity(vert->omega,grid,prop,phys,comm,myproc);
	// ISendRecvWData(vert->omega,grid,myproc,comm);
  //     }

        // Wpredictor calculate w^*
        // now calculate omega^* for the new generalized vertical coordinate
        if(prop->vertcoord!=1){
          // recalculate uc and vc for the predictor velocity field
	  if(DEBUG) printf("Debugging - Main time loop: ComputeUC (WPredictor)...\n");	  
          ComputeUC(phys->uc, phys->vc, phys,grid, myproc, prop->interp,prop->kinterp,prop->subgrid);
          // now we have uc^* and vc^*
          // compute ul^* and vl^* at each layer top
	  if(DEBUG) printf("Debugging - Main time loop: ComputeUl (WPredictor)...\n");	  	  
          ComputeUl(grid, prop, phys, myproc);
          // compute zc gradient
          // zf is calculate in ComputeZc function
	  if(DEBUG) printf("Debugging - Main time loop: ComputeCellAveragedHorizontalGradient (WPredictor)...\n");
	  //	  ComputeCellAveragedHorizontalGradient(vert->dzdx, 0, vert->zf, grid, prop, phys, myproc);
	  //	  ComputeCellAveragedHorizontalGradient(vert->dzdy, 1, vert->zf, grid, prop, phys, myproc);
	  ComputeCellAveragedHorizontalGradientCell(vert->dzdx, 0, vert->zl, 1, grid, prop, phys, myproc);
	  ComputeCellAveragedHorizontalGradientCell(vert->dzdy, 1, vert->zl, 1, grid, prop, phys, myproc);
	  /*
	  int iptr;
	  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	    i=grid->cellp[iptr];
	    printf("%f %f\n",grid->xv[i],vert->dzdx[i][grid->Nk[i]]);
	  }
	  */
	  /*
	  for(j=0;j<grid->Ne;j++) {
	    if(grid->mark[j]==0)
	      printf("%f %f\n",grid->xe[j],vert->zf[j][grid->Nk[i]]);
	  }
	  */
	  
          // compute U3
	  if(DEBUG) printf("Debugging - Main time loop: ComputeOmega(WPredictor)...\n");	  
          ComputeOmega(grid, prop, phys,-1, myproc);
        }

        // Source term for the pressure-Poisson equation is in phys->stmp
	if(DEBUG) printf("Debugging - Main time loop: ComputeQSource...\n");
	/*
	for(j=0;j<grid->Ne;j++) {
	  nc1=grid->grad[2*j];
	  nc2=grid->grad[2*j+1];
	  if(nc1==-1) nc1=nc2;
	  if(nc2==-1) nc2=nc1;
	  for(k=0;k<grid->Nkmax;k++) {
	    grid->dzf[j][k]=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
	  }
	}
	*/
	//	TvdFluxHeight(grid, phys, prop, vert->dzfmeth,comm, myproc);
	// Set flux height based on u^{im} (previous was based on u^{im*}
	
	//comment out
	//SetFluxHeight(grid,phys,prop,vert->dzfmeth,comm,myproc);
	//ISendRecvEdgeData3D(grid->dzf,grid,myproc,comm);

	ComputeQSource(phys->stmp,grid,phys,prop,myproc,numprocs);

        // Solve for the nonhydrostatic pressure.  
        // phys->stmp2/qc contains the initial guess
        // phys->stmp contains the source term
        // phys->stmp3 is used for temporary storage
	if(DEBUG) printf("Debugging - Main time loop: CGSolveQ...\n");

	if(prop->include_slope_terms)
	  BiCGSolveQ(phys->qc,phys->stmp,phys->stmp3,grid,phys,prop,myproc,numprocs,comm);
	else
	  CGSolveQ(phys->qc,phys->stmp,phys->stmp3,grid,phys,prop,myproc,numprocs,comm);
	
        // Correct the nonhydrostatic velocity field with the nonhydrostatic pressure
        // correction field phys->stmp2/qc.  This will correct phys->u so that it is now
        // the volume-conserving horizontal velocity field.  
        // phys->w is not corrected since
        // it is obtained via continuity.  
        // Also, update the total nonhydrostatic pressure
        // with the pressure correction.
	if(DEBUG) printf("Debugging - Main time loop: Corrector...\n");
        Corrector(phys->qc,grid,phys,prop,myproc,numprocs,comm);

        // Send/recv the horizontal velocity data after it has been corrected.
        ISendRecvEdgeData3D(phys->u,grid,myproc,comm);

        // Send q to the boundary cells now that it has been updated
        ISendRecvCellData3D(phys->q,grid,myproc,comm);
      } else if(!(prop->interp == PEROT)) {
        // support quadratic interpolation work
        // Send/recv the horizontal velocity data for use with more complex interpolation
	if(DEBUG) printf("Debugging - Main time loop: ISendRecvEdgeData3D (if(interp!=PEROT))...\n");
	ISendRecvEdgeData3D(phys->u,grid,myproc,comm);
      }
      t_nonhydro+=Timer()-t0;

      // apply continuity via Eqn 82
      int old_nonhydrostatic=prop->nonhydrostatic;
      //      prop->nonhydrostatic=0;
      
  //     if(prop->vertcoord==1)
  //     {
  //       // w_im is calculated
	// if(DEBUG) printf("Debugging - Main time loop: Continuity (after q)...\n");
  //       Continuity(phys->w,grid,phys,prop);
  //       ISendRecvWData(phys->w,grid,myproc,comm);
  //     } else {
  //       // NOW vert->omega stores the omega^n+1 with nonhydrostatic pressure correction
  //       if(!prop->nonhydrostatic || prop->vertcoord==5)   
  //       {
	//   if(DEBUG) printf("Debugging - Main time loop: LayerAveragedContinuity...\n");	  
  //         LayerAveragedContinuity(vert->omega,grid,prop,phys,comm,myproc);
  //         ISendRecvWData(vert->omega,grid,myproc,comm);          
  //       }
  //     }
        /* //commented above back in...


        // if nonhydrostatic=1, w is solved in the corrector function
        // no need to solve from omega
        // w is recalculated by omega only for hydrostatic case
        if(!prop->nonhydrostatic)
        {  
          // recalculate uc and vc for the predictor velocity field
	  if(DEBUG) printf("Debugging - Main time loop: ComputeUC (hydrostatic)...\n");	  	  
          ComputeUC(phys->uc, phys->vc, phys,grid, myproc, prop->interp,prop->kinterp,prop->subgrid);
          // now we have uc^* and vc^*
          // compute ul^* and vl^* at each layer top
	  if(DEBUG) printf("Debugging - Main time loop: ComputeUl (hydrostatic)...\n");	  	  	  
          ComputeUl(grid, prop, phys, myproc);
          // compute zc gradient
          // zf is calculate in ComputeZc function
	  if(DEBUG) printf("Debugging - Main time loop: ComputcellAveragedHorizontalGradient (hydrostatic)...\n");
	  ComputeCellAveragedHorizontalGradientCell(vert->dzdx, 0, vert->zl, 1, grid, prop, phys, myproc);
	  ComputeCellAveragedHorizontalGradientCell(vert->dzdy, 1, vert->zl, 1, grid, prop, phys, myproc); 	  
	  //	  ComputeCellAveragedHorizontalGradient(vert->dzdx, 0, vert->zf, grid, prop, phys, myproc);
	  //          ComputeCellAveragedHorizontalGradient(vert->dzdy, 1, vert->zf, grid, prop, phys, myproc); 

        }

        if(!prop->nonhydrostatic || prop->vertcoord==5)   {
          // compute w from omega
	  if(DEBUG) printf("Debugging - Main time loop: ComputeOmega (hydrostatic)...\n");	  
          ComputeOmega(grid, prop, phys,0, myproc);
	}

        ISendRecvWData(phys->w,grid,myproc,comm);

        // update U3 with the new w
	if(DEBUG) printf("Debugging - Main time loop: ComputeU3...\n");	  	
        ComputeOmega(grid, prop, phys,-1, myproc);
        ISendRecvWData(vert->U3,grid,myproc,comm);
      }
      prop->nonhydrostatic=old_nonhydrostatic;
      */
      ISendRecvWData(phys->w,grid,myproc,comm);         
      
      // Set scalar and wind stress boundary values at new time step n; 
      // (n-1) is old time step.
      // BoundaryVelocities and OpenBoundaryFluxes were called in UPredictor to set the
      // boundary velocities to the new time step values for use in the 
      // free surface calculation.
      if(DEBUG) printf("Debugging - Main time loop: BoundaryScalars...\n");	        
      BoundaryScalars(grid,phys,prop,myproc,comm);

      if(DEBUG) printf("Debugging - Main time loop: WindStress...\n");	              
      WindStress(grid,phys,prop,met,myproc);

      // dz change so hmarshleft and marshtop may change
      if(prop->marshmodel) {
      if(DEBUG) printf("Debugging - Main time loop: SetMarsh...\n");	              	
        SetMarshTop(grid,phys,myproc);
      }

      if(prop->beta || prop->gamma) {
	if(DEBUG) printf("Debugging - Main time loop: SetDensity...\n");	              		
        SetDensity(grid,phys,prop);
      }

      // calculate any variable values
      if(DEBUG) printf("Debugging - Main time loop: UserDefinedFunction...\n");	              	      
      UserDefinedFunction(grid,phys,prop,myproc);

      // u now contains velocity on all edges at the new time step
      if(DEBUG) printf("Debugging - Main time loop: ComputeUC...\n");	              	            
      ComputeUC(phys->uc, phys->vc, phys,grid, myproc, prop->interp,prop->kinterp,prop->subgrid);
      //printf("Done (%d).\n",myproc);

      // now send interprocessor data
      ISendRecvCellData3D(phys->uc,grid,myproc,comm);
      ISendRecvCellData3D(phys->vc,grid,myproc,comm);
    }

    // Adjust the velocity field in the new cells if the newcells variable is set 
    // to 1 in suntans.dat.  Once this is done, send the interprocessor 
    // u-velocities to the neighboring processors.
    if(prop->newcells) {
      if(DEBUG) printf("Debugging - Main time loop: NewCells...\n");	              	                  
      NewCells(grid,phys,prop);
      ISendRecvEdgeData3D(phys->u,grid,myproc,comm);
    }

    // Compute average
    if(prop->calcaverage){
      if(DEBUG) printf("Debugging - Main time loop: Averages...\n");
      UpdateAverageVariables(grid,average,phys,met,prop,comm,myproc); 
      UpdateAverageScalars(grid,average,phys,met,prop,comm,myproc); 
    }

    // Check whether or not run is blowing up
    t0=Timer();
    if(DEBUG) printf("Debugging - Main time loop: Check blowup...\n");
    blowup=(Check(grid,phys,prop,myproc,numprocs,comm) || blowup);
    t_check+=Timer()-t0;
    
    // Output data based on ntout specified in suntans.dat
    t0=Timer();
    if (prop->outputNetcdf==0){
      // Write to binary
      if(DEBUG) printf("Debugging - Main time loop: OutputPhysicalVariables...\n");
      OutputPhysicalVariables(grid,phys,prop,myproc,numprocs,blowup,comm); 
      // subgrid
      if(prop->subgrid) {
	if(DEBUG) printf("Debugging - Main time loop: OutputSubgridVariables...\n");	
        OutputSubgridVariables(grid, prop, myproc, numprocs, comm);
      }
      if(prop->vertcoord!=1) {
	if(DEBUG) printf("Debugging - Main time loop: OutputVertCoordinate...\n");		
        OutputVertCoordinate(grid,prop,myproc,numprocs,comm);
      }
    }else {
      // Output data to netcdf
      if(DEBUG) printf("Debugging - Main time loop: WriteOutputNC...\n");
      WriteOutputNC(prop, grid, phys, met, blowup, myproc);
    }

    // Output the average arrays
    if(prop->calcaverage){
      if(DEBUG) printf("Debugging - Main time loop: Output Averages...\n");      
      if(prop->mergeArrays){
        WriteAverageNCmerge(prop,grid,average,phys,met,blowup,numprocs,comm,myproc);
      }else{
        WriteAverageNC(prop,grid,average,phys,met,blowup,comm,myproc);
      }
    }
    InterpData(grid,phys,prop,comm,numprocs,myproc);

    t_io+=Timer()-t0;
    // Output progress
    if(DEBUG) printf("Debugging - Main time loop: Progress...\n");
    Progress(prop,myproc,numprocs);
    if(blowup)
      break;

    //Close all open netcdf file
    /*
    if(prop->n==prop->nsteps+prop->nstart) {
      if(prop->outputNetcdf==1){
        //printf("Closing output netcdf file on processor: %d\n",myproc);
      	MPI_NCClose(prop->outputNetcdfFileID);
      }
      if(prop->netcdfBdy==1){
        //printf("Closing boundary netcdf file on processor: %d\n",myproc);
      	MPI_NCClose(prop->netcdfBdyFileID);
      }
      if(prop->readinitialnc==1){
        //printf("Closing initial netcdf file on processor: %d\n",myproc);
      	MPI_NCClose(prop->initialNCfileID );
      }
      if(prop->metmodel>0){
        //printf("Closing met netcdf file on processor: %d\n",myproc);
      	MPI_NCClose(prop->metncid);
      }
    }
    */
  }

  if(prop->mergeArrays) {
    if(VERBOSE>2 && myproc==0) printf("Freeing merging arrays...\n");
    if(DEBUG) printf("Debugging - Main time loop: FreeMergingArrays...\n");    
    FreeMergingArrays(grid,myproc);
  }
}

/*
 * Function: StoreVariables
 * Usage: StoreVariables(grid,phys);
 * ---------------------------------
 * Store the old values of s, u, and w into stmp3, u_old, and w_old,
 * respectively.
 *
 */
static void StoreVariables(gridT *grid, physT *phys, propT *prop) {
  int i, j, k, iptr, jptr;
  REAL fab1, fab2, fab3;
  
 // Adams Bashforth coefficients
  if(prop->readOldVelocity){
  //if(1){
    fab1=prop->exfac1;
    fab2=prop->exfac2;
    fab3=prop->exfac3;

  } else if(prop->n==1 || prop->wetdry) {
    fab1=1;
    fab2=fab3=0;

    for(j=0;j<grid->Ne;j++)
      for(k=0;k<grid->Nke[j];k++)
        phys->Cn_U[j][k]=phys->Cn_U2[j][k]=0;
  } else if(prop->n==2) {
    fab1=3.0/2.0;
    fab2=-1.0/2.0;
    fab3=0;
  } else {
    fab1=prop->exfac1;
    fab2=prop->exfac2;
    fab3=prop->exfac3;   
  }
  
  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i]+1;k++) {
      phys->w_ex[i][k]=fab1*phys->w[i][k]+fab2*phys->w_old[i][k]+fab3*phys->w_old2[i][k];
      phys->stmp3[i][k]=phys->s[i][k];
      phys->w_old2[i][k]=phys->w_old[i][k];
      phys->w_old[i][k]=phys->w[i][k];      
    }

  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=0;
    k=0;
    for(k=0;k<grid->Nke[j];k++)
    {
      phys->u_ex[j][k]=fab1*phys->u[j][k]+fab2*phys->u_old[j][k]+fab3*phys->u_old2[j][k];
      phys->u_old2[j][k]=phys->u_old[j][k];
      phys->utmp[j][k]=phys->u_old[j][k]=phys->u[j][k];
    }
  }
}

/*
 * Function: HorizontalSourceTerms_n1
 * Usage: HorizontalSource(grid,phys,prop,myproc,numprocs);
 * --------------------------------------------------------
 * Compute the horizontal source term that is used to obtain the free surface.
 *
 * This function adds the following to the horizontal source term:
 *
 * 1) Old nonhydrostatic pressure gradient with theta method
 * 2) Coriolis terms with AB2
 * 3) Baroclinic term with AB2
 * 4) Horizontal and vertical advection of horizontal momentum with AB2
 * 5) Horizontal laminar+turbulent diffusion of horizontal momentum
 *
 * Cn_U contains the Adams-Bashforth terms at time step n-1.
 * If wetting and drying is employed, no advection is computed in 
 * the upper cell.
 *
 */
static void HorizontalSourceTerms_n1(gridT *grid, physT *phys, propT *prop,
  int myproc, int numprocs, MPI_Comm comm) 
{
int i, ib, iptr, boundary_index, nf, j, jptr, k, nc, nc1, nc2, ne, 
k0, kmin, kmax;
//REAL *a, *b, *c, sum, def1, def2, dgf, Cz, tempu,Ac,f_sum, ke1,ke2,Vm; //AB3
// additions to test divergence averaging for w in momentum calc
//REAL wedge[3], lambda[3], wik, tmp_x, tmp_y,tmp;
REAL *a, *b, *c, sum, def1, def2, dgf, Cz, tempu,Ac,f_sum, ke1,ke2,Vm; //AB3
int aneigh;

REAL **rho_t1, **rho_t2, **s_old2, **dzz_old2, **zc_old, **zc_old2, z;
//REAL **rho_t1, **rho_t2;

a = phys->a;
b = phys->b;
c = phys->c;

char str[BUFFERLENGTH], filename[BUFFERLENGTH];
FILE *fid;

//test one timestep forward.

  //allocate space for rho1, rho2, s_old2
  rho_t1 = (REAL **)malloc(grid->Nc*sizeof(REAL *));
  rho_t2 = (REAL **)malloc(grid->Nc*sizeof(REAL *));
  s_old2 = (REAL **)malloc(grid->Nc*sizeof(REAL *));
  dzz_old2 = (REAL **)malloc(grid->Nc*sizeof(REAL *));
  zc_old = (REAL **)malloc(grid->Nc*sizeof(REAL *));
  zc_old2 = (REAL **)malloc(grid->Nc*sizeof(REAL *));
  for(i=0;i<grid->Nc;i++) {
    rho_t1[i] = (REAL *)malloc(grid->Nk[i]*sizeof(REAL));
    rho_t2[i] = (REAL *)malloc(grid->Nk[i]*sizeof(REAL));
    s_old2[i] = (REAL *)malloc(grid->Nk[i]*sizeof(REAL));
    dzz_old2[i] = (REAL *)malloc(grid->Nk[i]*sizeof(REAL));
    zc_old[i] = (REAL *)malloc(grid->Nk[i]*sizeof(REAL));
    zc_old2[i] = (REAL *)malloc(grid->Nk[i]*sizeof(REAL));
  } 

  if((int)MPI_GetValue(DATAFILE,"init_st_from_file","ReadGrid",myproc)){
    MPI_GetFile(filename,DATAFILE,"s_t2_init_file","HorizontalSource",myproc);  
    sprintf(str,"%s.%d",filename,myproc);  
    fid = MPI_FOpen(str,"r","HorizontalSource",myproc);

    for(j=0;j<grid->Nc;j++) {
      fread(s_old2[j],sizeof(REAL),grid->Nkmax,fid);
    }
    fclose(fid);
  } else{
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        s_old2[i][k] = phys->s[i][k];
      }
    }
  } 

  //add in to read in dzz_old2
  if(prop->vertcoord==4 && prop->readOldVelocity==1){
    MPI_GetFile(filename,DATAFILE,"dzz_t2_init_file","HorizontalSource",myproc);  
    sprintf(str,"%s.%d",filename,myproc);  
    fid = MPI_FOpen(str,"r","HorizontalSource",myproc);

    for(j=0;j<grid->Nc;j++) {
      fread(dzz_old2[j],sizeof(REAL),grid->Nkmax,fid);
    }
    fclose(fid);
  } else{
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        dzz_old2[i][k] = grid->dzzold[i][k];
      }
    }
  }

//initialize source terms as 0
for(j=0;j<grid->Ne;j++){
  for(k=0;k<grid->Nke[j];k++){
    phys->Cn_U[j][k]=phys->Cn_U2[j][k]=0;
  }
}

//calculate densities at previos time steps
for(i=0;i<grid->Nc;i++) {
for(k=grid->ctop[i];k<grid->Nk[i];k++) {
  REAL p =0; //our state eq doesnt use p, but should really be more general here 
  rho_t1[i][k]=StateEquation(prop,phys->s_old[i][k],phys->T[i][k],p);
  rho_t2[i][k]=StateEquation(prop,s_old2[i][k],phys->T[i][k],p);
}
}

//calculate zcs
for(i=0;i<grid->Nc;i++)
{
  //printf("zb[i] = %f, zBold[i] = %f, zBold2[i] = %f for i = %d \n", phys->zB[i], phys->zBold[i], phys->zBold2[i], i);
  z=-grid->dv[i]+phys->zBold[i];
  //printf("dv[i] =%f, z=%f, z should be %f \n", grid->dv[i], z, -grid->dv[i]+phys->zBold[i]);
  zc_old[i][grid->Nk[i]-1]=z+grid->dzzold[i][grid->Nk[i]-1]/2;
  //printf("zc_old[i][grid->Nk[i]-1]=%f \n", zc_old[i][grid->Nk[i]-1]);
  for(k=grid->Nk[i]-2;k>=grid->ctop[i];k--){
    zc_old[i][k]=zc_old[i][k+1]+grid->dzzold[i][k+1]/2+grid->dzzold[i][k]/2;
  } 

  z=-grid->dv[i]+phys->zBold2[i];
  zc_old2[i][grid->Nk[i]-1]=z+dzz_old2[i][grid->Nk[i]-1]/2;
  k=grid->Nk[i]-1;
  // printf("dzzold[i][grid->Nk[i]-1] = %f, dzz_old2[i][grid->Nk[i]-1] = %f, dzz[i][grid->Nk[i]-1] = %f \n", grid->dzzold[i][grid->Nk[i]-1], dzz_old2[i][grid->Nk[i]-1], grid->dzz[i][grid->Nk[i]-1]);
  // printf("zc[i][k] = %f, zc_old[i][k] = %f, zc_old2[i][k]=%f for i=%i, k=%i \n", vert->zc[i][k], zc_old[i][k], zc_old2[i][k], i, grid->Nk[i]-1);
  // printf("dv[i] =%f \n", grid->dv[i]);

  for(k=grid->Nk[i]-2;k>=grid->ctop[i];k--){
    zc_old2[i][k]=zc_old2[i][k+1]+dzz_old2[i][k+1]/2+dzz_old2[i][k]/2;
    //printf("zc[i][k] = %f, zc_old[i][k] = %f, zc_old2[i][k]=%f for i=%i, k=%i \n", vert->zc[i][k], zc_old[i][k], zc_old2[i][k], i, k);
  }
  
  // for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--){
  //   printf("zc[i][k] = %f for i = %d, k = %d \n", vert->zc[i][k], i, k);
  // }

  //test last timestep...
  // for(k=grid->Nk[i]-1;k>=0;k--){
  //   rho_t2[i][k]=rho_t1[i][k];
  //   rho_t1[i][k]=phys->rho[i][k];

  //   dzz_old2[i][k]=grid->dzzold[i][k];
  //   zc_old2[i][k]=zc_old[i][k];
  //   zc_old[i][k]=vert->zc[i][k];
  // }

  
}

ComputeUCPerot(phys->u_old2,phys->uold2,phys->vold2,phys->h,prop->kinterp,prop->subgrid,grid); //get these for the later part
ComputeUCPerot(phys->u_old,phys->uold,phys->vold,phys->h,prop->kinterp,prop->subgrid,grid); //get these for the later part


// //send and receive rhos 
ISendRecvCellData3D(rho_t1,grid,myproc,comm);
ISendRecvCellData3D(rho_t2,grid,myproc,comm); 
ISendRecvCellData3D(dzz_old2,grid,myproc,comm); 
ISendRecvCellData3D(zc_old,grid,myproc,comm); 
ISendRecvCellData3D(zc_old2,grid,myproc,comm); 
ISendRecvCellData3D(phys->uold2,grid,myproc,comm); 
ISendRecvCellData3D(phys->vold2,grid,myproc,comm); 
ISendRecvCellData3D(phys->uold,grid,myproc,comm); 
ISendRecvCellData3D(phys->vold,grid,myproc,comm); 




//start w montgomery maybe?
if(prop->vertcoord==2) {
    // Montgomery potential is in stmp3
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];
      
      phys->stmp3[i][grid->ctop[i]]=0;
      REAL *zk=phys->a;
      zk[grid->Nk[i]-1]=-grid->dv[i];
      for(k=grid->Nk[i]-1;k>grid->ctop[i];k--)
	      zk[k-1]=zk[k]+grid->dzzold[i][k];
      for(k=grid->ctop[i];k<grid->Nk[i]-1;k++) 
	      phys->stmp3[i][k+1]=phys->stmp3[i][k]+(rho_t1[i][k+1]-rho_t1[i][k])*prop->grav*zk[k];
    

      phys->stmp2[i][grid->ctop[i]]=0;
      zk[grid->Nk[i]-1]=-grid->dv[i];
      for(k=grid->Nk[i]-1;k>grid->ctop[i];k--)
	      zk[k-1]=zk[k]+dzz_old2[i][k];
      for(k=grid->ctop[i];k<grid->Nk[i]-1;k++) 
	      phys->stmp2[i][k+1]=phys->stmp2[i][k]+(rho_t2[i][k+1]-rho_t2[i][k])*prop->grav*zk[k];
    
      }

    ISendRecvCellData3D(phys->stmp3,grid,myproc,comm); //old
    ISendRecvCellData3D(phys->stmp2,grid,myproc,comm); //old2


  }

    

//calculate source terms, assume dzz does not change in time for now
for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
j = grid->edgep[jptr];

nc1 = grid->grad[2*j];
nc2 = grid->grad[2*j+1];

//changed dzzolds to dzz, uold to u, u_old2 to uold

if(prop->vertcoord==2) {
  if(grid->etop[j]<grid->Nke[j]-1) {
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
	    phys->Cn_U[j][k]-=prop->dt*(phys->stmp3[nc1][k]-phys->stmp3[nc2][k])/grid->dg[j];
      phys->Cn_U2[j][k]-=prop->dt*(phys->stmp2[nc1][k]-phys->stmp2[nc2][k])/grid->dg[j];
	  }
  }
} else {
if(grid->etop[j]<grid->Nke[j]-1) 
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      // from the top of the water surface to the bottom edge layer for the given step
      for(k0=grid->etop[j];k0<k;k0++) {
        // for the Cn_U at the particular layer, integrate density gradient over the depth
        // using standard integration to the half cell (extra factor of 1/2 in last term)
        // this is probably only 1st order accurate for the integration
        phys->Cn_U[j][k]-=prop->grav*prop->dt*
          (rho_t1[nc1][k0]*grid->dzzold[nc1][k0]-rho_t1[nc2][k0]*
          grid->dzzold[nc2][k0])/grid->dg[j];

        phys->Cn_U2[j][k]-=prop->grav*prop->dt*
          (rho_t2[nc1][k0]*dzz_old2[nc1][k0]-rho_t2[nc2][k0]*
          dzz_old2[nc2][k0])/grid->dg[j];
      }
      phys->Cn_U[j][k]-=0.5*prop->grav*prop->dt*
        (rho_t1[nc1][k]*grid->dzzold[nc1][k]-rho_t1[nc2][k]*grid->dzzold[nc2][k])/grid->dg[j];

      phys->Cn_U2[j][k]-=0.5*prop->grav*prop->dt*
        (rho_t2[nc1][k]*dzz_old2[nc1][k]-rho_t2[nc2][k]*dzz_old2[nc2][k])/grid->dg[j];
    }
  for(k=grid->etop[j];k<grid->Nke[j];k++) {
    phys->Cn_U[j][k]-=prop->dt*prop->grav*InterpToFace(j,k,rho_t1,phys->u_old,grid)*
      (zc_old[nc1][k]-zc_old[nc2][k])/grid->dg[j];

    phys->Cn_U2[j][k]-=prop->dt*prop->grav*InterpToFace(j,k,rho_t2,phys->u_old2,grid)*
      (zc_old2[nc1][k]-zc_old2[nc2][k])/grid->dg[j];
  }
}  
}

//that was just the beginning part.. now need advective and diffusive
//first for Cn_U

  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

      // now compute for no slip regions for type 4 boundary conditions
  for (jptr = grid->edgedist[4]; jptr < grid->edgedist[5]; jptr++)
  {
    // get index for edge pointers
    j = grid->edgep[jptr];
    //    ib=grid->grad[2*j];
    boundary_index = jptr-grid->edgedist[2];

    // get neighbor indices
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // check to see which of the neighboring cells is the ghost cell
    if (nc1 == -1)  // indicating boundary
      nc = nc2;
    else
      nc = nc1;

    // loop over the entire depth using a ghost cell for the ghost-neighbor
    // on type 4 boundary condition
    for (k=grid->ctop[nc]; k<grid->Nk[nc]; k++){
      
        phys->stmp[nc][k]  += -2.0*prop->nu_H*(
          phys->boundary_u[boundary_index][k] - phys->uold[nc][k])/grid->dg[j]*
          grid->df[j]/grid->Ac[nc];
        phys->stmp2[nc][k] += -2.0*prop->nu_H*(
          phys->boundary_v[boundary_index][k] - phys->vold[nc][k])/grid->dg[j]*
          grid->df[j]/grid->Ac[nc];

    }
  }

  // Now add on horizontal diffusion to stmp and stmp2
  // for the computational cells
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2]; 

    for(k=kmin;k<grid->Nke[j];k++) {
      // Eqn 58 and Eqn 59
      // seems like nu_lax should be distance weighted as in Eqn 60
      a[k]=(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->uold[nc2][k]-phys->uold[nc1][k])*grid->df[j]/grid->dg[j];
      b[k]=(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->vold[nc2][k]-phys->vold[nc1][k])*grid->df[j]/grid->dg[j];

        phys->stmp[nc1][k]-=a[k]/grid->Ac[nc1];
        phys->stmp[nc2][k]+=a[k]/grid->Ac[nc2];
        phys->stmp2[nc1][k]-=b[k]/grid->Ac[nc1];
        phys->stmp2[nc2][k]+=b[k]/grid->Ac[nc2];

    }

    // compute wall-drag for diffusion of momentum BC (only on side walls)
    // Eqn 64 and Eqn 65 and Eqn 58
    for(k=Max(grid->Nke[j],grid->ctop[nc1]);k<grid->Nk[nc1];k++) {

        phys->stmp[nc1][k]+=
          prop->CdW*fabs(phys->uold[nc1][k])*phys->uold[nc1][k]*grid->df[j]/grid->Ac[nc1];
        phys->stmp2[nc1][k]+=
          prop->CdW*fabs(phys->vold[nc1][k])*phys->vold[nc1][k]*grid->df[j]/grid->Ac[nc1];

    }

    for(k=Max(grid->Nke[j],grid->ctop[nc2]);k<grid->Nk[nc2];k++) {

        phys->stmp[nc2][k]+=
          prop->CdW*fabs(phys->uold[nc2][k])*phys->uold[nc2][k]*grid->df[j]/grid->Ac[nc2];
        phys->stmp2[nc2][k]+=
          prop->CdW*fabs(phys->vold[nc2][k])*phys->vold[nc2][k]*grid->df[j]/grid->Ac[nc2]; 

    }

  }

    // Send/recv stmp and stmp2 to account for advective fluxes in ghost cells at
  // interproc boundaries.
  ISendRecvCellData3D(phys->stmp,grid,myproc,comm);
  ISendRecvCellData3D(phys->stmp2,grid,myproc,comm);


  // type 2 boundary condition (specified flux in)
  for(jptr=grid->edgedist[2];jptr<0*grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    i = grid->grad[2*j];
    // zero existing calculations for type 2 edge
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

    sum=0;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1) {
           sum+=grid->Ac[nc];
        for(k=grid->ctop[nc];k<grid->Nk[nc];k++) {
          if(!prop->subgrid)
            Ac=grid->Ac[nc];
          else
            Ac=subgrid->Acceff[nc][k]*grid->dzz[nc][k];
          // get fluxes from other non-boundary cells
          phys->stmp[i][k]+=Ac*phys->stmp[nc][k];
          phys->stmp2[i][k]+=Ac*phys->stmp2[nc][k];
        }
      }
    }
    sum=1/sum;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      if(prop->subgrid)
      {
        sum=0;
        for(nf=0;nf<grid->nfaces[i];nf++)
          if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1)  
            if(k>=grid->ctop[nc])     
              sum+=subgrid->Acceff[nc][k]*grid->dzz[nc][k];
        sum=1/sum;
      } 
      // area-averaged calc
      phys->stmp[i][k]*=sum;
      phys->stmp2[i][k]*=sum;
    }
  }

  // computational cells
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    // Note that dgf==dg only when the cells are orthogonal!
    Return_def(&def1,&def2,nc1,nc2,j,grid);	   	    
    dgf = def1+def2;

    if(grid->ctop[nc1]>grid->ctop[nc2])
      k0=grid->ctop[nc1];
    else
      k0=grid->ctop[nc2];

    // compute momentum advection and diffusion contributions to Cn_U, note
    // the minus sign (why we needed it for the no-slip boundary condition)
    // for each face compute Cn_U performing averaging operation such as in
    // Eqn 41 and Eqn 42 and Eqn 55 etc.
    // the two equations correspond to the two adjacent cells to the edge and 
    // their contributions
    for(k=k0;k<grid->Nk[nc1];k++) 
      phys->Cn_U[j][k]-=def1/dgf
        *prop->dt*(phys->stmp[nc1][k]*grid->n1[j]+phys->stmp2[nc1][k]*grid->n2[j]);
    for(k=k0;k<grid->Nk[nc2];k++) 
      phys->Cn_U[j][k]-=def2/dgf
        *prop->dt*(phys->stmp[nc2][k]*grid->n1[j]+phys->stmp2[nc2][k]*grid->n2[j]);
  }

  // Now add on stmp and stmp2 from the boundaries 
  // for type 3 boundary condition
  for(jptr=grid->edgedist[3];jptr<grid->edgedist[4];jptr++) {
    j = grid->edgep[jptr]; 

    nc1 = grid->grad[2*j];
    k0=grid->ctop[nc1];

    for(nf=0;nf<grid->nfaces[nc1];nf++) {
      if((nc2=grid->neigh[nc1*grid->maxfaces+nf])!=-1) {
        ne=grid->face[nc1*grid->maxfaces+nf];
        for(k=k0;k<grid->Nk[nc1];k++) {
          phys->Cn_U[ne][k]-=
            grid->def[nc1*grid->maxfaces+nf]/grid->dg[ne]*
            prop->dt*(
                phys->stmp[nc2][k]*grid->n1[ne]+phys->stmp2[nc2][k]*grid->n2[ne]);
        }
      }
    }
  }

//that was just the beginning part.. now need advective and diffusive
//Now for Cn_u2!

  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

      // now compute for no slip regions for type 4 boundary conditions
  for (jptr = grid->edgedist[4]; jptr < grid->edgedist[5]; jptr++)
  {
    // get index for edge pointers
    j = grid->edgep[jptr];
    //    ib=grid->grad[2*j];
    boundary_index = jptr-grid->edgedist[2];

    // get neighbor indices
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // check to see which of the neighboring cells is the ghost cell
    if (nc1 == -1)  // indicating boundary
      nc = nc2;
    else
      nc = nc1;

    // loop over the entire depth using a ghost cell for the ghost-neighbor
    // on type 4 boundary condition
    for (k=grid->ctop[nc]; k<grid->Nk[nc]; k++){
      
        phys->stmp[nc][k]  += -2.0*prop->nu_H*(
          phys->boundary_u[boundary_index][k] - phys->uold2[nc][k])/grid->dg[j]*
          grid->df[j]/grid->Ac[nc];
        phys->stmp2[nc][k] += -2.0*prop->nu_H*(
          phys->boundary_v[boundary_index][k] - phys->vold2[nc][k])/grid->dg[j]*
          grid->df[j]/grid->Ac[nc];

    }
  }

  // Now add on horizontal diffusion to stmp and stmp2
  // for the computational cells
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2]; 

    for(k=kmin;k<grid->Nke[j];k++) {
      // Eqn 58 and Eqn 59
      // seems like nu_lax should be distance weighted as in Eqn 60
      a[k]=(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->uold2[nc2][k]-phys->uold2[nc1][k])*grid->df[j]/grid->dg[j];
      b[k]=(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->vold2[nc2][k]-phys->vold2[nc1][k])*grid->df[j]/grid->dg[j];

        phys->stmp[nc1][k]-=a[k]/grid->Ac[nc1];
        phys->stmp[nc2][k]+=a[k]/grid->Ac[nc2];
        phys->stmp2[nc1][k]-=b[k]/grid->Ac[nc1];
        phys->stmp2[nc2][k]+=b[k]/grid->Ac[nc2];

    }

    // compute wall-drag for diffusion of momentum BC (only on side walls)
    // Eqn 64 and Eqn 65 and Eqn 58
    for(k=Max(grid->Nke[j],grid->ctop[nc1]);k<grid->Nk[nc1];k++) {

        phys->stmp[nc1][k]+=
          prop->CdW*fabs(phys->uold2[nc1][k])*phys->uold2[nc1][k]*grid->df[j]/grid->Ac[nc1];
        phys->stmp2[nc1][k]+=
          prop->CdW*fabs(phys->vold2[nc1][k])*phys->vold2[nc1][k]*grid->df[j]/grid->Ac[nc1];

    }

    for(k=Max(grid->Nke[j],grid->ctop[nc2]);k<grid->Nk[nc2];k++) {

        phys->stmp[nc2][k]+=
          prop->CdW*fabs(phys->uold2[nc2][k])*phys->uold2[nc2][k]*grid->df[j]/grid->Ac[nc2];
        phys->stmp2[nc2][k]+=
          prop->CdW*fabs(phys->vold2[nc2][k])*phys->vold2[nc2][k]*grid->df[j]/grid->Ac[nc2]; 

    }

  }

    // Send/recv stmp and stmp2 to account for advective fluxes in ghost cells at
  // interproc boundaries.
  ISendRecvCellData3D(phys->stmp,grid,myproc,comm);
  ISendRecvCellData3D(phys->stmp2,grid,myproc,comm);


  // type 2 boundary condition (specified flux in)
  for(jptr=grid->edgedist[2];jptr<0*grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    i = grid->grad[2*j];
    // zero existing calculations for type 2 edge
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

    sum=0;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1) {
           sum+=grid->Ac[nc];
        for(k=grid->ctop[nc];k<grid->Nk[nc];k++) {
          if(!prop->subgrid)
            Ac=grid->Ac[nc];
          else
            Ac=subgrid->Acceff[nc][k]*grid->dzz[nc][k];
          // get fluxes from other non-boundary cells
          phys->stmp[i][k]+=Ac*phys->stmp[nc][k];
          phys->stmp2[i][k]+=Ac*phys->stmp2[nc][k];
        }
      }
    }
    sum=1/sum;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      if(prop->subgrid)
      {
        sum=0;
        for(nf=0;nf<grid->nfaces[i];nf++)
          if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1)  
            if(k>=grid->ctop[nc])     
              sum+=subgrid->Acceff[nc][k]*grid->dzz[nc][k];
        sum=1/sum;
      } 
      // area-averaged calc
      phys->stmp[i][k]*=sum;
      phys->stmp2[i][k]*=sum;
    }
  }

  // computational cells
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    // Note that dgf==dg only when the cells are orthogonal!
    Return_def(&def1,&def2,nc1,nc2,j,grid);	   	    
    dgf = def1+def2;

    if(grid->ctop[nc1]>grid->ctop[nc2])
      k0=grid->ctop[nc1];
    else
      k0=grid->ctop[nc2];

    // compute momentum advection and diffusion contributions to Cn_U, note
    // the minus sign (why we needed it for the no-slip boundary condition)
    // for each face compute Cn_U performing averaging operation such as in
    // Eqn 41 and Eqn 42 and Eqn 55 etc.
    // the two equations correspond to the two adjacent cells to the edge and 
    // their contributions
    for(k=k0;k<grid->Nk[nc1];k++) 
      phys->Cn_U2[j][k]-=def1/dgf
        *prop->dt*(phys->stmp[nc1][k]*grid->n1[j]+phys->stmp2[nc1][k]*grid->n2[j]);
    for(k=k0;k<grid->Nk[nc2];k++) 
      phys->Cn_U2[j][k]-=def2/dgf
        *prop->dt*(phys->stmp[nc2][k]*grid->n1[j]+phys->stmp2[nc2][k]*grid->n2[j]);
  }

  // Now add on stmp and stmp2 from the boundaries 
  // for type 3 boundary condition
  for(jptr=grid->edgedist[3];jptr<grid->edgedist[4];jptr++) {
    j = grid->edgep[jptr]; 

    nc1 = grid->grad[2*j];
    k0=grid->ctop[nc1];

    for(nf=0;nf<grid->nfaces[nc1];nf++) {
      if((nc2=grid->neigh[nc1*grid->maxfaces+nf])!=-1) {
        ne=grid->face[nc1*grid->maxfaces+nf];
        for(k=k0;k<grid->Nk[nc1];k++) {
          phys->Cn_U2[ne][k]-=
            grid->def[nc1*grid->maxfaces+nf]/grid->dg[ne]*
            prop->dt*(
                phys->stmp[nc2][k]*grid->n1[ne]+phys->stmp2[nc2][k]*grid->n2[ne]);
        }
      }
    }
  }


for(i=0;i<grid->Nc;i++) {
  free(rho_t1[i]);
  free(rho_t2[i]);
  free(s_old2[i]);
  free(dzz_old2[i]);
  free(zc_old[i]);
  free(zc_old2[i]);
}
free(rho_t1);
free(rho_t2);
free(s_old2);
free(dzz_old2);
free(zc_old);
free(zc_old2);



}


/*
 * Function: HorizontalSource
 * Usage: HorizontalSource(grid,phys,prop,myproc,numprocs);
 * --------------------------------------------------------
 * Compute the horizontal source term that is used to obtain the free surface.
 *
 * This function adds the following to the horizontal source term:
 *
 * 1) Old nonhydrostatic pressure gradient with theta method
 * 2) Coriolis terms with AB2
 * 3) Baroclinic term with AB2
 * 4) Horizontal and vertical advection of horizontal momentum with AB2
 * 5) Horizontal laminar+turbulent diffusion of horizontal momentum
 *
 * Cn_U contains the Adams-Bashforth terms at time step n-1.
 * If wetting and drying is employed, no advection is computed in 
 * the upper cell.
 *
 */
static void HorizontalSource(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, MPI_Comm comm) 
{
  int i, ib, iptr, boundary_index, nf, j, jptr, k, nc, nc1, nc2, ne, 
  k0, kmin, kmax;
  REAL *a, *b, *c, fab1, fab2, fab3, sum, def1, def2, dgf, Cz, tempu,Ac,f_sum, ke1,ke2,Vm; //AB3
  // additions to test divergence averaging for w in momentum calc
  REAL wedge[3], lambda[3], wik, tmp_x, tmp_y,tmp;
  int aneigh;

  a = phys->a;
  b = phys->b;
  c = phys->c;

  // get the Adams-Bashforth multi-step integration started
  // fab is 1 for a forward Euler calculation on the first time step,
  // for which Cn_U is 0.  Otherwise, fab=3/2 and Cn_U contains the
  // Adams-Bashforth terms at time step n-1

 // Adams Bashforth coefficients
  if(!prop->readOldVelocity){
    //if(1){
  if(prop->n==1 || prop->wetdry) {
    fab1=1;
    fab2=fab3=0;

    for(j=0;j<grid->Ne;j++)
      for(k=0;k<grid->Nke[j];k++)
        phys->Cn_U[j][k]=phys->Cn_U2[j][k]=0;
  } else if(prop->n==2) {
    fab1=3.0/2.0;
    fab2=-1.0/2.0;
    fab3=0;
  } else {
      fab1=prop->exfac1;
      fab2=prop->exfac2;
      fab3=prop->exfac3;   
    }
  } else{
    fab1=prop->exfac1;
    fab2=prop->exfac2;
    fab3=prop->exfac3;   
  }

  // Set utmp and ut to zero since utmp will store the source term of the
  // horizontal momentum equation
  for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nke[j];k++) {
      phys->utmp[j][k]=0;
      phys->ut[j][k]=0;
    }
  }

  // Update with old AB term
  // correct velocity based on non-hydrostatic pressure
  // over all computational edges
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    //AB3
    // for each edge over depth
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      // Equation 39: U_j,k^n+1 = U_j,k^* - dt*(qc_G2j,k - qc_G1j,k)/Dj
      // note that EE is used for the first time step
      phys->utmp[j][k]=fab2*phys->Cn_U[j][k]+fab3*phys->Cn_U2[j][k]+phys->u[j][k];
      phys->Cn_U2[j][k]=phys->Cn_U[j][k];
      phys->Cn_U[j][k]=0;
    }
  }

  // Only include old pressure gradient for correction method
  if(prop->pressureMethod==1) {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	phys->utmp[j][k]-=prop->dt/grid->dg[j]*(phys->q[nc1][k]-phys->q[nc2][k]);
      }
    }
  }

  // Add on explicit term to boundary edges (type 4 BCs)
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr]; 

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->utmp[j][k]=fab3*phys->Cn_U2[j][k]+fab2*phys->Cn_U[j][k]+phys->u[j][k];
      phys->Cn_U2[j][k]=phys->Cn_U[j][k];
      phys->Cn_U[j][k]=0;
    }
  }

  // note that the above lines appear to allow use to "flush" Cn_U so
  // that it's values are utilizes and then it is free for additional
  // computations

  // Add on a momentum source to momentum equation
  // currently covers the sponge layer and can be used for Coriolis for the 
  // 2D problem
  MomentumSource(phys->utmp,grid,phys,prop);
  
  // 3D Coriolis terms
  // note that this uses linear interpolation to the faces from the cell centers
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) 
  {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    for(k=grid->etop[j];k<grid->Nke[j];k++)
    {
      // if not z-level add relative voricity due to momentum advection
      if(prop->vertcoord!=1 && prop->nonlinear && prop->wetdry)  
        f_sum=prop->Coriolis_f+vert->f_re[j][k];
      else
        f_sum=prop->Coriolis_f;
      phys->Cn_U[j][k]+=prop->dt*f_sum*(
          InterpToFace(j,k,phys->vc,phys->u,grid)*grid->n1[j]-
          InterpToFace(j,k,phys->uc,phys->u,grid)*grid->n2[j]);
    }
  }

  if(prop->vertcoord==2) {
    // Montgomery potential is in stmp3
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];
      
      phys->stmp3[i][grid->ctop[i]]=0;
      REAL *zk=phys->a;
      zk[grid->Nk[i]-1]=-grid->dv[i];
      for(k=grid->Nk[i]-1;k>grid->ctop[i];k--)
	zk[k-1]=zk[k]+grid->dzz[i][k];
      for(k=grid->ctop[i];k<grid->Nk[i]-1;k++) 
	phys->stmp3[i][k+1]=phys->stmp3[i][k]+(phys->rho[i][k+1]-phys->rho[i][k])*prop->grav*zk[k];
    }
  } else {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];
      
      phys->stmp3[i][grid->ctop[i]]=0.5*grid->dzz[i][grid->ctop[i]]*phys->rho[i][grid->ctop[i]];
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	phys->stmp3[i][k]=phys->stmp3[i][k-1]+0.5*(grid->dzz[i][k-1]*phys->rho[i][k-1]+
						   grid->dzz[i][k]*phys->rho[i][k]);
      }
    }
  }    
  ISendRecvCellData3D(phys->stmp3,grid,myproc,comm);
    
  // Baroclinic term
  // over computational cells
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // why only has baroclinic pressure when Nke=1
    if(prop->vertcoord==1 || prop->vertcoord==5){
      if(grid->etop[j]<grid->Nke[j]-1) 
        for(k=grid->etop[j];k<grid->Nke[j];k++) {
          // from the top of the water surface to the bottom edge layer for the given step
          for(k0=Max(grid->ctop[nc1],grid->ctop[nc2]);k0<k;k0++) {
            // for the Cn_U at the particular layer, integrate density gradient over the depth
            // using standard integration to the half cell (extra factor of 1/2 in last term)
            // this is probably only 1st order accurate for the integration
            phys->Cn_U[j][k]-=0.5*prop->grav*prop->dt*
              (phys->rho[nc1][k0]-phys->rho[nc2][k0])*
              (grid->dzz[nc1][k0]+grid->dzz[nc2][k0])/grid->dg[j];
          }
          phys->Cn_U[j][k]-=0.25*prop->grav*prop->dt*
            (phys->rho[nc1][k]-phys->rho[nc2][k])*
            (grid->dzz[nc1][k]+grid->dzz[nc2][k])/grid->dg[j];
        }
    } else if(prop->vertcoord==2) {
      if(grid->etop[j]<grid->Nke[j]-1) {
        for(k=grid->etop[j];k<grid->Nke[j];k++) {
	  phys->Cn_U[j][k]-=prop->dt*(phys->stmp3[nc1][k]-phys->stmp3[nc2][k])/grid->dg[j];
	}
      }
    } else {
      /*
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	phys->Cn_U[j][k]-=prop->dt*prop->grav*(0.5*(phys->rho[nc1][k]+phys->rho[nc2][k])*
					       (vert->zc[nc1][k]-vert->zc[nc2][k])/grid->dg[j]+
					       (phys->stmp3[nc1][k]-phys->stmp3[nc2][k])/grid->dg[j]);
      }
      */
      if(grid->etop[j]<grid->Nke[j]-1) 
        for(k=grid->etop[j];k<grid->Nke[j];k++) {
          // from the top of the water surface to the bottom edge layer for the given step
          for(k0=grid->etop[j];k0<k;k0++) {
            // for the Cn_U at the particular layer, integrate density gradient over the depth
            // using standard integration to the half cell (extra factor of 1/2 in last term)
            // this is probably only 1st order accurate for the integration
            phys->Cn_U[j][k]-=prop->grav*prop->dt*
              (phys->rho[nc1][k0]*grid->dzz[nc1][k0]-phys->rho[nc2][k0]*
	       grid->dzz[nc2][k0])/grid->dg[j];
          }
          phys->Cn_U[j][k]-=0.5*prop->grav*prop->dt*
            (phys->rho[nc1][k]*grid->dzz[nc1][k]-phys->rho[nc2][k]*grid->dzz[nc2][k])/grid->dg[j];
        }
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	phys->Cn_U[j][k]-=prop->dt*prop->grav*InterpToFace(j,k,phys->rho,phys->u,grid)*
	  (vert->zc[nc1][k]-vert->zc[nc2][k])/grid->dg[j];
      }
    }      
  }

  // Set stmp and stmp2 to zero since these are used as temporary variables for advection and
  // diffusion.
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

  // new scheme for mometum advection if not z-level
  // keep the advection term conservative while add the additionterm u/J*dJ/dt
  // comment out since it is solved by implicit method, see function Upredictor

  if(prop->nonlinear && prop->vertcoord!=1)
  {
    if(!prop->wetdry)
    {
      if(vert->dJdtmeth==1)
      {
        for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) 
        {
           j = grid->edgep[jptr];
           nc1 = grid->grad[2*j];
           nc2 = grid->grad[2*j+1];
	   Return_def(&def1,&def2,nc1,nc2,j,grid);	   
           dgf = def1+def2;
           for(k=grid->etop[j];k<grid->Nke[j];k++) 
             phys->Cn_U[j][k]-=phys->u[j][k]*
             (def2/dgf*(1-grid->dzzold[nc1][k]/grid->dzz[nc1][k])+def1/dgf*(1-grid->dzzold[nc2][k]/grid->dzz[nc2][k]));
        }
      }
    } else {
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) 
      {
        j = grid->edgep[jptr];
        nc1 = grid->grad[2*j];
        nc2 = grid->grad[2*j+1];
	Return_def(&def1,&def2,nc1,nc2,j,grid);	   	
        for(k=grid->etop[j];k<grid->Nke[j];k++) 
        { 
	  /* 
          Vm=def1*grid->dzz[nc1][k]+def2*grid->dzz[nc2][k];
          ke1=(phys->uc[nc1][k]*phys->uc[nc1][k]+phys->vc[nc1][k]*phys->vc[nc1][k])*grid->dzz[nc1][k]/2;
          ke2=(phys->uc[nc2][k]*phys->uc[nc2][k]+phys->vc[nc2][k]*phys->vc[nc2][k])*grid->dzz[nc2][k]/2;
          phys->Cn_U[j][k]-=prop->dt*(ke1-ke2)/Vm;
	  */
          ke1=0.5*(pow(phys->uc[nc1][k],2)+pow(phys->vc[nc1][k],2));
          ke2=0.5*(pow(phys->uc[nc2][k],2)+pow(phys->vc[nc2][k],2));
	  phys->Cn_U[j][k]-=prop->dt*(ke1-ke2)/grid->dg[j];
        }  
      }
    }
  }
  // Compute Eulerian advection of momentum (nonlinear!=0)
  if(prop->nonlinear && (prop->vertcoord==1 || (prop->vertcoord!=1 && !prop->wetdry))) 
  {

    // Interpolate uc to faces and place into ut
    GetMomentumFaceValues(phys->ut,phys->uc,phys->boundary_u,phys->u,grid,phys,prop,comm,myproc,prop->nonlinear, prop->TVDmomentum);

    // Conservative method assumes ut is a flux
    if(prop->conserveMomentum)
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) 
      {
        j=grid->edgep[jptr];

        for(k=grid->etop[j];k<grid->Nke[j];k++){
          phys->ut[j][k]*=grid->dzf[j][k];
        }
      }

    // Now compute the cell-centered source terms and put them into stmp
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];
      // Store dzz in a since for conservative scheme need to divide by depth (since ut is a flux)
      if(prop->conserveMomentum) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          a[k]=grid->dzz[i][k];
      } else {
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          a[k]=1.0;
      }

      // for each face
      for(nf=0;nf<grid->nfaces[i];nf++) {

        // get the edge pointer
        ne = grid->face[i*grid->maxfaces+nf];

        // for all but the top cell layer
        for(k=grid->ctop[i];k<grid->Nk[i];k++) 
        {
          // this is basically Eqn 50, u-component
          if(!prop->subgrid || prop->wetdry)
            Ac=grid->Ac[i];
          else
            Ac=subgrid->Acceff[i][k];
          phys->stmp[i][k]+=
            phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/(a[k]*Ac);         
        }  

        // Top cell is filled with momentum from neighboring cells
        if(prop->conserveMomentum)
          for(k=grid->etop[ne];k<grid->ctop[i];k++) 
          {
            if(!prop->subgrid || prop->wetdry)
              Ac=grid->Ac[i];
            else
              Ac=subgrid->Acceff[i][k];	  
            phys->stmp[i][grid->ctop[i]]+=
            phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/(a[grid->ctop[i]]*Ac);
          }
      }
    }

    // Interpolate vc to faces and place into ut
    GetMomentumFaceValues(phys->ut,phys->vc,phys->boundary_v,phys->u,grid,phys,prop,comm,myproc,prop->nonlinear,prop->TVDmomentum);

    // Conservative method assumes ut is a flux
    if(prop->conserveMomentum)
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
	     j=grid->edgep[jptr];
	     for(k=grid->etop[j];k<grid->Nke[j];k++)
	       phys->ut[j][k]*=grid->dzf[j][k];
      }

    // Now compute the cell-centered source terms and put them into stmp.
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];

      for(k=0;k<grid->Nk[i];k++) 
        phys->stmp2[i][k]=0;

      // Store dzz in a since for conservative scheme need to divide by depth (since ut is a flux)
      if(prop->conserveMomentum) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          a[k]=grid->dzz[i][k];
      } else {
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          a[k]=1.0;
      }

      for(nf=0;nf<grid->nfaces[i];nf++) {

        ne = grid->face[i*grid->maxfaces+nf];

        for(k=grid->ctop[i];k<grid->Nk[i];k++)
        {
          if(!prop->subgrid || prop->wetdry)
            Ac=grid->Ac[i];
          else
            Ac=subgrid->Acceff[i][k];
          // Eqn 50, v-component
          phys->stmp2[i][k]+=
            phys->ut[ne][k]*phys->u[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/(a[k]*Ac);
        }
        // Top cell is filled with momentum from neighboring cells
        if(prop->conserveMomentum)
          for(k=grid->etop[ne];k<grid->ctop[i];k++) 
          {
            if(!prop->subgrid || prop->wetdry)
              Ac=grid->Ac[i];
            else
              Ac=subgrid->Acceff[i][k];
            phys->stmp2[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->u[ne][k]*
              grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/(a[k]*Ac);
          }
      }
    }

    // stmp2 holds the v-component of advection of momentum and
    // note that this could also be sped up by pulling out common
    // factors to reduce redundant calculations

    /* Vertical advection of momentum calc */
    // Only if thetaM<0, otherwise use implicit scheme in UPredictor()

    if(prop->thetaM<0 && prop->vertcoord==1) {
      // Now do vertical advection of momentum
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i=grid->cellp[iptr];
        switch(prop->nonlinear) {
      	  case 1:
            for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
              a[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->uc[i][k]+
                (phys->w[i][k]-fabs(phys->w[i][k]))*phys->uc[i][k-1]);
              b[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->vc[i][k]+
                (phys->w[i][k]-fabs(phys->w[i][k]))*phys->vc[i][k-1]);
            }
            break;
          case 2: case 5:
            for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
              a[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k]+
                grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k-1]));
              b[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k]+
                grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k-1]));
            }
            break;
        	case 4:
            for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
              Cz = 2.0*phys->w[i][k]*prop->dt/(grid->dzz[i][k]+grid->dzz[i][k-1]);
              a[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k]+
                grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k-1])-
                0.5*Cz*(phys->uc[i][k-1]-phys->uc[i][k]));
              b[k] = phys->w[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k]+
                grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k-1])-
                0.5*Cz*(phys->vc[i][k-1]-phys->vc[i][k]));
            }
            break;
        	default:
            for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
              a[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->uc[i][k]+
                (phys->w[i][k]-fabs(phys->w[i][k]))*phys->uc[i][k-1]);
              b[k] = 0.5*((phys->w[i][k]+fabs(phys->w[i][k]))*phys->vc[i][k]+
                (phys->w[i][k]-fabs(phys->w[i][k]))*phys->vc[i][k-1]);
            }
            break;
        }
        // Always do first-order upwind in bottom cell if partial stepping is on
        if(prop->stairstep==0) {
          k = grid->Nk[i]-1;
          a[k] = 0.5*(
            (phys->w[i][k]+fabs(phys->w[i][k]))*phys->uc[i][k]+
            (phys->w[i][k]-fabs(phys->w[i][k]))*phys->uc[i][k-1]);
          b[k] = 0.5*(
            (phys->w[i][k]+fabs(phys->w[i][k]))*phys->vc[i][k]+
            (phys->w[i][k]-fabs(phys->w[i][k]))*phys->vc[i][k-1]);
        }
        a[grid->ctop[i]]=phys->w[i][grid->ctop[i]]*phys->uc[i][grid->ctop[i]];
        b[grid->ctop[i]]=phys->w[i][grid->ctop[i]]*phys->vc[i][grid->ctop[i]];
        a[grid->Nk[i]]=0;
        b[grid->Nk[i]]=0;

        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          if(prop->subgrid)
          {
            // if wetting and drying happens use Ac insteady of acceff
            if(!prop->wetdry){
              phys->stmp[i][k]+=(a[k]*subgrid->Acveff[i][k]-a[k+1]*subgrid->Acveff[i][k+1])/grid->dzz[i][k]/subgrid->Acceff[i][k];
              phys->stmp2[i][k]+=(b[k]*subgrid->Acveff[i][k]-b[k+1]*subgrid->Acveff[i][k+1])/grid->dzz[i][k]/subgrid->Acceff[i][k];                
            } else {
              phys->stmp[i][k]+=(a[k]*subgrid->Acveff[i][k]-a[k+1]*subgrid->Acveff[i][k+1])/grid->dzz[i][k]/grid->Ac[i];
              phys->stmp2[i][k]+=(b[k]*subgrid->Acveff[i][k]-b[k+1]*subgrid->Acveff[i][k+1])/grid->dzz[i][k]/grid->Ac[i];                  
            }
          } else {
            phys->stmp[i][k]+=(a[k]-a[k+1])/grid->dzz[i][k];
            phys->stmp2[i][k]+=(b[k]-b[k+1])/grid->dzz[i][k];
          }
        }
      }
    } 
   
    // add explicit form for vertical momentum advection
    if(prop->thetaM<0 && prop->vertcoord!=1) {
      // Now do vertical advection of momentum
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i=grid->cellp[iptr];
        switch(prop->nonlinear) {
          case 1:
            for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
              a[k] = 0.5*((vert->omega_old[i][k]+fabs(vert->omega_old[i][k]))*phys->uc[i][k]+
                (vert->omega_old[i][k]-fabs(vert->omega_old[i][k]))*phys->uc[i][k-1]);
              b[k] = 0.5*((vert->omega_old[i][k]+fabs(vert->omega_old[i][k]))*phys->vc[i][k]+
                (vert->omega_old[i][k]-fabs(vert->omega_old[i][k]))*phys->vc[i][k-1]);
            }
            break;
          case 2: case 5:
            for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
              a[k] = vert->omega_old[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k]+
                grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k-1]));
              b[k] = vert->omega_old[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k]+
                grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k-1]));
            }
            break;

          case 4:
            for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
              Cz = 2.0*vert->omega_old[i][k]*prop->dt/(grid->dzz[i][k]+grid->dzz[i][k-1]);
              a[k] = vert->omega_old[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k]+
                grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->uc[i][k-1])-
                0.5*Cz*(phys->uc[i][k-1]-phys->uc[i][k]));
              b[k] = vert->omega_old[i][k]*((grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k]+
                grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])*phys->vc[i][k-1])-
                0.5*Cz*(phys->vc[i][k-1]-phys->vc[i][k]));
            }
            break;
          default:
            for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
              a[k] = 0.5*((vert->omega_old[i][k]+fabs(vert->omega_old[i][k]))*phys->uc[i][k]+
                (vert->omega_old[i][k]-fabs(vert->omega_old[i][k]))*phys->uc[i][k-1]);
              b[k] = 0.5*((vert->omega_old[i][k]+fabs(vert->omega_old[i][k]))*phys->vc[i][k]+
                (vert->omega_old[i][k]-fabs(vert->omega_old[i][k]))*phys->vc[i][k-1]);
            }
            break;
        }
        
        // Always do first-order upwind in bottom cell if partial stepping is on
        if(prop->stairstep==0) {
          k = grid->Nk[i]-1;
          a[k] = 0.5*(
                      (vert->omega_old[i][k]+fabs(vert->omega_old[i][k]))*phys->uc[i][k]+
                      (vert->omega_old[i][k]-fabs(vert->omega_old[i][k]))*phys->uc[i][k-1]);
          b[k] = 0.5*(
                      (vert->omega_old[i][k]+fabs(vert->omega_old[i][k]))*phys->vc[i][k]+
                      (vert->omega_old[i][k]-fabs(vert->omega_old[i][k]))*phys->vc[i][k-1]);
        }
        
        a[grid->ctop[i]]=vert->omega_old[i][grid->ctop[i]]*phys->uc[i][grid->ctop[i]];
        b[grid->ctop[i]]=vert->omega_old[i][grid->ctop[i]]*phys->vc[i][grid->ctop[i]];
        a[grid->Nk[i]]=0;
        b[grid->Nk[i]]=0;
        
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          if(prop->subgrid)
          {
            // if wetting and drying happens use Ac insteady of acceff
            if(!prop->wetdry){
              phys->stmp[i][k]+=(a[k]*subgrid->Acveff[i][k]-a[k+1]*subgrid->Acveff[i][k+1])/grid->dzz[i][k]/subgrid->Acceff[i][k];
              phys->stmp2[i][k]+=(b[k]*subgrid->Acveff[i][k]-b[k+1]*subgrid->Acveff[i][k+1])/grid->dzz[i][k]/subgrid->Acceff[i][k];                
            } else {
              phys->stmp[i][k]+=(a[k]*subgrid->Acveff[i][k]-a[k+1]*subgrid->Acveff[i][k+1])/grid->dzz[i][k]/grid->Ac[i];
              phys->stmp2[i][k]+=(b[k]*subgrid->Acveff[i][k]-b[k+1]*subgrid->Acveff[i][k+1])/grid->dzz[i][k]/grid->Ac[i];                  
            }

          } else {
            phys->stmp[i][k]+=(a[k]-a[k+1])/grid->dzz[i][k];
            phys->stmp2[i][k]+=(b[k]-b[k+1])/grid->dzz[i][k];
          }
        }
      }
    } // end of nonlinear computation
  }

  // stmp and stmp2 just store the summed C_H and C_V values for horizontal
  // advection (prior to utilization via Eqn 41 or Eqn 47)

  /* Horizontal diffusion calculations */

  // now compute for no slip regions for type 4 boundary conditions
  for (jptr = grid->edgedist[4]; jptr < grid->edgedist[5]; jptr++)
  {
    // get index for edge pointers
    j = grid->edgep[jptr];
    //    ib=grid->grad[2*j];
    boundary_index = jptr-grid->edgedist[2];

    // get neighbor indices
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // check to see which of the neighboring cells is the ghost cell
    if (nc1 == -1)  // indicating boundary
      nc = nc2;
    else
      nc = nc1;

    // loop over the entire depth using a ghost cell for the ghost-neighbor
    // on type 4 boundary condition
    for (k=grid->ctop[nc]; k<grid->Nk[nc]; k++){
      if(!prop->subgrid){
        phys->stmp[nc][k]  += -2.0*prop->nu_H*(
          phys->boundary_u[boundary_index][k] - phys->uc[nc][k])/grid->dg[j]*
          grid->df[j]/grid->Ac[nc];
        phys->stmp2[nc][k] += -2.0*prop->nu_H*(
          phys->boundary_v[boundary_index][k] - phys->vc[nc][k])/grid->dg[j]*
          grid->df[j]/grid->Ac[nc];
      } else {
        phys->stmp[nc][k]  += -2.0*prop->nu_H*(
          phys->boundary_u[boundary_index][k] - phys->uc[nc][k])/grid->dg[j]*grid->dzf[j][k]*
          grid->df[j]/subgrid->Acceff[nc][k]/grid->dzz[nc][k];
        phys->stmp2[nc][k] += -2.0*prop->nu_H*(
          phys->boundary_v[boundary_index][k] - phys->vc[nc][k])/grid->dg[j]*grid->dzf[j][k]*
          grid->df[j]/subgrid->Acceff[nc][k]/grid->dzz[nc][k];        
      }
    }
  }

  // Now add on horizontal diffusion to stmp and stmp2
  // for the computational cells
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2]; 

    for(k=kmin;k<grid->Nke[j];k++) {
      // Eqn 58 and Eqn 59
      // seems like nu_lax should be distance weighted as in Eqn 60
      a[k]=(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->uc[nc2][k]-phys->uc[nc1][k])*grid->df[j]/grid->dg[j];
      b[k]=(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->vc[nc2][k]-phys->vc[nc1][k])*grid->df[j]/grid->dg[j];

      if(!prop->subgrid)
      {
        phys->stmp[nc1][k]-=a[k]/grid->Ac[nc1];
        phys->stmp[nc2][k]+=a[k]/grid->Ac[nc2];
        phys->stmp2[nc1][k]-=b[k]/grid->Ac[nc1];
        phys->stmp2[nc2][k]+=b[k]/grid->Ac[nc2];
      }else{
        phys->stmp[nc1][k]-=a[k]/subgrid->Acceff[nc1][k]/grid->dzz[nc1][k]*grid->dzf[j][k];
        phys->stmp[nc2][k]+=a[k]/subgrid->Acceff[nc2][k]/grid->dzz[nc2][k]*grid->dzf[j][k];
        phys->stmp2[nc1][k]-=b[k]/subgrid->Acceff[nc1][k]/grid->dzz[nc1][k]*grid->dzf[j][k];
        phys->stmp2[nc2][k]+=b[k]/subgrid->Acceff[nc2][k]/grid->dzz[nc2][k]*grid->dzf[j][k];
      }
    }

    // compute wall-drag for diffusion of momentum BC (only on side walls)
    // Eqn 64 and Eqn 65 and Eqn 58
    for(k=Max(grid->Nke[j],grid->ctop[nc1]);k<grid->Nk[nc1];k++) {
      if(!prop->subgrid)
      {
        phys->stmp[nc1][k]+=
          prop->CdW*fabs(phys->uc[nc1][k])*phys->uc[nc1][k]*grid->df[j]/grid->Ac[nc1];
        phys->stmp2[nc1][k]+=
          prop->CdW*fabs(phys->vc[nc1][k])*phys->vc[nc1][k]*grid->df[j]/grid->Ac[nc1];
      } else {
        phys->stmp[nc1][k]+=
          prop->CdW*fabs(phys->uc[nc1][k])*phys->uc[nc1][k]*grid->df[j]*grid->dzf[j][k]/subgrid->Acceff[nc1][k]/grid->dzz[nc1][k];
        phys->stmp2[nc1][k]+=
          prop->CdW*fabs(phys->vc[nc1][k])*phys->vc[nc1][k]*grid->df[j]*grid->dzf[j][k]/subgrid->Acceff[nc1][k]/grid->dzz[nc1][k];
      }
    }

    for(k=Max(grid->Nke[j],grid->ctop[nc2]);k<grid->Nk[nc2];k++) {
      if(!prop->subgrid)
      {
        phys->stmp[nc2][k]+=
          prop->CdW*fabs(phys->uc[nc2][k])*phys->uc[nc2][k]*grid->df[j]/grid->Ac[nc2];
        phys->stmp2[nc2][k]+=
          prop->CdW*fabs(phys->vc[nc2][k])*phys->vc[nc2][k]*grid->df[j]/grid->Ac[nc2]; 
      } else {
        phys->stmp[nc2][k]+=
          prop->CdW*fabs(phys->uc[nc2][k])*phys->uc[nc2][k]*grid->df[j]*grid->dzf[j][k]/subgrid->Acceff[nc2][k]/grid->dzz[nc2][k];
        phys->stmp2[nc2][k]+=
          prop->CdW*fabs(phys->vc[nc2][k])*phys->vc[nc2][k]*grid->df[j]*grid->dzf[j][k]/subgrid->Acceff[nc2][k]/grid->dzz[nc2][k];         
      }
    }

  }

  // Check to make sure integrated fluxes are 0 for conservation
  // This will not be conservative if CdW or nu_H are nonzero!
  if(WARNING && prop->CdW==0 && prop->nu_H==0) {
    sum=0;
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        if(!prop->subgrid)
          sum+=grid->Ac[i]*phys->stmp[i][k]*grid->dzz[i][k];
        else
          sum+=subgrid->Acceff[i][k]*phys->stmp[i][k]*grid->dzz[i][k];      
    }
    if(fabs(sum)>CONSERVED) printf("Warning, not U-momentum conservative!\n");

    sum=0;
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        if(!prop->subgrid)
          sum+=grid->Ac[i]*phys->stmp2[i][k]*grid->dzz[i][k];
        else
          sum+=subgrid->Acceff[i][k]*phys->stmp2[i][k]*grid->dzz[i][k];
    }
    if(fabs(sum)>CONSERVED) printf("Warning, not V-momentum conservative!\n");
  }

  // Send/recv stmp and stmp2 to account for advective fluxes in ghost cells at
  // interproc boundaries.
  ISendRecvCellData3D(phys->stmp,grid,myproc,comm);
  ISendRecvCellData3D(phys->stmp2,grid,myproc,comm);

  // type 2 boundary condition (specified flux in)
  for(jptr=grid->edgedist[2];jptr<0*grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    i = grid->grad[2*j];
    // zero existing calculations for type 2 edge
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

    sum=0;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1) {
           sum+=grid->Ac[nc];
        for(k=grid->ctop[nc];k<grid->Nk[nc];k++) {
          if(!prop->subgrid)
            Ac=grid->Ac[nc];
          else
            Ac=subgrid->Acceff[nc][k]*grid->dzz[nc][k];
          // get fluxes from other non-boundary cells
          phys->stmp[i][k]+=Ac*phys->stmp[nc][k];
          phys->stmp2[i][k]+=Ac*phys->stmp2[nc][k];
        }
      }
    }
    sum=1/sum;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      if(prop->subgrid)
      {
        sum=0;
        for(nf=0;nf<grid->nfaces[i];nf++)
          if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1)  
            if(k>=grid->ctop[nc])     
              sum+=subgrid->Acceff[nc][k]*grid->dzz[nc][k];
        sum=1/sum;
      } 
      // area-averaged calc
      phys->stmp[i][k]*=sum;
      phys->stmp2[i][k]*=sum;
    }
  }

  // computational cells
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    // Note that dgf==dg only when the cells are orthogonal!
    Return_def(&def1,&def2,nc1,nc2,j,grid);	   	    
    dgf = def1+def2;

    if(grid->ctop[nc1]>grid->ctop[nc2])
      k0=grid->ctop[nc1];
    else
      k0=grid->ctop[nc2];

    // compute momentum advection and diffusion contributions to Cn_U, note
    // the minus sign (why we needed it for the no-slip boundary condition)
    // for each face compute Cn_U performing averaging operation such as in
    // Eqn 41 and Eqn 42 and Eqn 55 etc.
    // the two equations correspond to the two adjacent cells to the edge and 
    // their contributions
    for(k=k0;k<grid->Nk[nc1];k++) 
      phys->Cn_U[j][k]-=def1/dgf
        *prop->dt*(phys->stmp[nc1][k]*grid->n1[j]+phys->stmp2[nc1][k]*grid->n2[j]);
    for(k=k0;k<grid->Nk[nc2];k++) 
      phys->Cn_U[j][k]-=def2/dgf
        *prop->dt*(phys->stmp[nc2][k]*grid->n1[j]+phys->stmp2[nc2][k]*grid->n2[j]);
  }

  // Now add on stmp and stmp2 from the boundaries 
  // for type 3 boundary condition
  for(jptr=grid->edgedist[3];jptr<grid->edgedist[4];jptr++) {
    j = grid->edgep[jptr]; 

    nc1 = grid->grad[2*j];
    k0=grid->ctop[nc1];

    for(nf=0;nf<grid->nfaces[nc1];nf++) {
      if((nc2=grid->neigh[nc1*grid->maxfaces+nf])!=-1) {
        ne=grid->face[nc1*grid->maxfaces+nf];
        for(k=k0;k<grid->Nk[nc1];k++) {
          phys->Cn_U[ne][k]-=
            grid->def[nc1*grid->maxfaces+nf]/grid->dg[ne]*
            prop->dt*(
                phys->stmp[nc2][k]*grid->n1[ne]+phys->stmp2[nc2][k]*grid->n2[ne]);
        }
      }
    }
  }

  // note that we now basically have the term dt*F_j,k in Equation 33
  // update utmp 
  // this will complete the adams-bashforth time stepping
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr]; 

    for(k=grid->etop[j];k<grid->Nke[j];k++){
      phys->utmp[j][k]+=fab1*phys->Cn_U[j][k];
    }  
  }
  
//  // check to make sure we don't have a blow-up
//  for(j=0;j<grid->Ne;j++) 
//    for(k=grid->etop[j];k<grid->Nke[j];k++) 
//      if(phys->utmp[j][k]!=phys->utmp[j][k]) {
//        printf("0Error in utmp at j=%d k=%d (U***=nan)\n",j,k);
//        exit(1);
//      }
}

/*
 * Function: NewCells
 * Usage: NewCells(grid,phys,prop);
 * --------------------------------
 * Adjust the velocity in the new cells that were previously dry.
 * This function is required for an Eulerian advection scheme because
 * it is difficult to compute the finite volume form of advection in the
 * upper cells when wetting and drying is employed.  In this function
 * the velocity in the new cells is set such that the quantity u*dz is
 * conserved from one time step to the next.  This works well without
 * wetting and drying.  When wetting and drying is employed, it is best
 * to extrapolate from the lower cells to obtain the velocity in the new
 * cells.
 *
 */
static void NewCells(gridT *grid, physT *phys, propT *prop) {
  int j, jptr, k;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    if(grid->etop[j]<grid->Nke[j]-1) {    // Only update new cells if there is more than one layer

      // If a new cell is created above the old cell then set the velocity in the new cell equal to that in
      // the old one.
      //
      // Otherwise if the free surface jumps up more than one cell then set the velocity equal to the cell beneath
      // the old one.

      if(grid->etop[j]==grid->etopold[j]-1) {
        phys->u[j][grid->etop[j]]=phys->u[j][grid->etopold[j]];
      } else { 
        for(k=grid->etop[j];k<=grid->etopold[j];k++)
          phys->u[j][k]=phys->u[j][grid->etopold[j]+1];
      }
    }
  }
}

static void WSourceTerms(gridT *grid, physT *phys, propT *prop,
  int myproc, int numprocs, MPI_Comm comm) {
int i, ib, iptr, j, jptr, k, ne, nf, nc, nc1, nc2, kmin, boundary_index;
REAL *a, *b, *c, Cz;

a = phys->a;
b = phys->b;
c = phys->c;

//read old dzz 
REAL **dzz_old2, **dzf_old, **dzf_old2;

char str[BUFFERLENGTH], filename[BUFFERLENGTH];
FILE *fid;


  //allocate space dzzold2, dzfold and old2
  dzz_old2 = (REAL **)malloc(grid->Nc*sizeof(REAL *));
  dzf_old = (REAL **)malloc(grid->Ne*sizeof(REAL *));
  dzf_old2 = (REAL **)malloc(grid->Ne*sizeof(REAL *));
  for(i=0;i<grid->Nc;i++) {
    dzz_old2[i] = (REAL *)malloc(grid->Nk[i]*sizeof(REAL));
  } 
  for(j=0;j<grid->Ne;j++) {
    dzf_old[j] = (REAL *)malloc(grid->Nke[j]*sizeof(REAL));
    dzf_old2[j] = (REAL *)malloc(grid->Nke[j]*sizeof(REAL));
  } 
  

  //read in
  if(prop->vertcoord==4 && prop->readOldVelocity==1){
    MPI_GetFile(filename,DATAFILE,"dzz_t2_init_file","HorizontalSource",myproc);  
    sprintf(str,"%s.%d",filename,myproc);  
    fid = MPI_FOpen(str,"r","HorizontalSource",myproc);

    for(j=0;j<grid->Nc;j++) {
      fread(dzz_old2[j],sizeof(REAL),grid->Nkmax,fid);
    }
    fclose(fid);
  } else{
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
        dzz_old2[i][k] = grid->dzzold[i][k];
  }

  ISendRecvCellData3D(dzz_old2,grid,myproc,comm);

  //now get old dzfs
  SetOldFluxHeight(dzf_old, dzf_old2, dzz_old2, grid, phys, prop, comm, myproc);

  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      phys->Cn_W[i][k]=phys->Cn_W2[i][k]=0;

  //Cn_W here 
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=0;

  //calculate Cn_W here 
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2];

    for(k=kmin;k<grid->Nke[j];k++) {
      a[k]=.5*(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->w[nc2][k]-phys->w[nc1][k]+phys->w[nc2][k+1]-phys->w[nc1][k+1])*grid->df[j]/grid->dg[j];
      
      if(!prop->subgrid)
      {
        phys->stmp[nc1][k]-=a[k]/grid->Ac[nc1];
        phys->stmp[nc2][k]+=a[k]/grid->Ac[nc2];
      } else {
        phys->stmp[nc1][k]-=a[k]/subgrid->Acceff[nc1][k]/grid->dzz[nc1][k]*grid->dzf[j][k];
        phys->stmp[nc2][k]+=a[k]/subgrid->Acceff[nc2][k]/grid->dzz[nc2][k]*grid->dzf[j][k];        
      }
      //      phys->stmp[nc1][k]-=prop->nu_H*(phys->w[nc2][k]-phys->w[nc1][k])*grid->df[j]/grid->dg[j]/grid->Ac[nc1];
      //      phys->stmp[nc2][k]-=prop->nu_H*(phys->w[nc1][k]-phys->w[nc2][k])*grid->df[j]/grid->dg[j]/grid->Ac[nc2];
    }
    for(k=grid->Nke[j];k<grid->Nk[nc1];k++) 
      if(!prop->subgrid)
        phys->stmp[nc1][k]+=0.25*prop->CdW*fabs(phys->w[nc1][k]+phys->w[nc1][k+1])*
        (phys->w[nc1][k]+phys->w[nc1][k+1])*grid->df[j]/grid->Ac[nc1];
      else
        phys->stmp[nc1][k]+=0.25*prop->CdW*fabs(phys->w[nc1][k]+phys->w[nc1][k+1])*
        (phys->w[nc1][k]+phys->w[nc1][k+1])*grid->df[j]*grid->dzf[j][k]/grid->dzz[nc1][k]/subgrid->Acceff[nc1][k];
    for(k=grid->Nke[j];k<grid->Nk[nc2];k++) 
      if(!prop->subgrid)
        phys->stmp[nc2][k]+=0.25*prop->CdW*fabs(phys->w[nc2][k]+phys->w[nc2][k+1])*
        (phys->w[nc2][k]+phys->w[nc2][k+1])*grid->df[j]/grid->Ac[nc2];
      else
        phys->stmp[nc2][k]+=0.25*prop->CdW*fabs(phys->w[nc2][k]+phys->w[nc2][k+1])*
        (phys->w[nc2][k]+phys->w[nc2][k+1])*grid->df[j]*grid->dzf[j][k]/grid->dzz[nc2][k]/subgrid->Acceff[nc2][k];      
  }

  // do the same for type 4 boundary conditions but utilize no-slip boundary condition
  for (jptr = grid->edgedist[4]; jptr < grid->edgedist[5]; jptr++){
    // get index for edge pointers
    j = grid->edgep[jptr];
    ib=grid->grad[2*j];
    boundary_index = jptr-grid->edgedist[2];

    // get neighbor indices
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // check to see which of the neighboring cells is the ghost cell
    if (nc1 == -1)  // indicating boundary
      nc = nc2;
    else
      nc = nc1;

    // loop over the entire depth
    for (k=grid->ctop[nc]; k<grid->Nke[nc]; k++){
      if(!prop->subgrid)
        phys->stmp[nc][k]  += -2.0*prop->nu_H*
          (phys->boundary_w[boundary_index][k] - 0.5*(phys->w[nc][k] + phys->w[nc][k+1]))/grid->dg[j]*
          grid->df[j]/grid->Ac[nc];
      else
        phys->stmp[nc][k]  += -2.0*prop->nu_H*
          (phys->boundary_w[boundary_index][k] - 0.5*(phys->w[nc][k] + phys->w[nc][k+1]))/grid->dg[j]*
          grid->df[j]/subgrid->Acceff[nc][k]*grid->dzf[j][k]/grid->dzz[nc][k];
    }
  }


  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 

    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      phys->Cn_W[i][k]-=prop->dt*(grid->dzz[i][k-1]*phys->stmp[i][k-1]+grid->dzz[i][k]*phys->stmp[i][k])/
        (grid->dzz[i][k-1]+grid->dzz[i][k]);
    }

    // Top flux advection consists only of top cell
    k=grid->ctop[i];
    phys->Cn_W[i][k]-=prop->dt*phys->stmp[i][k];

  }




//Cn_W2 here
for(i=0;i<grid->Nc;i++)
  for(k=0;k<grid->Nk[i];k++) 
    phys->stmp[i][k]=0;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2];

    for(k=kmin;k<grid->Nke[j];k++) {
      a[k]=.5*(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->w_old[nc2][k]-phys->w_old[nc1][k]+phys->w_old[nc2][k+1]-phys->w_old[nc1][k+1])*grid->df[j]/grid->dg[j];
      
      if(!prop->subgrid)
      {
        phys->stmp[nc1][k]-=a[k]/grid->Ac[nc1];
        phys->stmp[nc2][k]+=a[k]/grid->Ac[nc2];
      } else {
        phys->stmp[nc1][k]-=a[k]/subgrid->Acceff[nc1][k]/grid->dzzold[nc1][k]*dzf_old[j][k];
        phys->stmp[nc2][k]+=a[k]/subgrid->Acceff[nc2][k]/grid->dzzold[nc2][k]*dzf_old[j][k];        
      }
      //      phys->stmp[nc1][k]-=prop->nu_H*(phys->w[nc2][k]-phys->w[nc1][k])*grid->df[j]/grid->dg[j]/grid->Ac[nc1];
      //      phys->stmp[nc2][k]-=prop->nu_H*(phys->w[nc1][k]-phys->w[nc2][k])*grid->df[j]/grid->dg[j]/grid->Ac[nc2];
    }
    for(k=grid->Nke[j];k<grid->Nk[nc1];k++) 
      if(!prop->subgrid)
        phys->stmp[nc1][k]+=0.25*prop->CdW*fabs(phys->w_old[nc1][k]+phys->w_old[nc1][k+1])*
        (phys->w_old[nc1][k]+phys->w_old[nc1][k+1])*grid->df[j]/grid->Ac[nc1];
      else
        phys->stmp[nc1][k]+=0.25*prop->CdW*fabs(phys->w_old[nc1][k]+phys->w_old[nc1][k+1])*
        (phys->w_old[nc1][k]+phys->w_old[nc1][k+1])*grid->df[j]*dzf_old[j][k]/grid->dzzold[nc1][k]/subgrid->Acceff[nc1][k];
    for(k=grid->Nke[j];k<grid->Nk[nc2];k++) 
      if(!prop->subgrid)
        phys->stmp[nc2][k]+=0.25*prop->CdW*fabs(phys->w_old[nc2][k]+phys->w_old[nc2][k+1])*
        (phys->w_old[nc2][k]+phys->w_old[nc2][k+1])*grid->df[j]/grid->Ac[nc2];
      else
        phys->stmp[nc2][k]+=0.25*prop->CdW*fabs(phys->w_old[nc2][k]+phys->w_old[nc2][k+1])*
        (phys->w_old[nc2][k]+phys->w_old[nc2][k+1])*grid->df[j]*dzf_old[j][k]/grid->dzzold[nc2][k]/subgrid->Acceff[nc2][k];      
  }

  // do the same for type 4 boundary conditions but utilize no-slip boundary condition
  for (jptr = grid->edgedist[4]; jptr < grid->edgedist[5]; jptr++){
    // get index for edge pointers
    j = grid->edgep[jptr];
    ib=grid->grad[2*j];
    boundary_index = jptr-grid->edgedist[2];

    // get neighbor indices
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // check to see which of the neighboring cells is the ghost cell
    if (nc1 == -1)  // indicating boundary
      nc = nc2;
    else
      nc = nc1;

    // loop over the entire depth
    for (k=grid->ctop[nc]; k<grid->Nke[nc]; k++){
      if(!prop->subgrid)
        phys->stmp[nc][k]  += -2.0*prop->nu_H*
          (phys->boundary_w[boundary_index][k] - 0.5*(phys->w_old[nc][k] + phys->w_old[nc][k+1]))/grid->dg[j]*
          grid->df[j]/grid->Ac[nc];
      else
        phys->stmp[nc][k]  += -2.0*prop->nu_H*
          (phys->boundary_w[boundary_index][k] - 0.5*(phys->w_old[nc][k] + phys->w_old[nc][k+1]))/grid->dg[j]*
          grid->df[j]/subgrid->Acceff[nc][k]*dzf_old[j][k]/grid->dzzold[nc][k];
    }
  }


  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 

    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
      phys->Cn_W2[i][k]-=prop->dt*(grid->dzzold[i][k-1]*phys->stmp[i][k-1]+grid->dzzold[i][k]*phys->stmp[i][k])/
        (grid->dzzold[i][k-1]+grid->dzzold[i][k]);

    // Top flux diffusion consists only of top cell
    k=grid->ctop[i];
    phys->Cn_W2[i][k]-=prop->dt*phys->stmp[i][k];

  }

  //   for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
  //     phys->Cn_W2[i][k]-=prop->dt*(grid->dzz[i][k-1]*phys->stmp[i][k-1]+grid->dzz[i][k]*phys->stmp[i][k])/
  //       (grid->dzz[i][k-1]+grid->dzz[i][k]);
  //   }

  //   // Top flux advection consists only of top cell
  //   k=grid->ctop[i];
  //   phys->Cn_W2[i][k]-=prop->dt*phys->stmp[i][k];

  // }

  for(i=0;i<grid->Nc;i++) {
    free(dzz_old2[i]);
  }
  for(j=0;j<grid->Ne;j++) {
    free(dzf_old[j]);
    free(dzf_old2[j]);
  }
  free(dzf_old);
  free(dzf_old2);
  free(dzz_old2);

}

/*
 * Function: WPredictor
 * Usage: WPredictor(grid,phys,prop,myproc,numprocs,comm);
 * -------------------------------------------------------
 * This function updates the vertical predicted velocity field with:
 *
 * 1) Old nonhydrostatic pressure gradient with theta method
 * 2) Horizontal and vertical advection of vertical momentum with AB2
 * 3) Horizontal laminar+turbulent diffusion of vertical momentum with AB2
 * 4) Vertical laminar+turbulent diffusion of vertical momentum with theta method
 *
 * Cn_W contains the Adams-Bashforth terms at time step n-1.
 * If wetting and drying is employed, no advection is computed in 
 * the upper cell.
 *
 */
static void WPredictor(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, MPI_Comm comm) {
  int i, ib, iptr, j, jptr, k, ne, nf, nc, nc1, nc2, kmin, boundary_index, nonlinear=2, Wadv;
  REAL fab1,fab2,fab3, fac1,fac2,fac3, sum, *a, *b, *c, *d, *e, *r, *a0, *b0, *c0, *d0, *e0,
    Cz, u_im, u_im_km1, Q_im, *omega_eff;
  REAL w0 = 1.0;
  REAL d_old, div, divmax=-INFTY, err_max=-INFTY;

  REAL *J_eff_new = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"WPredictor");
  REAL *J_eff_old = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"WPredictor");

  // Momentum advection scheme selector.  nonlinear==0 in suntans.dat selects the new
  // scheme (Wadv follows W_MOMENTUM_ADVECTION); nonlinear!=0 selects the stock Eulerian
  // advection below, in which case every W_MOMENTUM_ADVECTION block must stand down so
  // that w is not advected twice.
  if(prop->nonlinear)
    Wadv = -1;
  else
    Wadv = W_MOMENTUM_ADVECTION;

  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;
  e = phys->e;
  r = phys->pent_source;
  omega_eff = phys->pent_a;

  a0 = phys->pent_b;
  b0 = phys->bp;
  c0 = phys->bm;  
  d0 = phys->pent_d;
  e0 = phys->pent_e; 
  REAL *ap = phys->ap;
  REAL *am = phys->am;
  
  // AB3
  if(!prop->readOldVelocity){
    //if(1){
  if(prop->n==1) {
    fab1=1;
    fab2=fab3=0;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
        phys->Cn_W[i][k]=phys->Cn_W2[i][k]=0;
  } else if(prop->n==2) {
    fab1=3.0/2.0;
    fab2=-1.0/2.0;
    fab3=0;
    }  else {
      fab1=prop->exfac1;
      fab2=prop->exfac2;
      fab3=prop->exfac3;   
    }
  } else{
    fab1=prop->exfac1;
    fab2=prop->exfac2;
    fab3=prop->exfac3;   
  }

  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++) 
      phys->stmp[i][k]=phys->stmp2[i][k]=0;

  if(Wadv==1 && prop->n==prop->nsteps+prop->nstart && CHECKCONSISTENCY) {
    prop->nu=prop->nu_H=0;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++) {
	phys->w_ex[i][k]=phys->w[i][k]=phys->w_old[i][k]=phys->w_old2[i][k]=w0;
	phys->Cn_W[i][k]=phys->Cn_W2[i][k]=phys->q[i][k]=0;
      }
  }

  // Interpolate w to faces and place into ut
  if(Wadv==1) {
    nonlinear=5;

    ISendRecvWData(phys->w_ex,grid,myproc,comm);   
    ISendRecvEdgeData3D(phys->utmp,grid,myproc,comm);   
    //was u_old 
    GetMomentumFaceValues(phys->ut,phys->w_ex,phys->boundary_w,phys->utmp,
			  grid,phys,prop,comm,myproc,nonlinear,prop->TVDmomentum);
  }
  if(DEBUG) printf("Debugging - Wpred, w_ex advected to edges...\n");
  
  // Compute horizontal advection and diffusion
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 

    if(Wadv==1) {
      // J_eff[k] stores 0.5*(J_{k-1} + J_{k}) except at top and bottom
      // where it stores J_{k}
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	J_eff_new[k] = 0.5*(grid->dzz[i][k-1]+grid->dzz[i][k]);
	J_eff_old[k] = 0.5*(grid->dzzold[i][k-1]+grid->dzzold[i][k]);
      }
      J_eff_new[grid->ctop[i]]=grid->dzz[i][grid->ctop[i]];
      J_eff_old[grid->ctop[i]]=grid->dzzold[i][grid->ctop[i]];    

      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	//phys->wtmp[i][k]=J_eff_old[k]/J_eff_new[k]*phys->w[i][k] +
	// (0.5*(J_eff_old[k]+J_eff_new[k]))/J_eff_new[k]*(fab2*phys->Cn_W[i][k] + fab3*phys->Cn_W2[i][k]);
	phys->wtmp[i][k]=J_eff_old[k]/J_eff_new[k]*phys->w[i][k] +
	fab2*phys->Cn_W[i][k] + fab3*phys->Cn_W2[i][k];
	phys->Cn_W2[i][k]=phys->Cn_W[i][k];
	phys->Cn_W[i][k]=0;
      }
    } else {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	phys->wtmp[i][k]=phys->w[i][k] +
	  fab2*phys->Cn_W[i][k] + fab3*phys->Cn_W2[i][k];
	phys->Cn_W2[i][k]=phys->Cn_W[i][k];
	phys->Cn_W[i][k]=0;
      }
    }

    // Only include old pressure gradient for correction method
    if(prop->pressureMethod==1) {    
      if(0){
	  for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
	    phys->wtmp[i][k]-=2.0*prop->dt/(grid->dzz[i][k-1]+grid->dzz[i][k])*
	      (phys->q[i][k-1]-phys->q[i][k])*(J_eff_old[k])/J_eff_new[k];
	  phys->wtmp[i][grid->ctop[i]]+=2.0*prop->dt/grid->dzz[i][grid->ctop[i]]*
	    phys->q[i][grid->ctop[i]]*(J_eff_old[k])/J_eff_new[k];
	}else{
	  for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
	    phys->wtmp[i][k]-=2.0*prop->dt/(grid->dzz[i][k-1]+grid->dzz[i][k])*
	  (phys->q[i][k-1]-phys->q[i][k]);
      phys->wtmp[i][grid->ctop[i]]+=2.0*prop->dt/grid->dzz[i][grid->ctop[i]]*
	phys->q[i][grid->ctop[i]];
	}
    }

    if(Wadv==1) {
      // Add horizontal advection of w and store it in stmp2
      for(nf=0;nf<grid->nfaces[i];nf++) {
	ne = grid->face[i*grid->maxfaces+nf];
	
	for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	  
	  // The top cell advection is given by momentum advection with a control volume equal
	  // to the top cell, not averaged over the cells above and below the w-face.
	  
	  // u_im = fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k];
    u_im = fac1*phys->utmp[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k];
	
	  if(k==grid->ctop[i]) {
	    Q_im = u_im*grid->dzf[ne][k];
	  } else {
      u_im_km1 = fac1*phys->utmp[ne][k-1]+fac2*phys->u_old[ne][k-1]+fac3*phys->u_old2[ne][k-1];	  
	    //u_im_km1 = fac1*phys->u[ne][k-1]+fac2*phys->u_old[ne][k-1]+fac3*phys->u_old2[ne][k-1];	  
	    Q_im = 0.5*(u_im_km1*grid->dzf[ne][k-1]+u_im*grid->dzf[ne][k]);
	  }
	  
	  phys->stmp2[i][k]-=
	    prop->dt/(grid->Ac[i]*J_eff_new[k])*
	    phys->ut[ne][k]*Q_im*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];
	}
      }
    }
  }
  if(DEBUG) printf("Debugging - Wpred, hoirz advection to stmp2...\n");

  // Stock (upstream) Eulerian advection of w, active when nonlinear!=0.
  // Accumulates into stmp, which is folded into Cn_W below alongside horizontal diffusion.
  if(prop->nonlinear && (prop->vertcoord==1 || (prop->vertcoord!=1 && !prop->wetdry))) {
    // First compute w at the cell centers (since w is defined at the faces)
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        phys->wc[i][k]=0.5*(phys->w[i][k]+phys->w[i][k+1]);
    }

    // Interpolate wc to faces and place into ut
    GetMomentumFaceValues(phys->ut,phys->wc,phys->boundary_w,phys->u_old,grid,phys,prop,comm,myproc,prop->nonlinear,prop->TVDmomentum);

    if(prop->conserveMomentum)
      for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
        j=grid->edgep[jptr];

        for(k=grid->etop[j];k<grid->Nke[j];k++)
          phys->ut[j][k]*=grid->dzf[j][k];
      }

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i=grid->cellp[iptr];

      // For conservative scheme need to divide by depth (since ut is a flux)
      if(prop->conserveMomentum) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          a[k]=grid->dzz[i][k];
      } else {
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          a[k]=1.0;
      }

      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          if(!prop->subgrid || prop->wetdry)
            phys->stmp[i][k]+=phys->ut[ne][k]*phys->u_old[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/(a[k]*grid->Ac[i]);
          else
            phys->stmp[i][k]+=phys->ut[ne][k]*phys->u_old[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/(a[k]*subgrid->Acceff[i][k]);

        // Top cell is filled with momentum from neighboring cells
        if(prop->conserveMomentum){
          for(k=grid->etop[ne];k<grid->ctop[i];k++)
            if(!prop->subgrid || prop->wetdry)
              phys->stmp[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->u_old[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/(a[k]*grid->Ac[i]);
            else
              phys->stmp[i][grid->ctop[i]]+=phys->ut[ne][k]*phys->u_old[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/(a[k]*subgrid->Acceff[i][k]);
        }
      }

      // Vertical advection; note that in this formulation first-order upwinding is not implemented.
      // if not z-level, the vertical momentum advection part is added later
      if(prop->nonlinear==1 || prop->nonlinear==2 || prop->nonlinear==5) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          if(prop->vertcoord==1)
          {
            if(!prop->subgrid || prop->wetdry)
              phys->stmp[i][k]+=(pow(phys->w[i][k],2)-pow(phys->w[i][k+1],2))/grid->dzz[i][k];
            else
              phys->stmp[i][k]+=(pow(phys->w[i][k],2)*subgrid->Acveffold[i][k]-pow(phys->w[i][k+1],2)*subgrid->Acveffold[i][k+1])/grid->dzz[i][k]/subgrid->Acceff[i][k];
          }else{
            if(!prop->subgrid)
              phys->stmp[i][k]+=(vert->omega_old[i][k]*phys->w[i][k]-vert->omega_old[i][k+1]*phys->w[i][k+1])/grid->dzz[i][k];
            else
              phys->stmp[i][k]+=(vert->omega_old[i][k]*phys->w[i][k]*subgrid->Acveffold[i][k]-
                vert->omega_old[i][k+1]*phys->w[i][k+1]*subgrid->Acveffold[i][k+1])/grid->dzz[i][k]/subgrid->Acceff[i][k];
          }
        }
      }
    }

    // Check to make sure integrated fluxes are 0 for conservation
    if(WARNING && prop->CdW==0 && prop->nu_H==0) {
      sum=0;
      for(i=0;i<grid->Nc;i++) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          if(!prop->subgrid || prop->wetdry)
            sum+=grid->Ac[i]*phys->stmp[i][k]*grid->dzz[i][k];
          else
            sum+=subgrid->Acceff[i][k]*phys->stmp[i][k]*grid->dzz[i][k];
      }
      if(fabs(sum)>CONSERVED)
        printf("Warning, not W-momentum conservative!\n");
    }
  }

  // Compute horizontal diffusion of cell-centered w and put it into phys->stmp[][]
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2];

    for(k=kmin;k<grid->Nke[j];k++) {
      a[k]=.5*(prop->nu_H+0.5*(phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]))*
        (phys->w[nc2][k]-phys->w[nc1][k]+phys->w[nc2][k+1]-phys->w[nc1][k+1])*grid->df[j]/grid->dg[j];
      
      phys->stmp[nc1][k]-=a[k]/grid->Ac[nc1];
      phys->stmp[nc2][k]+=a[k]/grid->Ac[nc2];
    }

    for(k=grid->Nke[j];k<grid->Nk[nc1];k++) {
      phys->stmp[nc1][k]+=0.25*prop->CdW*fabs(phys->w[nc1][k]+phys->w[nc1][k+1])*
	(phys->w[nc1][k]+phys->w[nc1][k+1])*grid->df[j]/grid->Ac[nc1];
    }

    for(k=grid->Nke[j];k<grid->Nk[nc2];k++) {
      phys->stmp[nc2][k]+=0.25*prop->CdW*fabs(phys->w[nc2][k]+phys->w[nc2][k+1])*
	(phys->w[nc2][k]+phys->w[nc2][k+1])*grid->df[j]/grid->Ac[nc2];
    }
  }

  // Add horizontal diffusion from the type-four boundaries (no-slip with velocity specified)
  for (jptr = grid->edgedist[4]; jptr < grid->edgedist[5]; jptr++){
    // get index for edge pointers
    j = grid->edgep[jptr];
    ib=grid->grad[2*j];
    boundary_index = jptr-grid->edgedist[2];

    // get neighbor indices
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // check to see which of the neighboring cells is the ghost cell
    if (nc1 == -1)  // indicating boundary
      nc = nc2;
    else
      nc = nc1;

    // loop over the entire depth
    for (k=grid->ctop[nc]; k<grid->Nke[nc]; k++){
        phys->stmp[nc][k]  += -2.0*prop->nu_H*
          (phys->boundary_w[boundary_index][k] - 0.5*(phys->w[nc][k] + phys->w[nc][k+1]))/grid->dg[j]*
          grid->df[j]/grid->Ac[nc];
    }
  }

  // Now use the cell-centered horizontal diffusion terms to update the diffusion at the faces
  // Then add horizontal diffusion and advection to wtmp
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 

    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
      phys->Cn_W[i][k]-=prop->dt*(grid->dzz[i][k-1]*phys->stmp[i][k-1]+grid->dzz[i][k]*phys->stmp[i][k])/
        (grid->dzz[i][k-1]+grid->dzz[i][k]);

    // Top flux diffusion consists only of top cell
    k=grid->ctop[i];
    phys->Cn_W[i][k]-=prop->dt*phys->stmp[i][k];

    // Stock: additional part for the new vertical coordinate (nonlinear!=0 only)
    if(prop->vertcoord!=1 && prop->nonlinear)
    {
      if(!prop->wetdry){
        if(vert->dJdtmeth==1)
        {
          for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
            phys->Cn_W[i][k]-=phys->w[i][k]*(grid->dzz[i][k]*(1-grid->dzzold[i][k-1]/grid->dzz[i][k-1])+
                grid->dzz[i][k-1]*(1-grid->dzzold[i][k]/grid->dzz[i][k]))/(grid->dzz[i][k]+grid->dzz[i][k-1]);
          k=grid->ctop[i];
          phys->Cn_W[i][k]-=phys->w[i][k]*(1-grid->dzzold[i][k]/grid->dzz[i][k]);
        }
      } else {
        for(k=grid->ctop[i]+1;k<grid->Nk[i];k++){
          // subtract the addition part of horizontal advection
          phys->Cn_W[i][k]-=prop->dt*(InterpToLayerTopFace(i,k,phys->uc,grid)*
                vert->dwdx[i][k]+InterpToLayerTopFace(i,k,phys->vc,grid)*vert->dwdy[i][k]);

          phys->Cn_W[i][k]-=prop->dt*vert->omega_old[i][k]*(phys->w[i][k-1]-phys->w[i][k+1])/(grid->dzz[i][k]+grid->dzz[i][k-1]);
        }

        k=grid->ctop[i];
        phys->Cn_W[i][k]-=prop->dt*(InterpToLayerTopFace(i,k,phys->uc,grid)*vert->dwdx[i][k]+
                  InterpToLayerTopFace(i,k,phys->vc,grid)*vert->dwdy[i][k]);
        phys->Cn_W[i][k]-=prop->dt*vert->omega_old[i][k]*(phys->w[i][k]-phys->w[i][k+1])/grid->dzz[i][k];
      }
    }
  }

  // Stock: vertical advection using Lax-Wendroff
  if(prop->nonlinear==4 && (prop->vertcoord==1 || (prop->vertcoord!=1 && !prop->wetdry)))
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      if(prop->vertcoord==1)
        for(k=grid->ctop[i]+1;k<grid->Nk[i]+1;k++) {
          Cz = 0.5*(phys->w[i][k-1]+phys->w[i][k])*prop->dt/grid->dzz[i][k-1];
          a[k]=0.5*(phys->w[i][k-1]+phys->w[i][k])*(0.5*(phys->w[i][k-1]+phys->w[i][k])-0.5*Cz*(phys->w[i][k-1]-phys->w[i][k]));
        }
      else
        for(k=grid->ctop[i]+1;k<grid->Nk[i]+1;k++) {
          Cz = 0.5*(vert->omega_old[i][k-1]+vert->omega_old[i][k])*prop->dt/grid->dzz[i][k-1];
          a[k]=0.5*(vert->omega_old[i][k-1]+vert->omega_old[i][k])*(0.5*(phys->w[i][k-1]+phys->w[i][k])-0.5*Cz*(phys->w[i][k-1]-phys->w[i][k]));
        }
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
        phys->Cn_W[i][k]-=2.0*prop->dt*(a[k]-a[k+1])/(grid->dzz[i][k]+grid->dzz[i][k+1]);
      }
    }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    if(0) {
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	J_eff_new[k] = 0.5*(grid->dzz[i][k-1]+grid->dzz[i][k]);
	J_eff_old[k] = 0.5*(grid->dzzold[i][k-1]+grid->dzzold[i][k]);
      }
      J_eff_new[grid->ctop[i]]=grid->dzz[i][grid->ctop[i]];
      J_eff_old[grid->ctop[i]]=grid->dzzold[i][grid->ctop[i]];
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	phys->wtmp[i][k]+=(fab1*phys->Cn_W[i][k]*(0.5*(J_eff_old[k]+J_eff_new[k]))/J_eff_new[k]+phys->stmp2[i][k]);
    }else{
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	phys->wtmp[i][k]+=(fab1*phys->Cn_W[i][k]+phys->stmp2[i][k]); 
    }
  }
  if(DEBUG) printf("Debugging - Wpred, source terms handled...\n");

  // wtmp now contains the right hand side without the vertical diffusion terms.  Now we
  // add the vertical advection diffusion terms to the explicit side and invert the tridiagonal for
  // the implicit terms.
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 

    // Vertical diffusion
    if(grid->Nk[i]-grid->ctop[i]>1) {
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) { // multiple layers
        a[k] = 2*(prop->nu+prop->laxWendroff_Vertical*phys->nu_lax[i][k-1]+
            phys->nu_tv[i][k-1])/grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
        b[k] = 2*(prop->nu+prop->laxWendroff_Vertical*phys->nu_lax[i][k]+
            phys->nu_tv[i][k])/grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
      }
      b[grid->ctop[i]]=(prop->nu+prop->laxWendroff_Vertical*phys->nu_lax[i][grid->ctop[i]]+
          phys->nu_tv[i][grid->ctop[i]])/pow(grid->dzz[i][grid->ctop[i]],2);
      a[grid->ctop[i]]=b[grid->ctop[i]];

      // Add on the explicit part of the vertical diffusion term
      // add the new implicit scheme
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
        phys->wtmp[i][k]+=prop->dt*(a[k]*(fac2*phys->w_old[i][k-1]+fac3*phys->w_old2[i][k-1])
            -(a[k]+b[k])*(fac2*phys->w_old[i][k]+fac3*phys->w_old2[i][k])
            +b[k]*(fac2*phys->w_old[i][k+1]+fac3*phys->w_old2[i][k+1]));
      phys->wtmp[i][grid->ctop[i]]+=prop->dt*(-(a[grid->ctop[i]]+b[grid->ctop[i]])*(fac2*phys->w_old[i][grid->ctop[i]]+fac3*phys->w_old2[i][grid->ctop[i]])
          +(a[grid->ctop[i]]+b[grid->ctop[i]])*(fac2*phys->w_old[i][grid->ctop[i]+1]+fac3*phys->w_old2[i][grid->ctop[i]+1]));

      // Now formulate the components of the tridiagonal inversion.
      // c is the diagonal entry, a is the lower diagonal, and b is the upper diagonal.
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        c[k]=1+prop->dt*fac1*(a[k]+b[k]);
        a[k]*=(-prop->dt*fac1);
        b[k]*=(-prop->dt*fac1);
      }
      b[grid->ctop[i]]+=a[grid->ctop[i]];

      // Vertical advection
      if(Wadv==1) {
	// J_eff[k] stores 0.5*(J_{k-1} + J_{k}) 
	for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	  J_eff_new[k] = 0.5*(grid->dzz[i][k-1]+grid->dzz[i][k]);
	  J_eff_old[k] = 0.5*(grid->dzzold[i][k-1]+grid->dzzold[i][k]);
	}
	J_eff_new[grid->ctop[i]]=grid->dzz[i][grid->ctop[i]];
	J_eff_old[grid->ctop[i]]=grid->dzzold[i][grid->ctop[i]];    	
	if(DEBUG) printf("Debugging - Wpred, Jacobains set ...\n");

	for(k=grid->ctop[i];k<grid->Nk[i]-1;k++) 
	  omega_eff[k] = 0.5*(vert->omega_im[i][k]+vert->omega_im[i][k+1]);
	omega_eff[grid->Nk[i]-1]=vert->omega_im[i][grid->Nk[i]-1]+0;
	if(DEBUG) printf("Debugging - omega eff set...\n");

  if(nonlinear == 2){
    //Coefficients based on central differencing
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
      a0[k] = -0.5/J_eff_new[k]*omega_eff[k-1];
      b0[k] = -0.5/J_eff_new[k]*(omega_eff[k-1]-omega_eff[k]);
      c0[k] = 0.5/J_eff_new[k]*omega_eff[k];
    }
    c0[grid->Nk[i]-1]=0;
    k=grid->ctop[i];
    a0[k] = 0;
    b0[k] = -1/J_eff_new[k]*vert->omega_im[i][k];
    c0[k] = 1/J_eff_new[k]*vert->omega_im[i][k+1];
  }else if(nonlinear==5 && prop->TVDmomentum>=5){
  //coefficients for quick and sharp
    GetPentDiagVert(a0, b0, c0, d0, e0, phys->wp,phys->wm,
      omega_eff, grid->dzz,phys->w_ex, i,grid->Nk[i],grid->ctop[i],prop->dt, prop->TVDmomentum,1,1);

    k=grid->ctop[i];
    a0[k] = 0;
    b0[k] = 0;
    c0[k] = -1*vert->omega_im[i][k];
    d0[k] = 1*vert->omega_im[i][k+1];
    e0[k] = 0;

  }else{
    //default to normal TVD
    GetApAmVert(ap, am, phys->wp,phys->wm, phys->Cp, phys->Cm, phys->rp, phys->rm,
      omega_eff, grid->dzz,phys->w_ex, i,grid->Nk[i],grid->ctop[i],prop->dt, prop->TVDmomentum,1,1);

      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++) {
        a0[k] = -1/J_eff_new[k]*am[k-1];
        b0[k] = -1/J_eff_new[k]*(ap[k-1]-am[k]);
        c0[k] = 1/J_eff_new[k]*ap[k];
      }

      c0[grid->Nk[i]-1]=0;
      k=grid->ctop[i];
      a0[k] = 0;
      b0[k] = -1/J_eff_new[k]*vert->omega_im[i][k];
      c0[k] = 1/J_eff_new[k]*vert->omega_im[i][k+1];

  }
      
  if(DEBUG) printf("Debugging - Wpred, got vertical coefficients...\n");
	// Explicit part of vertical advection

  if(nonlinear==2 || prop->TVDmomentum<5){
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      phys->wtmp[i][k]+=prop->dt*(a0[k]*(fac2*phys->w_old[i][k-1]+fac3*phys->w_old2[i][k-1])+
                b0[k]*(fac2*phys->w_old[i][k]+fac3*phys->w_old2[i][k])+
                c0[k]*(fac2*phys->w_old[i][k+1]+fac3*phys->w_old2[i][k+1]));
    }

    // Implicit part of vertical advection
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      a[k]-=prop->dt*fac1*a0[k];
      c[k]-=prop->dt*fac1*b0[k];
      b[k]-=prop->dt*fac1*c0[k];
    }
  } else{
    for(k=grid->ctop[i]+2;k<grid->Nk[i]-1;k++) {
      phys->wtmp[i][k]+=prop->dt/J_eff_new[k]*(a0[k]*(fac2*phys->w_old[i][k-2]+fac3*phys->w_old2[i][k-2])+
                b0[k]*(fac2*phys->w_old[i][k-1]+fac3*phys->w_old2[i][k-1])+
                c0[k]*(fac2*phys->w_old[i][k]+fac3*phys->w_old2[i][k])+
                d0[k]*(fac2*phys->w_old[i][k+1]+fac3*phys->w_old2[i][k+1])+
                e0[k]*(fac2*phys->w_old[i][k+2]+fac3*phys->w_old2[i][k+2]));
    }

    //lets try not calling any non existent ws even if its multiplied by a 0
    k=grid->ctop[i]; //there is w[0] but there is no w[-1] or w[-2]
    phys->wtmp[i][k]+=prop->dt/J_eff_new[k]*(c0[k]*(fac2*phys->w_old[i][k]+fac3*phys->w_old2[i][k])+
					     d0[k]*(fac2*phys->w_old[i][k+1]+fac3*phys->w_old2[i][k+1])+
					     e0[k]*(fac2*phys->w_old[i][k+2]+fac3*phys->w_old2[i][k+2]));

    k=grid->ctop[i]+1; //there is w[0] but there is no w[-1]
    phys->wtmp[i][k]+=prop->dt/J_eff_new[k]*(b0[k]*(fac2*phys->w_old[i][k-1]+fac3*phys->w_old2[i][k-1])+
					     c0[k]*(fac2*phys->w_old[i][k]+fac3*phys->w_old2[i][k])+
					     d0[k]*(fac2*phys->w_old[i][k+1]+fac3*phys->w_old2[i][k+1])+
					     e0[k]*(fac2*phys->w_old[i][k+2]+fac3*phys->w_old2[i][k+2]));

    k=grid->Nk[i]-1; //there is w[Nk] but no w[Nk+1]. Also w[Nk+1] needs to be 0, so we can just ignore that one all together
    phys->wtmp[i][k]+=prop->dt/J_eff_new[k]*(a0[k]*(fac2*phys->w_old[i][k-2]+fac3*phys->w_old2[i][k-2])+
					     b0[k]*(fac2*phys->w_old[i][k-1]+fac3*phys->w_old2[i][k-1])+
					     c0[k]*(fac2*phys->w_old[i][k]+fac3*phys->w_old2[i][k]));

    if(DEBUG) printf("Debugging - Wpred, explicit vertical through boundaries...\n");

    // Implicit part of vertical advection
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      e[k]=b[k]; //place holder
      b[k]=a[k]; //upper diagonal
      d[k]=e[k]; //lower diagonal
      e[k] = 0; //lower lower diag
      a[k] = 0; //upper upper diag


      a[k]-=prop->dt/J_eff_new[k]*fac1*a0[k];
      b[k]-=prop->dt/J_eff_new[k]*fac1*b0[k];
      c[k]-=prop->dt/J_eff_new[k]*fac1*c0[k];
      d[k]-=prop->dt/J_eff_new[k]*fac1*d0[k];
      e[k]-=prop->dt/J_eff_new[k]*fac1*e0[k];
      if(DEBUG) printf("Debugging - Wpred, set penta coeffs...\n");
    }
  }
  }
      
      // Pentadiagonal coefficients are only built in the Wadv==1 block, so any other
      // mode (stock advection, or Wadv==0) must fall through to the tridiagonal solve.
      if(Wadv!=1 || nonlinear==2 || prop->TVDmomentum<5){
      TriSolve(&(a[grid->ctop[i]]),&(c[grid->ctop[i]]),&(b[grid->ctop[i]]),
          &(phys->wtmp[i][grid->ctop[i]]),&(phys->w[i][grid->ctop[i]]),grid->Nk[i]-grid->ctop[i]);
      }else
        ReducePentadiag(&(a[grid->ctop[i]]),&(b[grid->ctop[i]]),&(c[grid->ctop[i]]),
        &(d[grid->ctop[i]]),&(e[grid->ctop[i]]),&(phys->wtmp[i][grid->ctop[i]]),&(phys->w[i][grid->ctop[i]]),grid->Nk[i]-grid->ctop[i]);
    } else { // one layer
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
        phys->w[i][k]=phys->wtmp[i][k];
    }
  }
  if(DEBUG) printf("Debugging - Wpred, through matrix inversion...\n");


  // Check whether the momentum advection scheme is CWC at the last time step
  // Only works if nu=nu_H=0
  if(Wadv==1 && prop->n==prop->nsteps+prop->nstart && CHECKCONSISTENCY) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	div=vert->omega_im[i][k]-vert->omega_im[i][k+1];
	for(nf=0;nf<grid->nfaces[i];nf++) {
	  ne = grid->face[i*grid->maxfaces+nf];
	  u_im = fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k];
	  
	  div+=1.0/grid->Ac[i]*
	    u_im*grid->dzf[ne][k]*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];
	}
	div=(grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt+div;
	if(fabs(div)>divmax)
	  divmax=fabs(div);
	  
	if(fabs((phys->w[i][k]-w0)/w0)>err_max)
	  err_max=fabs((phys->w[i][k]-w0)/w0);

	/*
	if(myproc==0 && i==floor(grid->Nc/2)) {
	  printf("k=%d, wnew=%.2e, |(wnew-w0)/w0|=%.4e, div=%.4e\n",
		 k,phys->w[i][k],fabs((phys->w[i][k]-w0)/w0),div);
	}
	*/
      }
    }

    if(Wadv==1 && prop->n==prop->nsteps+prop->nstart && CHECKCONSISTENCY) {
      printf("W adv: proc: %d, divmax = %.4e, err_max = %.4e\n",myproc,divmax,err_max);
    }
  }

  // W advection
  if(Wadv==0) {
    int i, j, iptr, jptr, k, nc1, nc2;
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      for(k=0;k<grid->Nk[i];k++) 
	phys->wtmp[i][k]=phys->stmp3[i][k]=phys->Cn_T[i][k]=0;
      phys->wtmp[i][grid->Nk[i]]=0;
    }
	
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=0;k<grid->ctopold[i];k++){
      	phys->vold[i][k]=phys->uold[i][k]=0;
        phys->uold2[i][k]=0;
      }
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        //first estimate with central differencing for sign?
	phys->vold[i][k]=phys->uold[i][k]=0.5*(phys->w_old[i][k]+phys->w_old[i][k+1]);
        phys->uold2[i][k]=0.5*(phys->w_old2[i][k]+phys->w_old2[i][k+1]);

                //instead, use quick. First do with uniform 
        if(k>0 && k<(grid->Nk[i]-1)){
          if(phys->vold[i][k]>0)
            phys->vold[i][k]=phys->uold[i][k]=0.5*(phys->w_old[i][k]+phys->w_old[i][k+1]) - (1/8)*(phys->w_old[i][k] - 2*phys->w_old[i][k+1]+ phys->w_old[i][k+2]);
          else if(phys->vold[i][k]<0)
            phys->vold[i][k]=phys->uold[i][k]=0.5*(phys->w_old[i][k]+phys->w_old[i][k+1]) - (1/8)*(phys->w_old[i][k+1] - 2*phys->w_old[i][k]+ phys->w_old[i][k-1]);

          if(phys->uold2[i][k]>0)
            phys->uold2[i][k]=0.5*(phys->w_old2[i][k]+phys->w_old2[i][k+1]) - (1/8)*(phys->w_old2[i][k] - 2*phys->w_old2[i][k+1]+ phys->w_old2[i][k+2]);
          else if(phys->uold2[i][k]<0)
            phys->uold2[i][k]=0.5*(phys->w_old2[i][k]+phys->w_old2[i][k+1]) - (1/8)*(phys->w_old2[i][k+1] - 2*phys->w_old2[i][k]+ phys->w_old2[i][k-1]);
        }

        //let be central if w=0 so it stays 0. also let it be at boundary 

        phys->uold2[i][k]=0.5*(phys->w[i][k]+phys->w[i][k+1]);
      }
      
    }
    ISendRecvCellData3D(phys->uold,grid,myproc,comm);
    ISendRecvCellData3D(phys->vold,grid,myproc,comm);

    ISendRecvCellData3D(phys->uold2,grid,myproc,comm);

    //not sure about the uold2's here, were s_old which also  isn't right. 
    if(prop->vertcoord==1) {
       UpdateScalars(grid,phys,prop,phys->w_im,phys->w_SfH_tm1,phys->w_SfH_tm2, phys->w_SfHv_t, phys->w_SfHv_tm1,phys->uold,phys->uold2,phys->boundary_w,2,2,phys->Cn_T,
		     0,0,NULL,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,1,1,prop->TVDmomentum);
      //UpdateScalars(grid,phys,prop,phys->w_im,phys->w_SfH_tm1,phys->w_SfH_tm2, phys->w_SfHv_t, phys->w_SfHv_tm1,phys->uold,phys->wc,phys->boundary_w,2,2,phys->Cn_T,
		  //  0,0,NULL,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,0,1,prop->TVDmomentum);
    } else {
      //UpdateScalars(grid,phys,prop,vert->omega_im,phys->w_SfH_tm1,phys->w_SfH_tm2, phys->w_SfHv_t, phys->w_SfHv_tm1,phys->uold,phys->wc,phys->boundary_w,2,2,phys->Cn_T,
		  //  0,0,NULL,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,0,1,prop->TVDmomentum);
       UpdateScalars(grid,phys,prop,vert->omega_im,phys->w_SfH_tm1,phys->w_SfH_tm2, phys->w_SfHv_t, phys->w_SfHv_tm1,phys->uold,phys->uold2,phys->boundary_w,2,2,phys->Cn_T,
		     0,0,NULL,prop->theta,NULL,NULL,NULL,NULL,0,0,comm,myproc,1,1,prop->TVDmomentum);
    }
    ISendRecvCellData3D(phys->uold,grid,myproc,comm);

    REAL delta, dzdt;
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)  {
	// No interpolation, just averaging of upper and lower layers
	//	delta = 0.5*(phys->uold[i][k-1]-phys->vold[i][k-1]+
	//		     phys->uold[i][k]-phys->vold[i][k]);
	// Linear interpolation to face where w is stored

    //old
  	// delta = ((phys->uold[i][k-1]-phys->vold[i][k-1])*grid->dzz[i][k]+
		//  (phys->uold[i][k]-phys->vold[i][k])*grid->dzz[i][k-1])/
	  // (grid->dzz[i][k-1]+grid->dzz[i][k]);
      
      //QUICK where we can
      if(k>1 && k<(grid->Nk[i]-1)){
        if(.5*(phys->uold[i][k]+phys->uold[i][k-1])>0){
          delta = 0.5*((phys->uold[i][k] - phys->vold[i][k])+(phys->uold[i][k-1] - phys->vold[i][k-1])) - (1/8)*((phys->uold[i][k-1] - phys->vold[i][k-1]) - 2*(phys->uold[i][k] - phys->vold[i][k])+ (phys->uold[i][k+1] - phys->vold[i][k+1]));
        }
        else if(.5*(phys->uold[i][k]+phys->uold[i][k-1])<0){
          delta = 0.5*((phys->uold[i][k] - phys->vold[i][k])+(phys->uold[i][k-1] - phys->vold[i][k-1])) - (1/8)*((phys->uold[i][k] - phys->vold[i][k]) - 2*(phys->uold[i][k-1] - phys->vold[i][k-1])+ (phys->uold[i][k-2] - phys->vold[i][k-2]));
      } else 
	delta = ((phys->uold[i][k-1]-phys->vold[i][k-1])*grid->dzz[i][k]+
		 (phys->uold[i][k]-phys->vold[i][k])*grid->dzz[i][k-1])/
	  (grid->dzz[i][k-1]+grid->dzz[i][k]);
      } else { //use old delta at limits
        delta = ((phys->uold[i][k-1]-phys->vold[i][k-1])*grid->dzz[i][k]+
		    (phys->uold[i][k]-phys->vold[i][k])*grid->dzz[i][k-1])/
	      (grid->dzz[i][k-1]+grid->dzz[i][k]);
      }
	  
    
	phys->w[i][k]+=delta;

      
      }

      k=grid->ctop[i];
      //      if(grid->dzz[i][k]>2.0*vert->vertdzmin)
      phys->w[i][k]+=(phys->uold[i][k]-phys->vold[i][k]);

      k=grid->Nk[i];
      //dzdt = (1.5*phys->zB[i] - 2.0*phys->zBold[i] + 0.5*phys->zBold2[i])/prop->dt;
      //phys->w[i][k]=dzdt
      //+0*phys->uc[i][k-1]*0.5*(3*vert->dzdx[i][k-1]-vert->dzdx[i][k-2])
      //+0*phys->vc[i][k-1]*0.5*(3*vert->dzdy[i][k-1]-vert->dzdy[i][k-2]);
      dzdt = (1.5*phys->zB[i] - 2.0*phys->zBold[i] + 0.5*phys->zBold2[i])/prop->dt;
      phys->w[i][k]=dzdt
        +phys->uc[i][k-1]*0.5*(3*vert->dzdx[i][k-1]-vert->dzdx[i][k-2])
        +phys->vc[i][k-1]*0.5*(3*vert->dzdy[i][k-1]-vert->dzdy[i][k-2]);
    }
  }  

  /*
  if(Wadv==1){ //still need to make w at bed be whatever is happening at bed !!
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      k=grid->Nk[i];
      REAL dzdt;

      dzdt = (1.5*phys->zB[i] - 2.0*phys->zBold[i] + 0.5*phys->zBold2[i])/prop->dt;
      phys->w[i][k]=dzdt
        +phys->uc[i][k-1]*0.5*(3*vert->dzdx[i][k-1]-vert->dzdx[i][k-2])
        +phys->vc[i][k-1]*0.5*(3*vert->dzdy[i][k-1]-vert->dzdy[i][k-2]);
    }

   ISendRecvWData(phys->w,grid,myproc,comm);
  }
  */

  for(i=0;i<grid->Nc;i++){
    for(k=0;k<grid->Nk[i];k++){
      if(phys->w[i][k]!=phys->w[i][k]){
	printf("w nan for i=%d, k=%d \n", i, k);
      }
    }
  }

  //Free Jacobians
  SunFree(J_eff_new,grid->Nkmax*sizeof(REAL),"Wpredictor");
  SunFree(J_eff_old,grid->Nkmax*sizeof(REAL),"Wpredictor");
}

/*
 * Function: EddyViscosity
 * Usage: EddyViscosity(grid,phys,prop,w,comm,myproc);
 * ---------------------------------------------------
 * This function is used to compute the eddy viscosity, the
 * shear stresses, and the drag coefficients at the upper and lower
 * boundaries.
 *
 */
static void EddyViscosity(gridT *grid, physT *phys, propT *prop, REAL **wnew, MPI_Comm comm, int myproc)
{
  int i, k;

  if(prop->turbmodel>=1){
    //my25(grid,phys,prop,wnew,phys->qT,phys->qT_old,phys->lT,phys->lT_old,phys->Cn_q,phys->Cn_l,phys->nu_tv,phys->kappa_tv,comm,myproc);
    my25(grid,phys,prop,wnew,phys->qT,phys->lT,phys->qT_old,phys->lT_old,phys->Cn_q,phys->Cn_l,phys->nu_tv,phys->kappa_tv,comm,myproc);
    //Smagorinsky(grid,phys,prop,comm,myproc);

  }else if(prop->turbmodel==-1){
    for(i=0;i<grid->Nc;i++) 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
	phys->nu_tv[i][k] = 0.052;
  }
}


/*
 * Function: GSSolve
 * Usage: GSSolve(grid,phys,prop,myproc,numprocs,comm);
 * ----------------------------------------------------
 * Solve the free surface equation with a Gauss-Seidell relaxation.
 * This function is used for debugging only.
 *
 */
// classic GS following nearly directly from Oliver's notes
static void GSSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, nf, ne, n, niters, *N;
  REAL *h, *hold, *D, *hsrc, myresid, resid, residold, tmp, relax, myNsrc, Nsrc, coef;

  h = phys->h;
  hold = phys->hold;
  D = phys->D;
  hsrc = phys->htmp;
  N = grid->normal;

  tmp = prop->grav*pow(prop->theta*prop->dt,2);

  // each processor must have all the boundary information 
  // for the starting values of h (as in Lh)
  ISendRecvCellData2D(h,grid,myproc,comm);

  relax = prop->relax;
  niters = prop->maxiters;
  resid=0;
  myresid=0;

  // debugging code since it's not used below
  //    myNsrc=0;
  //
  //    // for interior calculation cells 
  //    // this seems like debugging code
  //    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
  //      i = grid->cellp[iptr];
  //
  //      myNsrc+=pow(hsrc[i],2);
  //    }
  //    MPI_Reduce(&myNsrc,&(Nsrc),1,MPI_DOUBLE,MPI_SUM,0,comm);
  //    MPI_Bcast(&Nsrc,1,MPI_DOUBLE,0,comm);
  //    Nsrc=sqrt(Nsrc);

  for(n=0;n<niters;n++) {

    // hold = h;
    for(i=0;i<grid->Nc;i++) {
      hold[i] = h[i];
    }

    // for all the computational cells (since the boundary cells are 
    // already set in celldist[1] to celldist[...]
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      // get teh cell pointer
      i = grid->cellp[iptr];

      // right hand side (h_src)
      h[i] = hsrc[i];

      coef=1;
      for(nf=0;nf<grid->nfaces[i];nf++) 
        if(grid->neigh[i*grid->maxfaces+nf]!=-1) {
          ne = grid->face[i*grid->maxfaces+nf];

          // coef is the diagonal coefficient term
          coef+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
          h[i]+=relax*tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]*
            phys->h[grid->neigh[i*grid->maxfaces+nf]]/grid->Ac[i];
        }
      // divide by diagonal coefficient term
      h[i]/=coef;
    }

    // now need to compare against the residual term to 
    // determine when GS has converged
    residold=resid;
    myresid=0;
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      hold[i] = hsrc[i];

      coef=1;
      for(nf=0;nf<grid->nfaces[i];nf++) 
        if(grid->neigh[i*grid->maxfaces+nf]!=-1) {
          ne = grid->face[i*grid->maxfaces+nf];
          coef+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]/grid->Ac[i];
          hold[i]+=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne]*
            phys->h[grid->neigh[i*grid->maxfaces+nf]]/grid->Ac[i];
        }
      // compute the residual
      myresid+=pow(hold[i]/coef-h[i],2);
    }
    MPI_Reduce(&myresid,&(resid),1,MPI_DOUBLE,MPI_SUM,0,comm);
    // - is this line necessary?
    MPI_Bcast(&resid,1,MPI_DOUBLE,0,comm);
    resid=sqrt(resid);

    // send all final results out to each processor now
    ISendRecvCellData2D(h,grid,myproc,comm);
    // wait for communication to finish
    MPI_Barrier(comm);

    // if we have met the tolerance criteria
    if(fabs(resid)<prop->epsilon)
      break;
  }
  if(n==niters && myproc==0 && WARNING) 
    printf("Warning... Iteration not converging after %d steps! RES=%e\n",n,resid);

  for(i=0;i<grid->Nc;i++)
    if(h[i]!=h[i]) 
      printf("NaN h[%d] in gssolve!\n",i);
}

/*
 * Function: Continuity
 * Usage: Continuity(w,grid,phys,prop);
 * ------------------------------------
 * Compute the vertical velocity field that satisfies continuity.  Use
 * the upwind flux face heights to ensure consistency with continuity.
 *
 */
void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop)
//static void Continuity(REAL **w, gridT *grid, physT *phys, propT *prop)
{
  int i, k, nf, iptr, ne, nc1, nc2, j, jptr;
  REAL ap, am, dzfnew, theta=prop->theta,fac1,fac2,fac3,Ac,sum;
  
  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=0;k<grid->Nk[i]+1;k++)
      w[i][k] = 0;

    // Continuity is written in terms of flux-face heights at time level n
    // Even though w is updated to grid->ctop[i], it only uses dzf which is
    // 0 for new cells (i.e. if grid->ctop[i]<grid->ctopold[i]).
    // no flux into bottom cells
    w[i][grid->Nk[i]] = 0;
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      // w_old is the old value for w (basically at timestep n) where
      // w in this context is for n+1
      if(prop->subgrid)
        w[i][k] = (subgrid->Acveff[i][k+1]*w[i][k+1]- fac2/fac1*(subgrid->Acveffold[i][k]*phys->w_old[i][k]-subgrid->Acveffold[i][k+1]*phys->w_old[i][k+1])-\
          fac3/fac1*(phys->w_old2[i][k]*subgrid->Acveffold2[i][k]-subgrid->Acveffold2[i][k+1]*phys->w_old2[i][k+1]))/subgrid->Acveff[i][k];
      else
        w[i][k] = w[i][k+1]- fac2/fac1*(phys->w_old[i][k]-phys->w_old[i][k+1])-fac3/fac1*(phys->w_old2[i][k]-phys->w_old2[i][k+1]);

      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        // subgrid change
        if(prop->subgrid)
          Ac=subgrid->Acveff[i][k];
        else 
          Ac=grid->Ac[i];

        // set the flux vertical area is explicit
        if(k<grid->Nke[ne])
          w[i][k]-=(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k])*
            grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac/fac1*grid->dzf[ne][k];
      }
    }
  }

  // calculate w_im for update scalar
  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      if(!prop->subgrid)
        phys->w_im[i][k]=fac2*phys->w_old[i][k]+fac3*phys->w_old2[i][k]+fac1*w[i][k];
      else
        phys->w_im[i][k]=(fac2*phys->w_old[i][k]*subgrid->Acveffold[i][k]+
                fac3*phys->w_old2[i][k]*subgrid->Acveffold2[i][k]+
                fac1*w[i][k]*subgrid->Acveff[i][k])/subgrid->Acveff[i][k];
    phys->w_im[i][grid->Nk[i]]=0;
  }
}

/*
 * Function: ComputeConservatives
 * Usage: ComputeConservatives(grid,phys,prop,myproc,numprocs,comm);
 * -----------------------------------------------------------------
 * Compute the total mass, volume, and potential energy within the entire
 * domain and return a warning if the mass and volume are not conserved to within
 * the tolerance CONSERVED specified in suntans.h 
 *
 */
void ComputeConservatives(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs,
			  MPI_Comm comm)
{
  int i, iptr, k;
  REAL mass, volume, volh, height, Ep;

  if(myproc==0) phys->mass=0;
  if(myproc==0) phys->volume=0;
  if(myproc==0) phys->Ep=0;

  // volh is the horizontal integral of h+d, whereas vol is the
  // 3-d integral of dzz
  mass=0;
  volume=0;
  volh=0;
  Ep=0;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    height = 0;
    volh+=grid->Ac[i]*(grid->dv[i]+phys->h[i]);
    Ep+=0.5*prop->grav*grid->Ac[i]*(phys->h[i]+grid->dv[i])*(phys->h[i]-grid->dv[i]);
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      height += grid->dzz[i][k];
      volume+=grid->Ac[i]*grid->dzz[i][k];
      mass+=phys->s[i][k]*grid->Ac[i]*grid->dzz[i][k];
    }
  }

  // Comment out the volh reduce if that integral is desired.  The
  // volume integral is used since the volh integral is useful
  // only for debugging.
  MPI_Reduce(&mass,&(phys->mass),1,MPI_DOUBLE,MPI_SUM,0,comm);
  //MPI_Reduce(&volh,&(phys->volume),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&volume,&(phys->volume),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&Ep,&(phys->Ep),1,MPI_DOUBLE,MPI_SUM,0,comm);

  // Compare the quantities to the original values at the beginning of the
  // computation.  If prop->n==0 (beginning of simulation), then store the
  // starting values for comparison.
  if(myproc==0) {
    if(prop->n==0) {
      phys->volume0 = phys->volume;
      phys->mass0 = phys->mass;
      phys->Ep0 = phys->Ep;
    } else {
      if(fabs((phys->volume-phys->volume0)/phys->volume0)>CONSERVED && prop->volcheck)
        printf("Warning! Not volume conservative at step %d! V(0)=%e, V(t)=%e\n",prop->n,
            phys->volume0,phys->volume);
      if(fabs((phys->mass-phys->mass0)/phys->volume0)>CONSERVED && prop->masscheck) 
        printf("Warning! Not mass conservative at step %d! M(0)=%e, M(t)=%e\n", prop->n,
            phys->mass0,phys->mass);
    }
  }
}

/*
 * Function: ComputeUCPerot
 * Usage: ComputeUCPerot(u,uc,vc,grid);
 * -------------------------------------------
 * Compute the cell-centered components of the velocity vector and place them
 * into uc and vc.  This function estimates the velocity vector with
 *
 * u = 1/Area * Sum_{faces} u_{face} normal_{face} df_{face}*d_{ef,face}
 *
 */
void ComputeUCPerot(REAL **u, REAL **uc, REAL **vc, REAL *h, int kinterp, int subgridmodel, gridT *grid) {

  int k, n, ne, nf, iptr;
  REAL sum;

  // for each computational cell (non-stage defined)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    // get cell pointer transfering from boundary coordinates 
    // to grid coordinates
    n=grid->cellp[iptr];

    // initialize over all depths
    for(k=0;k<grid->Nk[n];k++) {
      uc[n][k]=0;
      vc[n][k]=0;
    }
    // over all interior cells
    for(k=grid->ctop[n];k<grid->Nk[n];k++) {
      // over each face
      for(nf=0;nf<grid->nfaces[n];nf++) {
        ne = grid->face[n*grid->maxfaces+nf];
        if(!(grid->smoothbot) || k<grid->Nke[ne]){
          uc[n][k]+=u[ne][k]*grid->n1[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];//*grid->dzf[ne][k];
          vc[n][k]+=u[ne][k]*grid->n2[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];//*grid->dzf[ne][k];
        }else{ 
          uc[n][k]+=u[ne][grid->Nke[ne]-1]*grid->n1[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];//*grid->dzf[ne][grid->Nke[ne]-1];
          vc[n][k]+=u[ne][grid->Nke[ne]-1]*grid->n2[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];//*grid->dzf[ne][grid->Nke[ne]-1];
        }
      }

      // In case of divide by zero (shouldn't happen)
      if (grid->dzz[n][k] > DRYCELLHEIGHT) {
	uc[n][k]/=grid->Ac[n];//(grid->Ac[n]*grid->dzz[n][k]);
	vc[n][k]/=grid->Ac[n];//(grid->Ac[n]*grid->dzz[n][k]);
      } else {
          uc[n][k] = 0;
          vc[n][k] = 0;
      }
    }

    //top cell only - don't account for depth
    /*
    k=grid->ctop[n];
    // over each face
    for(nf=0;nf<grid->nfaces[n];nf++) {
      ne = grid->face[n*grid->maxfaces+nf];
      if(!(grid->smoothbot) || k<grid->Nke[ne]){
        uc[n][k]+=u[ne][k]*grid->n1[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
        vc[n][k]+=u[ne][k]*grid->n2[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
      }else{ 
        uc[n][k]+=u[ne][grid->Nke[ne]-1]*grid->n1[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
        vc[n][k]+=u[ne][grid->Nke[ne]-1]*grid->n2[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];       
      }
    }
    uc[n][k]/=grid->Ac[n];
    vc[n][k]/=grid->Ac[n];
    */
  }
}


/*static void ComputeUCPerot(REAL **u, REAL **uc, REAL **vc, REAL *h, int kinterp, int subgridmodel, gridT *grid) {

  int i,k, n, ne, nf, iptr,nc1,nc2,dry=1;
  REAL alpha,d,V;
  REAL sum;

  // for each computational cell (non-stage defined)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    // get cell pointer transfering from boundary coordinates 
    // to grid coordinates
    n=grid->cellp[iptr];
    d=h[n]+grid->dv[n];

    // initialize over all depths
    for(k=0;k<grid->Nk[n];k++) {
      uc[n][k]=0;
      vc[n][k]=0;
    }
    
    // over the entire depth (cell depth)
    for(k=grid->ctop[n];k<grid->Nk[n];k++) {
      V=0;
      dry=1;
      for(nf=0;nf<grid->nfaces[n];nf++) {
        ne = grid->face[n*grid->maxfaces+nf];
        if(grid->dzf[ne][k]>0)
          dry=0;
        V+=grid->dzf[ne][k]*grid->df[ne]*grid->def[n*grid->maxfaces+nf];
      }        
      V/=2;
      // over each face
      for(nf=0;nf<grid->nfaces[n];nf++) {
        ne = grid->face[n*grid->maxfaces+nf];
        nc1=grid->grad[2*ne];
        nc2=grid->grad[2*ne+1];
        if(nc1==-1)
          nc1=nc2;
        if(nc2==-1)
          nc2=nc1;
      
        if(subgridmodel)
          if(V/subgrid->Acceff[n][k]/grid->dzz[n][k]<=1)
            alpha=grid->dzf[ne][k]/grid->dzz[n][k]*grid->Ac[n]/subgrid->Acceff[n][k];
          else
            alpha=grid->dzf[ne][k]*grid->Ac[n]/V;
        else{
          if(V<=grid->dzz[n][k]*grid->Ac[n])
            alpha=grid->dzf[ne][k]/grid->dzz[n][k];
          else
            alpha=grid->dzf[ne][k]*grid->Ac[n]/V;
        }

        if(dry)
          alpha=1;

        //most stable with different dzf
        // not best for varying dzf
        alpha=1;
        if(!(grid->smoothbot)|| k<grid->Nke[ne]){
          uc[n][k]+=alpha*u[ne][k]*grid->n1[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
          vc[n][k]+=alpha*u[ne][k]*grid->n2[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
        }
        else{	
          uc[n][k]+=alpha*u[ne][grid->Nke[ne]-1]*grid->n1[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
          vc[n][k]+=alpha*u[ne][grid->Nke[ne]-1]*grid->n2[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
        }

        if(grid->dzz[n][k]<=DRYCELLHEIGHT)
        {
          uc[n][k]=0;
          vc[n][k]=0;
        }  
      }

      uc[n][k]/=grid->Ac[n];
      vc[n][k]/=grid->Ac[n];
    }
  }
}*/

/*
 * Function: ReadProperties
 * Usage: ReadProperties(prop,grid,myproc);
 * -----------------------------------
 * This function reads in the properties specified in the suntans.dat
 * data file.  Note that if an entry does not exist, a default can be used.
 *
 */
void ReadProperties(propT **prop, gridT *grid, int myproc)
{ 
  // allocate memory
  *prop = (propT *)SunMalloc(sizeof(propT),"ReadProperties");

  // set values from suntans.dat file (DATAFILE)
  (*prop)->thetaramptime = MPI_GetValue(DATAFILE,"thetaramptime","ReadProperties",myproc);
  (*prop)->theta = MPI_GetValue(DATAFILE,"theta","ReadProperties",myproc);
  (*prop)->thetaS = MPI_GetValue(DATAFILE,"thetaS","ReadProperties",myproc);
  (*prop)->thetaB = MPI_GetValue(DATAFILE,"thetaB","ReadProperties",myproc);
  (*prop)->beta = MPI_GetValue(DATAFILE,"beta","ReadProperties",myproc);
  (*prop)->kappa_s = MPI_GetValue(DATAFILE,"kappa_s","ReadProperties",myproc);
  (*prop)->kappa_sH = MPI_GetValue(DATAFILE,"kappa_sH","ReadProperties",myproc);
  (*prop)->gamma = MPI_GetValue(DATAFILE,"gamma","ReadProperties",myproc);
  (*prop)->kappa_T = MPI_GetValue(DATAFILE,"kappa_T","ReadProperties",myproc);
  (*prop)->kappa_TH = MPI_GetValue(DATAFILE,"kappa_TH","ReadProperties",myproc);
  (*prop)->nu = MPI_GetValue(DATAFILE,"nu","ReadProperties",myproc);
  (*prop)->nu_H = MPI_GetValue(DATAFILE,"nu_H","ReadProperties",myproc);
  (*prop)->tau_T = MPI_GetValue(DATAFILE,"tau_T","ReadProperties",myproc);
  (*prop)->z0T = MPI_GetValue(DATAFILE,"z0T","ReadProperties",myproc);
  (*prop)->z0B = MPI_GetValue(DATAFILE,"z0B","ReadProperties",myproc);
  (*prop)->BUFFERHEIGHT = MPI_GetValue(DATAFILE,"BUFFERHEIGHT","ReadProperties",myproc);  
  (*prop)->Intz0B= MPI_GetValue(DATAFILE,"Intz0B","ReadProperties",myproc);
  (*prop)->Intz0T= MPI_GetValue(DATAFILE,"Intz0T","ReadProperties",myproc); 
  (*prop)->cwave = MPI_GetValue(DATAFILE,"cwave","ReadProperties",myproc); 
  (*prop)->wavelength = MPI_GetValue(DATAFILE,"wavelength","ReadProperties",myproc);

   // user_def_var for output 
   (*prop)->output_user_var= MPI_GetValue(DATAFILE,"outputuservar","ReadProperties",myproc); 

  if((*prop)->Intz0B==1)
    MPI_GetFile((*prop)->INPUTZ0BFILE,DATAFILE,"inputz0Bfile","ReadFileNames",myproc);
  if((*prop)->Intz0T==1)
    MPI_GetFile((*prop)->INPUTZ0TFILE,DATAFILE,"inputz0Tfile","ReadFileNames",myproc);

  (*prop)->CdT = MPI_GetValue(DATAFILE,"CdT","ReadProperties",myproc);
  (*prop)->CdB = MPI_GetValue(DATAFILE,"CdB","ReadProperties",myproc);
  (*prop)->CdW = MPI_GetValue(DATAFILE,"CdW","ReadProperties",myproc);
  (*prop)->grav= MPI_GetValue(DATAFILE,"grav","ReadProperties",myproc);
  (*prop)->turbmodel = (int)MPI_GetValue(DATAFILE,"turbmodel","ReadProperties",myproc);
  (*prop)->dt = MPI_GetValue(DATAFILE,"dt","ReadProperties",myproc);
  (*prop)->Cmax = MPI_GetValue(DATAFILE,"Cmax","ReadProperties",myproc);
  (*prop)->nsteps = (int)MPI_GetValue(DATAFILE,"nsteps","ReadProperties",myproc);
  (*prop)->ntout = (int)MPI_GetValue(DATAFILE,"ntout","ReadProperties",myproc);
  (*prop)->ntoutStore = (int)MPI_GetValue(DATAFILE,"ntoutStore","ReadProperties",myproc);

  if((*prop)->ntoutStore==0)
    (*prop)->ntoutStore=(*prop)->nsteps;

  (*prop)->ntprog = (int)MPI_GetValue(DATAFILE,"ntprog","ReadProperties",myproc);
  (*prop)->ntconserve = (int)MPI_GetValue(DATAFILE,"ntconserve","ReadProperties",myproc);
  (*prop)->nonhydrostatic = (int)MPI_GetValue(DATAFILE,"nonhydrostatic","ReadProperties",myproc);
  (*prop)->cgsolver = (int)MPI_GetValue(DATAFILE,"cgsolver","ReadProperties",myproc);
  (*prop)->include_slope_terms = (int)MPI_GetValue(DATAFILE,"include_slope_terms","ReadProperties",myproc);
  (*prop)->movingBed = (int)MPI_GetValue(DATAFILE,"movingBed","ReadProperties",myproc);
  (*prop)->pressureMethod = (int)MPI_GetValue(DATAFILE,"pressureMethod","ReadProperties",myproc);    
  (*prop)->advectSalinity = (int)MPI_GetValue(DATAFILE,"advectSalinity","ReadProperties",myproc);
  (*prop)->advectTemperature = (int)MPI_GetValue(DATAFILE,"advectTemperature","ReadProperties",myproc);    
  (*prop)->maxiters = (int)MPI_GetValue(DATAFILE,"maxiters","ReadProperties",myproc);
  (*prop)->qmaxiters = (int)MPI_GetValue(DATAFILE,"qmaxiters","ReadProperties",myproc);
  (*prop)->qprecond = (int)MPI_GetValue(DATAFILE,"qprecond","ReadProperties",myproc);
  (*prop)->epsilon = MPI_GetValue(DATAFILE,"epsilon","ReadProperties",myproc);
  (*prop)->qepsilon = MPI_GetValue(DATAFILE,"qepsilon","ReadProperties",myproc);
  (*prop)->resnorm = MPI_GetValue(DATAFILE,"resnorm","ReadProperties",myproc);
  (*prop)->relax = MPI_GetValue(DATAFILE,"relax","ReadProperties",myproc);
  (*prop)->amp = MPI_GetValue(DATAFILE,"amp","ReadProperties",myproc);
  (*prop)->omega = MPI_GetValue(DATAFILE,"omega","ReadProperties",myproc);
  (*prop)->timescale = MPI_GetValue(DATAFILE,"timescale","ReadProperties",myproc);
  (*prop)->flux = MPI_GetValue(DATAFILE,"flux","ReadProperties",myproc);
  (*prop)->volcheck = MPI_GetValue(DATAFILE,"volcheck","ReadProperties",myproc);
  (*prop)->masscheck = MPI_GetValue(DATAFILE,"masscheck","ReadProperties",myproc);
  (*prop)->nonlinear = MPI_GetValue(DATAFILE,"nonlinear","ReadProperties",myproc);
  (*prop)->wetdry = MPI_GetValue(DATAFILE,"wetdry","ReadProperties",myproc);
  (*prop)->Coriolis_f = MPI_GetValue(DATAFILE,"Coriolis_f","ReadProperties",myproc);
  (*prop)->sponge_distance = MPI_GetValue(DATAFILE,"sponge_distance","ReadProperties",myproc);
  (*prop)->sponge_decay = MPI_GetValue(DATAFILE,"sponge_decay","ReadProperties",myproc);
  (*prop)->readSalinity = MPI_GetValue(DATAFILE,"readSalinity","ReadProperties",myproc);
  (*prop)->readOldVelocity = MPI_GetValue(DATAFILE,"readOldVelocity","ReadProperties",myproc);
  (*prop)->readTemperature = MPI_GetValue(DATAFILE,"readTemperature","ReadProperties",myproc);
  (*prop)->TVDsalt = MPI_GetValue(DATAFILE,"TVDsalt","ReadProperties",myproc);
  (*prop)->TVDtemp = MPI_GetValue(DATAFILE,"TVDtemp","ReadProperties",myproc);
  (*prop)->TVDturb = MPI_GetValue(DATAFILE,"TVDturb","ReadProperties",myproc);
  (*prop)->stairstep = MPI_GetValue(DATAFILE,"stairstep","ReadProperties",myproc);
  (*prop)->TVDmomentum = MPI_GetValue(DATAFILE,"TVDmomentum","ReadProperties",myproc); 
  (*prop)->conserveMomentum = MPI_GetValue(DATAFILE,"conserveMomentum","ReadProperties",myproc); 
  (*prop)->thetaM = MPI_GetValue(DATAFILE,"thetaM","ReadProperties",myproc); 
  (*prop)->newcells = MPI_GetValue(DATAFILE,"newcells","ReadProperties",myproc); 
  (*prop)->mergeArrays = MPI_GetValue(DATAFILE,"mergeArrays","ReadProperties",myproc); 
  (*prop)->computeSediments = MPI_GetValue(DATAFILE,"computeSediments","ReadProperties",myproc); 
  (*prop)->subgrid = MPI_GetValue(DATAFILE,"subgrid","ReadProperties",myproc); 
  (*prop)->marshmodel = MPI_GetValue(DATAFILE,"marshmodel","ReadProperties",myproc);
  (*prop)->wavemodel = MPI_GetValue(DATAFILE,"wavemodel","ReadProperties",myproc);
  (*prop)->culvertmodel = MPI_GetValue(DATAFILE,"culvertmodel","ReadProperties",myproc);
  (*prop)->vertcoord = MPI_GetValue(DATAFILE,"vertcoord","ReadProperties",myproc);
  (*prop)->ex = MPI_GetValue(DATAFILE,"ex","ReadProperties",myproc); //AB3
  (*prop)->im = MPI_GetValue(DATAFILE,"im","ReadProperties",myproc); //implicit method for momentum eqn.
  (*prop)->botLayers = (int)MPI_GetValue(DATAFILE,"botLayers","ReadProperties",myproc); 
  (*prop)->midLayers = (int)MPI_GetValue(DATAFILE,"midLayers","ReadProperties",myproc); 
  (*prop)->cuttoffDepth = MPI_GetValue(DATAFILE,"cuttoffDepth","ReadProperties",myproc); 
  
  // setup the factors for implicit and explicit schemes
  if((*prop)->ex==1) // AX2
  {
    (*prop)->exfac1=7.0/4.0;
    (*prop)->exfac2=-1.0;
    (*prop)->exfac3=1.0/4.0;
  } else if((*prop)->ex==2) { // AB2
    (*prop)->exfac1=1.5;
    (*prop)->exfac2=-0.5;
    (*prop)->exfac3=0;
  } else if((*prop)->ex==3) { //AB3
    (*prop)->exfac1=23.0/12.0;
    (*prop)->exfac2=-4.0/3.0;
    (*prop)->exfac3=5.0/12.0;
  } else { // Forward Euler
    (*prop)->exfac1=1;
    (*prop)->exfac2=0;
    (*prop)->exfac3=0;
  }

  if((*prop)->im==0) // theta
  {
    (*prop)->imfac1=(*prop)->theta;
    (*prop)->imfac2=1.0-(*prop)->theta;
    (*prop)->imfac3=0;
  } else if((*prop)->im==1) { // AM2
    (*prop)->imfac1=0.75;
    (*prop)->imfac2=0;
    (*prop)->imfac3=0.25;
  } else { //AI2
    (*prop)->imfac1=5.0/4.0;
    (*prop)->imfac2=-1.0;
    (*prop)->imfac3=3.0/4.0;
  }

  // When wetting and drying is desired:
  // -Do nonconservative momentum advection (conserveMomentum=0)
  // -Use backward Euler for vertical advection of horizontal momentum (thetaM=1)
  // -Update new cells (newcells=1)
  if((*prop)->wetdry) {
    (*prop)->conserveMomentum = 0;
    (*prop)->thetaM = 1;//Fully implicit
    //(*prop)->thetaM = 0.5;
    if((*prop)->vertcoord==1 || (*prop)->vertcoord==5)
      (*prop)->newcells = 1;
  }
  
  (*prop)->calcage = MPI_GetValue(DATAFILE,"calcage","ReadProperties",myproc);
  (*prop)->agemethod = MPI_GetValue(DATAFILE,"agemethod","ReadProperties",myproc);
  (*prop)->calcaverage = MPI_GetValue(DATAFILE,"calcaverage","ReadProperties",myproc);
  if ((*prop)->calcaverage)
      (*prop)->ntaverage = (int)MPI_GetValue(DATAFILE,"ntaverage","ReadProperties",myproc);
  (*prop)->latitude = MPI_GetValue(DATAFILE,"latitude","ReadProperties",myproc);
  (*prop)->gmtoffset = MPI_GetValue(DATAFILE,"gmtoffset","ReadProperties",myproc);
  (*prop)->metmodel = (int)MPI_GetValue(DATAFILE,"metmodel","ReadProperties",myproc);
  (*prop)->varmodel = (int)MPI_GetValue(DATAFILE,"varmodel","ReadProperties",myproc);
  (*prop)->nugget = MPI_GetValue(DATAFILE,"nugget","ReadProperties",myproc);
  (*prop)->sill = MPI_GetValue(DATAFILE,"sill","ReadProperties",myproc);
  (*prop)->range = MPI_GetValue(DATAFILE,"range","ReadProperties",myproc);
  (*prop)->outputNetcdf = (int)MPI_GetValue(DATAFILE,"outputNetcdf","ReadProperties",myproc);
  (*prop)->netcdfBdy = (int)MPI_GetValue(DATAFILE,"netcdfBdy","ReadProperties",myproc);
  (*prop)->readinitialnc = (int)MPI_GetValue(DATAFILE,"readinitialnc","ReadProperties",myproc);
  (*prop)->Lsw = MPI_GetValue(DATAFILE,"Lsw","ReadProperties",myproc);
  (*prop)->Cda = MPI_GetValue(DATAFILE,"Cda","ReadProperties",myproc);
  (*prop)->Ce = MPI_GetValue(DATAFILE,"Ce","ReadProperties",myproc);
  (*prop)->Ch = MPI_GetValue(DATAFILE,"Ch","ReadProperties",myproc);
  if((*prop)->outputNetcdf > 0 || (*prop)->netcdfBdy > 0 || (*prop)->readinitialnc > 0){ 
      MPI_GetString((*prop)->starttime,DATAFILE,"starttime","ReadProperties",myproc);
      MPI_GetString((*prop)->basetime,DATAFILE,"basetime","ReadProperties",myproc);

      (*prop)->nstepsperncfile=(int)MPI_GetValue(DATAFILE,"nstepsperncfile","ReadProperties",myproc);
      (*prop)->ncfilectr=(int)MPI_GetValue(DATAFILE,"ncfilectr","ReadProperties",myproc);
  }
  
  if((*prop)->nonlinear==2) {
    (*prop)->laxWendroff = MPI_GetValue(DATAFILE,"laxWendroff","ReadProperties",myproc);
    if((*prop)->laxWendroff!=0)
      (*prop)->laxWendroff_Vertical = MPI_GetValue(DATAFILE,"laxWendroff_Vertical","ReadProperties",myproc);
    else
      (*prop)->laxWendroff_Vertical = 0;
  } else {
    (*prop)->laxWendroff = 0;
    (*prop)->laxWendroff_Vertical = 0;
  }

  (*prop)->hprecond = MPI_GetValue(DATAFILE,"hprecond","ReadProperties",myproc);

  // addition for interpolation methods
  switch((int)MPI_GetValue(DATAFILE,"interp","ReadProperties",myproc)) {
    case 0: //Perot
      (*prop)->interp = PEROT;
      break;
    case 1: //Quad
      (*prop)->interp = QUAD;
      break;
    case 2: //Least-squares
      (*prop)->interp = LSQ;
      break;
 
    default:
      printf("ERROR: Specification of interpolation type is incorrect!\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
      break;
  }
  (*prop)->kinterp=(int)MPI_GetValue(DATAFILE,"kinterp","ReadProperties",myproc);

  if((int)MPI_GetValue(DATAFILE,"kinterp","ReadProperties",myproc))
  {  
    (*prop)->interp=PEROT;
    //if(myproc==0)
      //printf("kinterp is used, so interp is set as perot automatically\n");
  }
  
  if((*prop)->interp==QUAD && grid->maxfaces>DEFAULT_NFACES) {
    //printf("Warning in ReadProperties...interp set to PEROT for use with quad or hybrid grid.\n");
    (*prop)->interp=PEROT;
  }

  // additional data for pretty plot methods
  (*prop)->prettyplot = MPI_GetValue(DATAFILE,"prettyplot","ReadProperties",myproc);
  if((*prop)->prettyplot!=0 && grid->maxfaces>DEFAULT_NFACES) {
    printf("Warning in ReadProperties...prettyplot set to zero for use with quad or hybrid grid.\n");
    (*prop)->prettyplot=0;
  }

  // addition for linearized free surface where dzz=dz
  (*prop)->linearFS = (int)MPI_GetValue(DATAFILE,"linearFS","ReadProperties",myproc);
}

/*
 * Function: InterpToFace
 * Usage: uface = InterpToFace(j,k,phys->uc,u,grid);
 * -------------------------------------------------
 * Linear interpolation of a Voronoi-centered value to the face, using the equation
 * 
 * uface = 1/Dj*(def1*u2 + def2*u1);
 *
 * Note that def1 and def2 are not the same as grid->def[] unless the
 * triangles have not been corrected.  This affects the Coriolis term as currently 
 * implemented.
 *
 */
REAL InterpToFace(int j, int k, REAL **phi, REAL **u, gridT *grid) {
  int nc1, nc2;
  REAL def1, def2, Dj;
  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  if(nc1==-1)
    nc1=nc2;
  if(nc2==-1)
    nc2=nc1;

  Dj = grid->dg[j];
  Return_def(&def1,&def2,nc1,nc2,j,grid);	   	      
  //  def1=grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
  //  def2=Dj-def1;

  if(def1==0 || def2==0) {
    return UpWind(u[j][k],phi[nc1][k],phi[nc2][k]);
  }
  else {
    return (phi[nc1][k]*def2+phi[nc2][k]*def1)/(def1+def2);
  }
}

/*
 * Function: UFaceFlux
 * Usage: UFaceFlux(j,k,phi,phys->u,grid,prop->dt,prop->nonlinear);
 * ---------------------------------------------------------------------
 * Interpolation to obtain the flux of the scalar field phi (of type REAL **) on 
 * face j, k;  method==2: Central-differencing, method==4: Lax-Wendroff.
 *
 */
// note that we should not compute def1 and def2 in this function as they don't change
// these should be computed in grid.c and just looked up for efficiency
static REAL UFaceFlux(int j, int k, REAL **phi, REAL **u, gridT *grid, REAL dt, int method) {
  int nc1, nc2;
  REAL def1, def2, Dj, C=0;
  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  if(nc1==-1) nc1=nc2;
  if(nc2==-1) nc2=nc1;
  Dj = grid->dg[j];
  Return_def(&def1,&def2,nc1,nc2,j,grid);	   	        
  //  def1=grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
  //  def2=Dj-def1;

  if(method==4) C = u[j][k]*dt/Dj;
  if(method==2) C = 0;
  if(def1==0 || def2==0 || method==1) {
    // this happens on a boundary cell
    return UpWind(u[j][k],phi[nc1][k],phi[nc2][k]);
  }
  else {
    // on an interior cell (orthogonal cell) we can just distance interpolation the 
    // value to the face)
    // basically Eqn 51
    return (phi[nc1][k]*def2+phi[nc2][k]*def1)/(def1+def2)
      -C/2*(phi[nc1][k]-phi[nc2][k]);
  }
}

/*
 * Function: SetDensity
 * Usage: SetDensity(grid,phys,prop);
 * ----------------------------------
 * Sets the values of the density in the density array rho and
 * at the boundaries.
 *
 */
void SetDensity(gridT *grid, physT *phys, propT *prop) {
  int i, j, k, jptr, ib;
  REAL z, p;

  for(i=0;i<grid->Nc;i++) {
    z=phys->h[i];
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      z+=0.5*grid->dzz[i][k];
      p=RHO0*prop->grav*z;
      phys->rho[i][k]=StateEquation(prop,phys->s[i][k],phys->T[i][k],p);
      z+=0.5*grid->dzz[i][k];
    }
  }

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j=grid->edgep[jptr];
    ib=grid->grad[2*j];
    z=phys->h[ib];
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      z+=0.5*grid->dzz[ib][k];
      p=RHO0*prop->grav*z;
      phys->boundary_rho[jptr-grid->edgedist[2]][k]=
        StateEquation(prop,phys->boundary_s[jptr-grid->edgedist[2]][k],
            phys->boundary_T[jptr-grid->edgedist[2]][k],p);
      z+=0.5*grid->dzz[ib][k];
    }
  }
}

/*
 * Function: SetFluxHeight
 * Usage: SetFluxHeight(grid,phys,prop);
 * -------------------------------------
 * Set the value of the flux height dzf at time step n for use
 * in continuity and scalar transport.
 *
 */
void SetFluxHeight(gridT *grid, physT *phys, propT *prop, int dzfmeth, MPI_Comm comm, int myproc) {
  int j, jptr, k, nc1, nc2;
  REAL u_im;

  /*    
  for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nkmax;k++) {
      grid->dzf[j][k]=0;
    }
  }
  */

  for(j=0;j<grid->Ne;j++) {                                                                                                                                         
    for(k=0;k<grid->Nkmax;k++) {
      nc1 = grid->grad[2*j]; 
      nc2 = grid->grad[2*j+1]; 
      //printf("before set dzf: dzz[nc1][k] is %f, dzz[nc2][k] is %f, dz[k] is %f, and dzf[j][k] is %f \n", grid->dzz[nc1][k], grid->dzz[nc2][k], grid->dz[k], grid->dzf[j][k]);
      //grid->dzf[j][k]=0;                                                                                                                                            
    }                                                                                                                                                               
  }   
  
  // Set flux height at boundaries as the height of interior point
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) { //changed frist from 2
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];
    
    for(k=0;k<grid->Nkmax;k++)
      grid->dzf[j][k]=grid->dzz[nc1][k];
  }

  int num_layers_bot = prop->botLayers; //40 (45) (4) (30). (40 for sandwaves most recently!)

  if(prop->vertcoord==2 || prop->vertcoord==4) {
    TvdFluxHeight(grid, phys, prop, dzfmeth, comm, myproc);  

    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      
      for(k=0;k<grid->Nkmax-num_layers_bot;k++){
	u_im = prop->imfac2*phys->u_old[j][k]+prop->imfac1*phys->u[j][k]+prop->imfac3*phys->u_old2[j][k];
		
	if(u_im>=0)
	  grid->dzf[j][k]=phys->SfHp[j][k];
	else
	  grid->dzf[j][k]=phys->SfHm[j][k];
      } for(k=grid->Nkmax-num_layers_bot;k<grid->Nkmax;k++){
          grid->dzf[j][k]=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
      }
    }
  } else {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];

      for(k=0;k<grid->Nkmax;k++) {
	  grid->dzf[j][k]=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
      }
    }
  }

  // for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
  //     j = grid->edgep[jptr];
      
  //     nc1 = grid->grad[2*j];
  //     nc2 = grid->grad[2*j+1];

  //     for(k=0;k<grid->Nkmax;k++) {
	//   grid->dzf[j][k]=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);
  //     }
  //   }

  /*
  for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nkmax;k++) {
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];

      if(nc1==-1)
	nc1=nc2;
      else if(nc2==-1)
	nc2=nc1;

      printf("for nc1 =%d, nc2=%d, k=%d: after set dzf: dzz[nc1][k] is %f, dzz[nc2][k] is %f, dz[k] is %f, and dzf[j][k] is %f \n", nc1, nc2, k, grid->dzz[nc1][k], grid->dzz[nc2][k], grid->dz[k], grid->dzf[j][k]);
      //grid->dzf[j][k]=0;                                                                                                                                          
                                                                                                                                                                    
    }
  }
  */
  
  /*for(j=0;j<grid->Ne;j++) {
    for(k=0;k<grid->Nkmax;k++) {
      if(grid->dzf[j][k]==0 && grid->mark[j]!=1)
	printf("DZF = 0 at k=%d, marker=%d\n",k,grid->mark[j]);
    }
  }
  */
  
}    

/*
 * Function: ComputeUC
 * Usage: ComputeUC(u,uc,vc,grid);
 * -------------------------------------------
 * Compute the cell-centered components of the velocity vector and place them
 * into uc and vc.  
 *
 */
//inline static void ComputeUC(physT *phys, gridT *grid, int myproc) {
void ComputeUC(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc, interpolation interp,int kinterp, int subgridmodel) {

  switch(interp) {
    case QUAD:
      // using Wang et al 2011 methods
      ComputeUCRT(ui, vi, phys,grid, myproc);
      break;
    case PEROT:
      ComputeUCPerot(phys->u,ui,vi,phys->h,kinterp,subgridmodel,grid);
      break;
    case LSQ:
      ComputeUCLSQ(phys->u,ui,vi,grid,phys);
      break;
    default:
      break;
  }

}

/*
 * Function: ComputeUCRT
 * Usage: ComputeUCRT(u,uc,vc,grid);
 * -------------------------------------------
 * Compute the cell-centered components of the velocity vector and place them
 * into uc and vc.  This function estimates the velocity vector with
 * methods outlined in Wang et al, 2011
 *
 */
static void ComputeUCRT(REAL **ui, REAL **vi, physT *phys, gridT *grid, int myproc) {

  int k, n, ne, nf, iptr;
  REAL sum;

  // first we need to reconstruct the nodal velocities using the RT0 basis functions
  //  if(myproc==0) printf("ComputeNodalVelocity\n");
  ComputeNodalVelocity(phys, grid, nRT2, myproc);

  // now we can get the tangential velocities from these results
  //  if(myproc==0) printf("ComputeTangentialVelocity\n");
  ComputeTangentialVelocity(phys, grid, nRT2, tRT2, myproc);

  // for each computational cell (non-stage defined)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    // get cell pointer transfering from boundary coordinates 
    // to grid coordinates
    n=grid->cellp[iptr];

    // initialize over all depths
    for(k=0;k<grid->Nk[n];k++) {
      phys->uc[n][k]=0;
      phys->vc[n][k]=0;
    }

    // over the entire depth (cell depth)
    // possible mistake here- make sure to the get the right value for the start
    // of the loop (etop?)
    //    if(myproc==0) printf("ComputeQuadraticInterp\n");
    if(grid->nfaces[n]==3){
      for(k=grid->ctop[n];k<grid->Nk[n];k++) {
        // now we can compute the quadratic interpolated velocity from these results
        ComputeQuadraticInterp(grid->xv[n], grid->yv[n], n, k, ui, 
            vi, phys, grid, nRT2, tRT2, myproc);
      }
    } else {
       // over the entire depth (cell depth)
       for(k=grid->ctop[n];k<grid->Nk[n];k++) {
        // over each face
        for(nf=0;nf<grid->nfaces[n];nf++) {
          ne = grid->face[n*grid->maxfaces+nf];
          if(!(grid->smoothbot) || k<grid->Nke[ne]){
            phys->uc[n][k]+=phys->u[ne][k]*grid->n1[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
            phys->vc[n][k]+=phys->u[ne][k]*grid->n2[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
          }
          else{	
            phys->uc[n][k]+=phys->u[ne][grid->Nke[ne]-1]*grid->n1[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
            phys->vc[n][k]+=phys->u[ne][grid->Nke[ne]-1]*grid->n2[ne]*grid->def[n*grid->maxfaces+nf]*grid->df[ne];
          }
        }

        phys->uc[n][k]/=grid->Ac[n];
        phys->vc[n][k]/=grid->Ac[n];
      }
    }
  } 
}

/*
 * Function: ComputeUCLSQ
 * Usage: ComputeUCRT(u,uc,vc,grid);
 * -------------------------------------------
 * Compute the cell-centered components of the velocity vector and place them
 * into uc and vc.  This function estimates the velocity vector with the least square method
 *
 */
static void ComputeUCLSQ(REAL **u, REAL **uc, REAL **vc, gridT *grid, physT *phys){
  int k, n, ne, nf, iptr;
  REAL sum;
  int ii,jj,kk;
  REAL **A = phys->A;
  REAL **AT = phys->AT;
  REAL **Apr = phys->Apr;
  REAL *bpr = phys->bpr;

  // for each computational cell (non-stage defined)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    // get cell pointer transfering from boundary coordinates 
    // to grid coordinates
    n=grid->cellp[iptr];

    // initialize over all depths
    for(k=0;k<grid->Nk[n];k++) {
      uc[n][k]=0;
      vc[n][k]=0;
    }
    for(nf=0;nf<grid->nfaces[n];nf++) {

        ne = grid->face[n*grid->maxfaces+nf];
        
        // Construct the normal array - A
        A[nf][0] = grid->n1[ne];
        A[nf][1] = grid->n2[ne];

        // Construct A transpose
        AT[0][nf] = A[nf][0];
        AT[1][nf] = A[nf][1];
    }

    for(k=grid->ctop[n];k<grid->Nk[n];k++) {
        
        // Multiply A' = A^T*A 
        // This can't be moved outside of this for loop because A' is
        // modified by linsolve()
        for(ii=0;ii<2;ii++){
            for(jj=0;jj<2;jj++){
                sum=0;
                for(nf=0;nf<grid->nfaces[n];nf++) {
                    sum += AT[ii][nf]*A[nf][jj];
                }
                Apr[ii][jj]=sum;
            }
        }
        // Compute b' = A^T * b
        for(ii=0;ii<2;ii++){
            sum=0;
            for(nf=0;nf<grid->nfaces[n];nf++) {
                ne = grid->face[n*grid->maxfaces+nf];
                if(!(grid->smoothbot) || k<grid->Nke[ne]){
                    sum += AT[ii][nf]*u[ne][k];
                    //sum += AT[ii][nf]*u[ne][k]*grid->dzf[ne][k];
                } else{
                    sum += AT[ii][nf]*u[ne][grid->Nke[ne]-1];
                    //sum += AT[ii][nf]*u[ne][grid->Nke[ne]-1]*grid->dzf[ne][grid->Nke[ne]-1];
                }
            }
            bpr[ii]=sum;
        }
        // Now solve the problem (A'A)x = (A'b)
        linsolve(Apr,bpr,2);
       
        uc[n][k]=bpr[0];
        vc[n][k]=bpr[1];        
        //if (grid->dzz[n][k]  > 1e-3) {
        //    uc[n][k]=bpr[0]/grid->dzz[n][k];
        //    vc[n][k]=bpr[1]/grid->dzz[n][k];
        //}else{
        //    uc[n][k]=0;
        //    vc[n][k]=0;
        //}
    }

  }
}


/*
 * Function: ComputeQuadraticInterp
 * Usage: ComputeQuadraticInterp(U, V, un, grid, ninterp, tinterp, myproc)
 * -------------------------------------------
 * Compute the quadratic interpolation of the velocity based on
 * a choice for tinterp.
 * Two options presently exist based on choice of tinterp:
 *  1. tRT1
 *  2. tRT2
 * which are outlined in Wang et al, 2011.
 *
 */
static void  ComputeQuadraticInterp(REAL x, REAL y, int ic, int ik, REAL **uc, 
    REAL **vc, physT *phys, gridT *grid, interpolation ninterp, 
    interpolation tinterp, int myproc) {
  // this formulation considers a specific cell ic and a specific level ik
  // for computation of uc and vc

  // we need to have several pieces for this to work: 
  // we need to be able to quickly compute the area subsets in the triangle
  // A12, A13, A23; we need nodal velocity vectors (nRT2), and tangential
  // velocity vectors
  REAL points[2][grid->nfaces[ic]];
  REAL xt[3], yt[3];
  REAL SubArea[grid->nfaces[ic]];
  REAL nu[grid->nfaces[ic]], nv[grid->nfaces[ic]];
  REAL eu[grid->nfaces[ic]], ev[grid->nfaces[ic]];
  int np[grid->nfaces[ic]], ne[grid->nfaces[ic]], nf, ie, ip;
  const REAL TotalArea = grid->Ac[ic];

  // loop over each of the faces 
  // possible error if NUMEDGECOLUMNS != NFACES!!
  for( nf=0; nf < grid->nfaces[ic]; nf ++) {
    // get the index for the vertex of the cell
    np[nf] = grid->cells[grid->maxfaces*ic + nf];
    // get the index for the cell edge
    ne[nf] = grid->face[grid->maxfaces*ic + nf];

    // get the coordinates for the points for the vertex
    // to compute the area of the subtriangle
    points[0][nf] = grid->xp[np[nf]];
    points[1][nf] = grid->yp[np[nf]];
  }
  // now we have the points for each node: points[x/y][nf=0,1,2]
  // and also the edges via grid->face[NFACES*ic + nf=0,1,2] for
  // each cell ic

  // now all that remains is to store the velocities and then
  // get A12, A13, A23 and we can 
  // use the formula after getting the vector values for the 
  // velocities at each of the points

  for( nf=0; nf < grid->nfaces[ic]; nf++) {
    // get the points for each vertex of the triangle
    xt[0] = points[0][nf];
    yt[0] = points[1][nf];
    // wrap the next one if it goes out of bounds
    xt[1] = points[0][(nf+1)%grid->nfaces[ic]]; 
    yt[1] = points[1][(nf+1)%grid->nfaces[ic]];
    // the interpolated center
    xt[2] = x;
    yt[2] = y;

    // get the area from the points where 
    // nf=0 => A12, nf=1 => A23, nf=2 => A13 from Wang et al 2011
    // note that this is normalized
    SubArea[nf] = GetArea(xt, yt, 3)/TotalArea;
    //          printf("nf = %d Area=%e\n", nf, SubArea[nf]);
  }

  // get the nodal and edge velocities projected unto global x,y coords
  for (nf=0; nf < grid->nfaces[ic]; nf++) {
    // node and edge indexes
    ip = np[nf]; ie = ne[nf];
    // get the nodal velocity components 
    nu[nf] = phys->nRT2u[ip][ik];
    nv[nf] = phys->nRT2v[ip][ik];
    // get the tangential velocity components
    if(tinterp == tRT2) {
      eu[nf] = phys->u[ie][ik]*grid->n1[ie] + phys->tRT2[ie][ik]*grid->n2[ie];
      ev[nf] = phys->u[ie][ik]*grid->n2[ie] - phys->tRT2[ie][ik]*grid->n1[ie];
    }
    else if(tinterp == tRT1) {
      // need to have the specific cell neighbor here!  Current implementation 
      // is not correct.  This seems like the easiest way to make sure 
      // that this works
      eu[nf] = phys->u[ie][ik]*grid->n1[ie] + phys->tRT1[ie][ik]*grid->n2[ie];
      ev[nf] = phys->u[ie][ik]*grid->n2[ie] - phys->tRT1[ie][ik]*grid->n1[ie];
    }
  }

  // now perform interpolation from the results via Wang et al 2011 eq 12
  uc[ic][ik] = 
    (2*SubArea[1]-1)*SubArea[1]*nu[0] 
    + (2*SubArea[2]-1)*SubArea[2]*nu[1]
    + (2*SubArea[0]-1)*SubArea[0]*nu[2]
    + 4*SubArea[2]*SubArea[1]*eu[0]
    + 4*SubArea[0]*SubArea[2]*eu[1]
    + 4*SubArea[0]*SubArea[1]*eu[2];
  vc[ic][ik] = 
    (2*SubArea[1]-1)*SubArea[1]*nv[0] 
    + (2*SubArea[2]-1)*SubArea[2]*nv[1]
    + (2*SubArea[0]-1)*SubArea[0]*nv[2]
    + 4*SubArea[2]*SubArea[1]*ev[0]
    + 4*SubArea[0]*SubArea[2]*ev[1]
    + 4*SubArea[0]*SubArea[1]*ev[2];

}

/*
 * Function: ComputeTangentialVelocity
 * Usage: ComputeTangentialVelocity(phys, grid, ninterp, tinterp, myproc)
 * -------------------------------------------
 * Compute the tangential velocity from nodal velocities (interp determines type)
 * Two options presently exist based on choice of tinterp:
 *  1. tRT1
 *  2. tRT2
 * which are outlined in Wang et al, 2011.
 *
 */
static void  ComputeTangentialVelocity(physT *phys, gridT *grid, 
    interpolation ninterp, interpolation tinterp, int myproc) {
  int ie, in, ic, tn, ink;
  int nodes[2], cells[2];
  REAL tempu, tempv, tempA; 
  int tempnode, tempcell;


  if(tinterp == tRT2) {
    // for each edge, we must compute the value of its tangential velocity 
    for(ie = 0; ie < grid->Ne; ie++) {

      // get the nodes on either side of the edge
      nodes[0] = grid->edges[NUMEDGECOLUMNS*ie];
      nodes[1] = grid->edges[NUMEDGECOLUMNS*ie+1];
      //     printf("nodes = %d %d\n", nodes[0], nodes[1]);

      // for each layer at the edge
      for (ink = 0; ink < grid->Nkc[ie]; ink++) {
        // initialize temporary values
        tempu = tempv = tempA = 0;

        // for each node
        for (in = 0; in < 2; in++) {
          // get particular value for node 
          tempnode = nodes[in];

          // ensure that the node exists at the level 
          // before continuing
          if(ink < grid->Nkp[tempnode]) {
            // accumulate the velocity area-weighted values and 
            // cell area
            tempA += grid->Actotal[tempnode][ink];
            tempu += grid->Actotal[tempnode][ink]*phys->nRT2u[tempnode][ink];
            tempv += grid->Actotal[tempnode][ink]*phys->nRT2v[tempnode][ink];
          }
        }
        // area-average existing nodal velocities
        if(tempA == 0) {
          phys->tRT2[ie][ink] = 0;
        }
        else {
          tempu /= tempA;
          tempv /= tempA;
          // get the correct component in the tangential direction to store it
          phys->tRT2[ie][ink] = grid->n2[ie]*tempu - grid->n1[ie]*tempv;
        }
      }

    }

  }
}

/*
 * Function: ComputeNodalVelocity
 * Usage: ComputeNodalVelocity(phys, grid, interp, myproc)
 * -------------------------------------------
 * Compute the nodal velocity using RT0 basis functions.  
 * Two options presently exist based on choice of interp:
 *  1. nRT1
 *  2. nRT2
 * which are outlined in Wang et al, 2011.
 *
 */
static void ComputeNodalVelocity(physT *phys, gridT *grid, interpolation interp, int myproc) {
  //  int in, ink, e1, e2, cell, cp1, cp2;
  int in, ink, inpc, ie, intemp, cell, cp1, cp2,
      e1, e2, n1, n2, onode;
  REAL tempu, tempv, Atemp, tempAu, tempAv;


  /* compute the nodal velocity for RT1 elements */
  // for each node compute its nodal value
  for(in = 0; in < grid->Np; in++) {/*{{{*/
    // there are Nkp vertical layers for each node, so much compute over each 
    // of these

    for(ink = 0; ink < grid->Nkp[in]; ink++) {
      // there will be numpcneighs values for each node at a
      //particular layer so we must compute each separately

      if(interp == nRT2) { Atemp = tempAu = tempAv = 0; }
      for(inpc = 0; inpc < grid->numpcneighs[in]; inpc++) {

        // check to make sure that the cell under consideration exists at the given z-level
        if (ink < grid->Nk[grid->pcneighs[in][inpc]]) { 
          // if so get neighbors to the node
          e1 = grid->peneighs[in][2*inpc];
          e2 = grid->peneighs[in][2*inpc+1];

          // compute the RT0 reconstructed nodal value for the edges
          ComputeRT0Velocity(&tempu, &tempv, 
              grid->n1[e1], grid->n2[e1], grid->n1[e2], grid->n2[e2], phys->u[e1][ink], phys->u[e2][ink]);

          // store the computed values
          phys->nRT1u[in][ink][inpc] = tempu;
          phys->nRT1v[in][ink][inpc] = tempv;
          if(interp == nRT2) {
            // accumulate area and weighted values
            Atemp += grid->Ac[grid->pcneighs[in][inpc]];
            tempAu += grid->Ac[grid->pcneighs[in][inpc]]*tempu;
            tempAv += grid->Ac[grid->pcneighs[in][inpc]]*tempv;
          }

        } 
        else { // cell neighbor doesn't exist at this level
          phys->nRT1u[in][ink][inpc] = 0;
          phys->nRT1v[in][ink][inpc] = 0;
        }
      }
      // compute nRT2 from nRT1 values if desired  (area-weighted average of all cells)
      if(interp == nRT2) {
        if(Atemp == 0) {
          printf("Error as Atemp is 0 in nodal calc!! at in,ink,Nkp=%d,%d,%d\n",in, ink, grid->Nkp[in]);
          // print all nodal neighbors
          printf("cell neighbors = ");
          for(inpc = 0; inpc < grid->numpcneighs[in]; inpc++) {
            printf(" %d(%d)", grid->pcneighs[in][inpc], grid->Nk[grid->pcneighs[in][inpc]]);

          }
          printf("\n");

        }
        /* now we can compute the area-averaged values with nRT1 elements */
        phys->nRT2u[in][ink] = tempAu/Atemp;
        phys->nRT2v[in][ink] = tempAv/Atemp;
      }
    }
  }
}

/*
 * Function: ComputeRT0Velocity
 * Usage: ComputeRT0Velocity(REAL* tempu, REAL* tempv, int e1, int e2, physT* phys);
 * -------------------------------------------
 * Compute the nodal velocity using RT0 basis functions (Appendix B, Wang et al 2011)
 *
 */
static void ComputeRT0Velocity(REAL *tempu, REAL *tempv, REAL e1n1, REAL e1n2, 
    REAL e2n1, REAL e2n2, REAL Uj1, REAL Uj2) 
{
  const REAL det = e1n1*e2n2 - e1n2*e2n1;
  // compute using the basis functions directly with inverted 2x2 matrix
  *tempu = (e2n2*Uj1 - e1n2*Uj2)/det;
  *tempv = (e1n1*Uj2 - e2n1*Uj1)/det;
}

/*
 * Function: HFaceFlux
 * Usage: HFaceFlux(j,k,phi,phys->u,grid,prop->dt,prop->nonlinear);
 * ---------------------------------------------------------------------
 * Interpolation to obtain the flux of the scalar field phi (of type REAL **) on 
 * face j, k;  method==2: Central-differencing, method==4: Lax-Wendroff.
 *
 */
// note that we should not compute def1 and def2 in this function as they don't change
// these should be computed in grid.c and just looked up for efficiency
static REAL HFaceFlux(int j, int k, REAL *phi, REAL **u, gridT *grid, REAL dt, int method) {
  int nc1, nc2;
  REAL def1, def2, Dj, C=0;
  nc1 = grid->grad[2*j];
  nc2 = grid->grad[2*j+1];
  if(nc1==-1) nc1=nc2;
  if(nc2==-1) nc2=nc1;
  Dj = grid->dg[j];
  Return_def(&def1,&def2,nc1,nc2,j,grid);	   	          
  //  def1=grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
  //  def2=Dj-def1;

  if(method==4) C = u[j][k]*dt/Dj;

  if(def1==0 || def2==0) {
    // this happens on a boundary cell
    return UpWind(u[j][k],phi[nc1],phi[nc2]);
  }
  else {
    // on an interior cell (orthogonal cell) we can just distance interpolation the 
    // value to the face)
    // basically Eqn 51
    return (phi[nc1]*def2+phi[nc2]*def1)/(def1+def2)
      -C/2*(phi[nc1]-phi[nc2]);
  }
}

/*
 * Function: getTsurf
 * --------------------------------------------
 * Returns the temperature at the surface cell
 */
static void getTsurf(gridT *grid, physT *phys){
  int i, iptr, ktop;
  int Nc = grid->Nc;
  
  //for(i=0;i<Nc;i++) {
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    ktop = grid->ctop[i];
    phys->Tsurf[i] = phys->T[i][ktop];
  }
}

/*
 * Function: getchangeT
 * ----------------------------------------------------------------
 * Returns the change in temperature at the surface cell
 */
static void getchangeT(gridT *grid, physT *phys){
  int i, iptr, ktop;
  int Nc = grid->Nc;
  
  //for(i=0;i<Nc;i++) {
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    ktop = grid->ctop[i];
    phys->dT[i] = phys->T[i][ktop] - phys->Tsurf[i];
  }
}
 /*
 *
 */
static void GetMomentumFaceValues(REAL **uface, REAL **ui, REAL **boundary_ui, REAL **U, gridT *grid, physT *phys, propT *prop,
				  MPI_Comm comm, int myproc, int nonlinear, int TVD) {
  int i, iptr, j, jptr, nc, nf, ne, nc1, nc2, k, kmin;
  REAL tempu, u_im;

  //U holds phys->u most of time, but want to maintain the current time step U for W. 

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];
    i = grid->grad[2*j];

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      u_im = prop->imfac1*U[j][k]+prop->imfac2*phys->u_old[j][k]+prop->imfac3*phys->u_old2[j][k];
      if(u_im>0)
        uface[j][k]=boundary_ui[jptr-grid->edgedist[2]][k];
      else
        uface[j][k]=ui[i][k];
    }
  }

  // type 4 boundary conditions used for no slip 
  // but typically we think of no-slip boundary conditions as 
  // not having flux across them, so this is more general
  for(jptr=grid->edgedist[4];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      uface[j][k]=boundary_ui[jptr-grid->edgedist[2]][k];
  }

  if(nonlinear==5){ //use tvd for advection of momemtum
    HorizontalFaceScalars(grid,phys,prop,ui,U,boundary_ui,TVD,3,comm,myproc);
  }
  // over each of the "computational" cells
  // Compute the u-component fluxes at the faces
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // figure out which column adjacent to edge is limiter for the
    // top of the cell (cells much share face to have flux)
    if(grid->ctop[nc1]>grid->ctop[nc2])
      kmin = grid->ctop[nc1];
    else
      kmin = grid->ctop[nc2];
    
    for(k=0;k<kmin;k++)
      uface[j][k]=0;

    // compute mass flow rate by interpolating the velocity to the face flux
    // note that we can use UFaceFlux to just compute the value directly if 
    // we have quadratic upwinding available to us
    // basically Eqn 47
    for(k=kmin;k<grid->Nke[j];k++) {
      // compute upwinding data
      if(U[j][k]>0)
        nc=nc2;
      else
        nc=nc1;

      switch(nonlinear) {
        case 1:
          uface[j][k]=ui[nc][k];
          break;
        case 2:
          uface[j][k]=UFaceFlux(j,k,ui,U,grid,prop->dt,nonlinear);
          break;
        case 4:
          uface[j][k]=UFaceFlux(j,k,ui,U,grid,prop->dt,nonlinear);
          break;
        case 5:
          //u_im = prop->imfac2*phys->u_old[j][k]+prop->imfac1*phys->u[j][k]+prop->imfac3*phys->u_old2[j][k];
          u_im = prop->imfac1*U[j][k]+prop->imfac2*phys->u_old[j][k]+prop->imfac3*phys->u_old2[j][k];
          //if(U[j][k]>0)
          if(u_im>0)
            tempu=phys->SfHp[j][k];
          else
            tempu=phys->SfHm[j][k];
          uface[j][k]=tempu;
          break;
        default:
          uface[j][k]=ui[nc][k];
          break;	  
      }
      if(nonlinear==2)
      uface[j][k]=UFaceFlux(j,k,ui,U,grid,prop->dt,2);      
          
    }
  }

  // Faces on type 3 cells are always updated with first-order upwind
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    for(nf=0;nf<grid->nfaces[i];nf++) {
      if((ne=grid->neigh[i*grid->maxfaces+nf])!=-1) {
        j=grid->face[i*grid->maxfaces+nf];

        nc1 = grid->grad[2*j];
        nc2 = grid->grad[2*j+1];

        if(grid->ctop[nc1]>grid->ctop[nc2])
          kmin = grid->ctop[nc1];
        else
          kmin = grid->ctop[nc2];

        for(k=0;k<kmin;k++)
          uface[j][k]=0;
        for(k=kmin;k<grid->Nke[j];k++) {
          u_im = prop->imfac1*U[j][k]+prop->imfac2*phys->u_old[j][k]+prop->imfac3*phys->u_old2[j][k];
          if(u_im>0)
            nc=nc2;
          else
            nc=nc1;
          uface[j][k]=ui[nc][k];
        }
      }
    }
  }
}

/*
 * Function: UpdateBottomHeight
 * Usage: UpdateBottomHeight(zB,zBold,zBold2,grid,prop,myproc,comm);
 * -----------------------------------------------------------------
 * Sets the bottom height at time-step n+1 (zB) and step n time (zBold), and n-1 (zBold2).
 *
 * Set balance_volume=1 to adjust zB so that it conserves volume. This is not necessary
 * since any nonzero change in volume due to zB will result in an change in the free-surface
 * height to conserve volume.
 *
 */
static int UpdateBottomHeight(REAL *zB, REAL *zBold, REAL *zBold2, REAL *zBoffline, gridT *grid,
			      propT *prop, physT *phys, int myproc, MPI_Comm comm) {
  int i, iptr, balance_volume=0; 
  REAL zb_sum, Ac_sum, my_zb_sum, my_Ac_sum, L;

  if(prop->n==1+prop->nstart) {
    
           
    //normal sine
    for(i=0;i<grid->Nc;i++) {

      L = prop->wavelength;
      
      zB[i]=.1*sin(2.0*PI/L*grid->xv[i]);
  
      zBold[i]=zB[i];
      zBold2[i]=zB[i];

      zBoffline[i]=zB[i];
    }
      
    

    ISendRecvCellData2D(zB,grid,myproc,comm);  
    ISendRecvCellData2D(zBold,grid,myproc,comm); 
    ISendRecvCellData2D(zBold2,grid,myproc,comm);
    ISendRecvCellData2D(zBoffline,grid,myproc,comm);
 }



  

  //let's try the Chou way
  //calc w_s 
  int nf, ne, nc1, nc2;
  REAL d_s = 0.0003; //300 micrometers
  //d_s = 0.001; //1000 micrometers
  //lets try a new one
  
  REAL g = 9.81; //don't want to use fake grav here I think
  REAL p0 = 0.4;
  REAL w_s, H, ue, ve, umag, u_b;
  REAL rho_s = 2650; // kg/m^3 for sand
  REAL G = (rho_s*g)/(RHO0*g); //specific gravity
  REAL ustar, tau_c_star, psi, q, dstar, tau_star, tau_b, tau_c0, dzbdn, phi, q_s, q_sx, q_sy, k;
  REAL delq;
  REAL C_s = 8;
  REAL phi_r=30*PI/180; //30 degrees in radians for angle of repose (artificially low angle) 
  REAL alpha_s = 0.0001;
  REAL a = 0.3; 
  REAL b = 2.0; 
  //REAL alpha_s = 2.5; 
  REAL starttime = 2830; //kind of arbitrary (was 3000)
  //starttime = 500000;
  REAL ramptime=5000/1.7; //ramp up over 5000 m or so? Estimating cwave as 1.7 here  
  REAL rampcenter = 2830;
  REAL nu = 1E-6; //use physical nu here instead of prop->nu?
  //REAL nu = prop->nu;
  REAL lambda_s = 1.0;
  //calc dstar, tau_cstar
  dstar = d_s*pow((G-1)*g/pow(nu, 2),(1.0/3));
  tau_c_star = 0.3*exp(-dstar/3) + 0.047*(1-exp(-dstar/20));
  //printf("tau_c is %f for x = %f \n", tau_c_star, grid->xv[i]);

  starttime = 0;

  //loop over cells? maybe edges instead?
  if(prop->rtime>=starttime){
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      H=phys->h[i]+grid->dv[i]; //should be averaged over 2 edges. 

      delq = 0;

      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        nc1 = grid->grad[2*ne];
        nc2 = grid->grad[2*ne+1]; 

        if(nc1==-1)
          nc1=nc2;
        if(nc2==-1)
          nc2=nc1;

	//u_b = prop->imfac1*phys->u[ne][grid->Nke[ne]-1]+prop->imfac2*phys->u_old[ne][grid->Nke[ne]-1]+prop->imfac3*phys->u_old2[ne][grid->Nke[ne]-1];
        u_b = phys->u[ne][grid->Nke[ne]-1];
	ue = u_b*grid->n1[ne];
        ve = u_b*grid->n2[ne];
        umag=sqrt(ue*ue + ve*ve);

	umag = sqrt(u_b*u_b);
	if(fabs(ve)>0)
	  printf("nonzero y vel! \n");

	if(prop->CdB==-1){
	  tau_b=RHO0*prop->nu*(2*umag)/(0.5*(grid->dzz[nc1][grid->Nk[nc1]-1]+grid->dzz[nc2][grid->Nk[nc2]-1]));
	}else if(prop->z0B==0){
	  tau_b = RHO0*0.0025*umag*umag;
	}else{
	  tau_b = RHO0*phys->CdB[ne]*umag*umag;  
	}
        //tau_b = RHO0*phys->CdB[ne]*umag*umag;
	//tau_b = RHO0*0.0025*umag*umag;
	//tau_b=RHO0*prop->nu*(2*umag)/(0.5*(grid->dzz[nc1][grid->Nk[nc1]-1]+grid->dzz[nc2][grid->Nk[nc2]-1]));

        //calc tau_star 
        tau_star = tau_b/((G-1)*RHO0*g*d_s);

        //calc critical tau 
        tau_c0 = 0.3/(1+1.2*dstar) + 0.055*(1-exp(-0.02*dstar));


        //calc slope of bed at this edge 
        //dzbdn = (zB[nc2]-zB[nc1])/grid->dg[ne]*grid->normal[ne]; //do we need the normal here? because no variation in x, this should give 0 for y and something for x.
        dzbdn = (zBoffline[nc2]-zBoffline[nc1])/grid->dg[ne]; //do we need the normal here? because no variation in x, this should give 0 for y and something for x.
        phi=(atan(dzbdn)); //for now just have in x direction since its 2D

         if(fabs(phi)>phi_r){
           k=alpha_s*(phi-phi_r)/phi_r; //had an abs on the first phi, but think we want directionality here.
         }else{
           k=0;
         }

        //try old way of constant decrease for slope term?
        //k = a*pow(fabs(tau_b/RHO0), b)*alpha_s;
         //k = 0.0;
	 //k = alpha_s*100; //*10 iso, /10 sig
	 //k = alpha_s/10.0;
	 k=0.0;

	//update critical tau for bed slope at edge
        //tau_c_star = tau_c0*sin(phi_r + phi)/sin(phi_r);
	tau_c_star = tau_c0;
	//tau_c_star = 0.0;

        //calc dimensional bedload transport (Meyer-Peter and Muller)
        if(tau_star > tau_c_star){
          q_s = 8.0*pow((tau_star-tau_c_star), 1.5)*(d_s*sqrt((G-1)*g*d_s));
	  //q_s = 8.0*pow((tau_star-tau_c_star), 1.0)*(d_s*sqrt((G-1)*g*d_s)); 
	} else {
          q_s = 0;
        }

	k = lambda_s * q_s;  //test just the other stuff

        if(fabs(umag)>0){
          q_sx = q_s*ue/umag;
          q_sy = q_s*ve/umag;

	  //q_sx = q_s;
        } else{
          q_sx = 0;
          q_sy = 0;
        }

	//q_sx = 0; //test just downslope transport!

	//ramp up qb? //turn off for now. 
	// q_sx=q_sx*.5*erf((prop->rtime-rampcenter)/ramptime);
	// q_sy=q_sy*.5*erf((prop->rtime-rampcenter)/ramptime);

        //we are going to ignore the n1, n2 part for now because we know its 2d but in theory should add component 1 *n1 and compoennt 2 *n2
        //q = psi * w_s * d_s;
	
	//+ k because of the negative overall?
	delq+=1.0/grid->Ac[i]*(q_sx*grid->n1[ne]+q_sy*grid->n2[ne] + k*dzbdn)*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];

        //delq+=1.0/grid->Ac[i]*(q_sx*grid->n1[ne]+q_sy*grid->n2[ne])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];
        //delq+=1.0/grid->Ac[i]*(q_sx)*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];                                                

      }

      //update zb
      //going in, want zB at n+1. zB contains n, zBold contains n-1, zBold2 contains n-2.
      
      zBold2[i]=zBold[i];
      zBold[i]=zB[i];
      //for zBoffline...
      zB[i] = zB[i];
      //printf("got here \n");
      //zBoffline[i]-=prop->dt/(1.0-p0)*100*delq;
      zBoffline[i]+=prop->dt*(1/(1-p0))*100*(phys->Deposition[i]-phys->Erosion[i])/rho_s;
      //zBoffline[i]=0;
      
      //remove init
      //zBoffline[i]-=prop->dt*(1/(1-p0))*100*delq + prop->dt*zBoffline[i]*(1/(10*prop->dt)*(0.5*(1.0-tanh((prop->rtime-rampcenter)/ramptime))));
      //zBoffline[i]-=prop->dt*(1/(1-p0))*100*delq + zB[i]*(.5*tanh((prop->rtime-rampcenter)/-ramptime));

      //unoffline it
      zB[i] = zBoffline[i];


      //now, zB contains n, zBold contains n, zBold2 contains n-1.
      
      //usual update step
      //zB[i]-=prop->dt*(1/(1-p0))*100*delq;
      
      //dzdt = (1.5*phys->zB[i] - 2.0*phys->zBold[i] + 0.5*phys->zBold2[i])/prop->dt;
      //zB[i] = (2.0*phys->zBold[i] - 0.5*phys->zBold2[i] - prop->dt*(1/(1-p0))*100*delq)/1.5;
      //zB[i]=zBold2[i] - 2*prop->dt*100*(1/(1-p0))*delq;
    }
  } else{
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      zBold2[i]=zBold[i];
      zBold[i]=zB[i];
      zB[i]=zB[i];

      zBoffline[i]=zB[i];
    }
  }
  

  //try borsje 

  //u_cr = 0.19*pow(d_s, 0.1)*log10(4*(phys->dv)/d_s);


  /*
  int nf, ne, nc1, nc2;
  REAL *dzbdx, *dzbdy;
  REAL zbe, umag;

  REAL ue, ve, dzz_bot, dzbdx_e, dzbdy_e, taubx, tauby, taub, taubn, qbx, qby, qbn, delq, delq_a, delq_t;
  REAL a = 0.3; 
  REAL b = 2.0; 
  REAL alpha_s = 2.5; 
  REAL p0=0.4;
  REAL phi_r=PI/30;
  REAL phi_r=30*PI/180;
  REAL phi; 
  REAL k, dhdx, dzbdn, dzbdn_x, dzbdn_y;

  //REAL rampcenter=5000/prop->cwave; //ramp up over 2500 m?
  //REAL ramptime=1000/prop->cwave; //ramp up over 2500 m?  
  //compute bed changes after first tidal period 
  //REAL starttime=(6*PI/prop->omega); //2 for test
  REAL starttime=(0); //start immediately
  //REAL starttime=10000/(prop->dt*prop->cwave)/2; //start after three periods
  //REAL ramptime=5000/prop->cwave; //ramp up over 2500 m? 
  //REAL rampcenter = starttime;

  if(prop->rtime>=starttime){
  //if(0){
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; //center pointer
    delq=0; //start sum at 0
    delq_a=0;
    delq_t=0;
        
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf]; //edge index of face
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1]; 

      if(nc1==-1)
        nc1=nc2;
      if(nc2==-1)
        nc2=nc1;

      //calculate edge velocities 
      // ue=0.5*(phys->uc[nc1][grid->Nk[nc1]-1]+phys->uc[nc2][grid->Nk[nc2]-1]); 
      // ve=0.5*(phys->vc[nc1][grid->Nk[nc1]-1]+phys->vc[nc2][grid->Nk[nc2]-1]); 

      ue = phys->u[ne][grid->Nke[ne]-1]*grid->n1[ne];
      ve = phys->u[ne][grid->Nke[ne]-1]*grid->n2[ne];

      umag=sqrt(ue*ue + ve*ve);
      if (umag != fabs(phys->u[ne][grid->Nke[ne]-1])){
        printf("umag = %f is not equal to mag phys->u = %f \n", umag, fabs(phys->u[ne][grid->Nke[ne]-1]));
      }

      //calculate dz at edge and interpolate slope terms onto edge
      //note that this is at previous time step...
      dzz_bot=0.5*(grid->dzz[nc1][grid->Nk[nc1]-1]+grid->dzz[nc2][grid->Nk[nc2]-1]);
     
      //calculate bottom stresses
      // taubx=prop->nu*(2*ue)/dzz_bot;
      // tauby=prop->nu*(2*ve)/dzz_bot;

      //calculate bottom stresses
      taubx=prop->CdB*ue*ue;
      tauby=prop->CdB*ve*ve;
      
      //taubx=phys->CdB[ne]*ue*fabs(ue);
      //tauby=phys->CdB[ne]*ve*fabs(ve);

      //taub = phys->CdB[ne]*phys->u[ne][grid->Nke[ne]-1]*fabs(phys->u[ne][grid->Nke[ne]-1]);

      taub = prop->CdB*phys->u[ne][grid->Nke[ne]-1]*phys->u[ne][grid->Nke[ne]-1];
      //taub = 2*prop->nu*fabs(phys->u[ne][grid->Nke[ne]-1])/dzz_bot; //total stress
      if(taub==0){
        taubn=0;
      }else{
        taubn=taub*phys->u[ne][grid->Nke[ne]-1]/umag;
      }

      //calculate qbx, qby
      // if(taubx==0){
      //   qbx=0;
      // }else{
      //   qbx=a*pow(fabs(taubn), b)*(taubx/fabs(taubn)-alpha_s*dzbdx_e);
      // }
      
      // if(tauby==0){
      //   qby=0;
      // }else{
      //   qby=a*pow(fabs(taubn), b)*(tauby/fabs(taubn)-alpha_s*dzbdy_e);      
      // }

      if(taubn==0){
            qby=0;
        qbx=0;
      }else{
        qbx=a*pow(fabs(taubn), b)*(taubx/fabs(taubn));
	      qby=a*pow(fabs(taubn), b)*(tauby/fabs(taubn));    
        qbn=a*pow(fabs(taubn), b)*phys->u[ne][grid->Nke[ne]-1]/umag;
      }

      dzbdn = (zB[nc1]-zB[nc2])/grid->dg[ne]*grid->normal[ne]; //do we need the normal here?
      dzbdn_x = dzbdn*fabs(ue)/umag;
      dzbdn_y = dzbdn*fabs(ve)/umag;
      phi=fabs(atan(dzbdn)); //for now just have in x direction
      //printf("phi is %f at %f = xv \n", phi, (grid->xv[nc1]+grid->xv[nc2])/2);

      if(fabs(phi)>phi_r){
        k=alpha_s*(fabs(phi)-phi_r)/phi_r;
      }else{
        k=0;
      }

      //ramp up qb
      //qbx=qbx*.5*erf((prop->rtime-rampcenter)/ramptime);
      //qby=qby*.5*erf((prop->rtime-rampcenter)/ramptime);

      //qbx=alpha_s*dzbdn;
      //qby=alpha_s*0;
      

      //qbx=-alpha_s*dzbdx_e;
      //qby=-alpha_s*dzbdy_e;

      //does this need a 1/dg?
      //alpha_s=0;
      //printf("nc1 is %d and nc2 is %d, which makes xv1 %f and xv2 %f \n", nc1, nc2, grid->xv[nc1], grid->xv[nc2]);
      //printf("zB[nc1] is %f, zB[nc2] is %f, dg is %f, normal is %d, dg is %f \n", zB[nc1], zB[nc2], grid->dg[ne], grid->normal[ne], grid->df[ne]);
      //delq+=alpha_s*(zB[nc1]-zB[nc2])*grid->normal[ne]*grid->df[ne];
      //delq+=alpha_s*(zB[nc1]-zB[nc2])/grid->dg[ne]*grid->normal[ne]*grid->df[ne];
      
      //add to sum for gradient term 
      //delq+=(qbx*grid->n1[ne]+qby*grid->n2[ne])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];
      //delq+=qbn*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];
      //delq+=qbn*grid->df[ne];
      
      delq_a+=a*pow(fabs(taubn), b)*alpha_s*(zB[nc1]-zB[nc2])/grid->dg[ne]*grid->normal[i*grid->maxfaces+nf]*grid->df[ne];
      
      //delq+=alpha_s*(zB[nc1]-zB[nc2])/grid->dg[ne]*grid->normal[i*grid->maxfaces+nf]*grid->df[ne]; //does this need a dg? I think yes. check if same with higher res if not
      //delq_a+=k*(zB[nc1]-zB[nc2])/grid->dg[ne]*grid->normal[i*grid->maxfaces+nf]*grid->df[ne]; 
      delq_t+=(qbx*grid->n1[ne]+qby*grid->n2[ne])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];
      //delq_t=0; 
    }
    delq_a/=grid->Ac[i];
    delq_t/=grid->Ac[i];
    //delq/=grid->Ac[i];
    delq=delq_t-delq_a;
    //delq=delq_t;

    //ramp it 
    //delq=delq*(1-exp(-prop->rtime/ramptime)); 
    
    //update zb
    zBold2[i]=zBold[i];
    zBold[i]=zB[i];
    zB[i]-=prop->dt*(1/(1-p0))*10*delq;
    
    
  }
  }else{
      zB[i]=zB[i];
      zBold2[i]=zBold[i];
      zBold[i]=zB[i];
      //printf("set new zb successfully \n");
    }
  */



  // for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
  //   i = grid->cellp[iptr];
    
  //   zBold2[i]=zBold[i];
  //   zBold[i]=zB[i];
  //   //zB[i]=prop->amp*exp(-pow((grid->yv[i]-50)/5,2))*sin(prop->omega*prop->rtime);    
  //   zB[i]=prop->amp*cos(2.0*PI*grid->xv[i]/1000)*sin(prop->omega*prop->rtime);
  // }

  if(balance_volume) {

    my_zb_sum = 0;
    my_Ac_sum = 0;
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      my_zb_sum+=grid->Ac[i]*zB[i];
      my_Ac_sum+=grid->Ac[i];
    }

    MPI_Reduce(&my_zb_sum,&(zb_sum),1,MPI_DOUBLE,MPI_SUM,0,comm);
    MPI_Bcast(&zb_sum,1,MPI_DOUBLE,0,comm);

    MPI_Reduce(&my_Ac_sum,&(Ac_sum),1,MPI_DOUBLE,MPI_SUM,0,comm);
    MPI_Bcast(&Ac_sum,1,MPI_DOUBLE,0,comm);        

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      zB[i]-=(zb_sum/Ac_sum);
    }
  }

  ISendRecvCellData2D(zBold2,grid,myproc,comm);      
  ISendRecvCellData2D(zBold,grid,myproc,comm);    
  ISendRecvCellData2D(zB,grid,myproc,comm);
  ISendRecvCellData2D(zBoffline,grid,myproc,comm);  
}