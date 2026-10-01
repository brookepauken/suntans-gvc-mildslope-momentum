/*
 * File: upredictor.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains the hydrostatic predictor that solves for the free surface.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "upredictor.h"
#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "boundaries.h"
#include "sediments.h"
#include "marsh.h"
#include "vertcoordinate.h"
#include "culvert.h"
#include "subgrid.h"

static void BiCGSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void CGSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void HPreconditioner(REAL *x, REAL *y, gridT *grid, physT *phys, propT *prop);
static void HCoefficients(REAL *coef, REAL *fcoef, gridT *grid, physT *phys, 
    propT *prop);
static REAL InnerProduct(REAL *x, REAL *y, gridT *grid, int myproc, int numprocs, 
    MPI_Comm comm);
static void OperatorH(REAL *x, REAL *y, REAL *coef, REAL *fcoef, gridT *grid, 
    physT *phys, propT *prop);

/*
 * Function: UPredictor 
 * Usage: UPredictor(grid,phys,prop,myproc,numprocs,comm);
 * -------------------------------------------------------
 * Predictor step for the horizontal velocity field.  This function
 * computes the free surface using the theta method and then uses 
 * it to update the predicted
 * velocity field in the absence of the nonhydrostatic pressure.
 *
 * Upon entry, phys->utmp contains the right hand side of the u-momentum equation
 *
 */
void UPredictor(gridT *grid, physT *phys, 
		propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, iptr, j, jptr, ne, nf, nf1, normal, nc1, nc2, k, n0, n1,iv,jv, botinterp,flag, Nkeb;
  REAL hmax,sum,sum0,sum1, dt=prop->dt, theta=prop->theta, h0, boundary_flag,fac1,fac2,fac3,tmp,tmp_x,tmp_y,tmp2;
  REAL *a, *b, *c, *d, *e1, **E, *a0, *b0, *c0, *d0, theta0, alpha,min,min2,def1,def2,dgf,l0,l1,zfb;
  REAL *dzm, Hf;
  
  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;
  e1 = phys->ap;
  E = phys->ut;
  dzm = phys->e;  

  a0 = phys->am;
  b0 = phys->bp;
  c0 = phys->bm;
  
  fac1=prop->imfac1; //theta (n+1/*)
  fac2=prop->imfac2; //(1-theta) (n)
  fac3=prop->imfac3; //third one (n-1)



  if(!prop->readOldVelocity){
  if(prop->n==1) {
    for(j=0;j<grid->Ne;j++)
      for(k=0;k<grid->Nke[j];k++){
        phys->u_old2[j][k]=phys->u[j][k];
      }
    for(i=0;i<grid->Nc;i++) 
      for(k=0;k<grid->Nk[i]+1;k++) 
      {
        phys->w_old2[i][k]=phys->w[i][k];
        phys->w_im[i][k]=0;
      }
    for(i=0;i<grid->Nc;i++)
      phys->h_old[i]=phys->h[i]; 
  }
  }else {
    if(prop->n==1) {
       for(i=0;i<grid->Nc;i++) 
         for(k=0;k<grid->Nk[i]+1;k++) 
         {
           phys->w_im[i][k]=fac2*phys->w_old[i][k]+fac3*phys->w_old2[i][k]+fac1*phys->w[i][k];;
         }

      for(i=0;i<grid->Nc;i++)
        phys->h_old[i]=phys->h[i]; 

      ISendRecvCellData2D(phys->h,grid,myproc,comm);
      ISendRecvCellData2D(phys->h_old,grid,myproc,comm);
    }
  }

  // Set D[j] = 0 
  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i]+1;k++) 
      phys->w_old[i][k]=phys->w[i][k];

  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=0;
    for(k=0;k<grid->Nke[j];k++)
      phys->u_old[j][k]=phys->u[j][k];
  }

  // note that phys->utmp is the horizontalsource term computed 
  // in HorizontalSource for 1/2(3F_j,k^n -F_j,k^n-1).
  // can be AB3

  // phys->u contains the velocity specified at the open boundaries
  // It is also the velocity at time step n. (type 2)
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    // transfer boundary flux velocities onto the horizontal source
    // term (exact as specified)
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->utmp[j][k]=phys->u[j][k];
  }

  // Update the velocity in the interior nodes with the old free-surface gradient
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    // Add the explicit part of the free-surface to create U**.
    // 5th term of Eqn 31
    for(k=grid->etop[j];k<grid->Nke[j];k++){
      phys->utmp[j][k]-=
        prop->grav*dt*(fac2*(phys->h[nc1]-phys->h[nc2])+fac3*(phys->h_old[nc1]-phys->h_old[nc2]))/grid->dg[j];
    }
  }

  // Drag term must be fully implicit
  theta0=theta;
  //theta=1; //trigger compile here

  // fac1=1;
  // fac2=0;
  // fac3=0;

  // Advection term for vertical momentum.  When alpha=1, first-order upwind,
  // alpha=0 is second-order central.  Always do first-order upwind when doing
  // vertically-implicit momentum advection
  //alpha=1;
  //alpha=0;
  if(prop->nonlinear==0)
    alpha=1;
  else
    alpha=0;

  // for each of the computational edges
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;
    
    def1 = grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
    def2 = grid->def[nc2*grid->maxfaces+grid->gradf[2*j+1]];

    // Flux heights for momentum (dzm) may not be the same as for volume (dzf) since
    // volume flux is upwinded whereas momentum is not necessarily (e.g. diffusion).
    Hf=0;
    for(k=0;k<grid->Nke[j];k++) {
      dzm[k]=grid->dzf[j][k];
      if(dzm[k]==0)
	dzm[k]=Max(grid->dzz[nc1][k],grid->dzz[nc2][k]);
      Hf+=dzm[k];
    }
    
    // Add the wind shear stress from the top cell
    phys->utmp[j][grid->etop[j]]+=2.0*dt*phys->tau_T[j]/
      (2.0*dzm[grid->etop[j]]);

    // define the bottom layer
    if(prop->vertcoord==1 || phys->CdB[j]==-1)
      Nkeb=grid->Nke[j]-1;
    else
      Nkeb=vert->Nkeb[j];

    // Create the tridiagonal entries and formulate U***
    // provided that we don't have a zero-depth top cell
    if(dzm[grid->etop[j]]!=0) {
      // initialize coefficients
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
        a[k]=0;
        b[k]=0;
        c[k]=0;
        d[k]=0;
      }

      // Vertical eddy-viscosity interpolated to faces since it is stored
      // at cell-centers.
      for(k=grid->etop[j]+1;k<grid->Nke[j];k++){ 
        c[k]=0.25*(phys->nu_tv[nc1][k-1]+phys->nu_tv[nc2][k-1]+
            phys->nu_tv[nc1][k]+phys->nu_tv[nc2][k]+
            prop->laxWendroff_Vertical*(phys->nu_lax[nc1][k-1]+phys->nu_lax[nc2][k-1]+
              phys->nu_lax[nc1][k]+phys->nu_lax[nc2][k]));
      }

      // Coefficients for the viscous terms.  Face heights are taken as
      // the average of the face heights on either side of the face (not upwinded).
      for(k=grid->etop[j]+1;k<grid->Nke[j];k++) 
        a[k]=2.0*(prop->nu+c[k])/(dzm[k]*(dzm[k-1]+dzm[k]));

      for(k=grid->etop[j];k<grid->Nke[j]-1;k++) {
        b[k]=2.0*(prop->nu+c[k+1])/(dzm[k]*(dzm[k]+dzm[k+1]));
      }
      // Added this on 030824 to correctly predict the bottom stress
      if(grid->etop[j]<grid->Nke[j]-1) {
	a[grid->etop[j]]=a[grid->etop[j]+1];
	b[grid->Nke[j]-1]=b[grid->Nke[j]-2];
      }

      // Coefficients for vertical momentum advection terms
      // d[] stores vertical velocity interpolated to faces vertically half-way between U locations
      // So d[k] contains w defined at the vertical w-location of cell k
      if(prop->vertcoord==1 && prop->nonlinear && 
        prop->thetaM>=0 && grid->Nke[j]-grid->etop[j]>1) {
        if(grid->ctop[nc1]>grid->ctop[nc2]) {
          n0=nc2;
          n1=nc1;
        } else {
          n0=nc1;
          n1=nc2;
        }
        // Don't do advection on vertical faces without water on both sides.
        for(k=0;k<grid->ctop[n1];k++)
          d[k]=0;
        for(k=grid->ctop[n1];k<grid->Nke[j];k++)
          d[k] = 0.5*(phys->w[n0][k]+phys->w[n1][k]);
        d[grid->Nke[j]]=0; // Assume w=0 at a corners (even if w is nonzero on one side of the face)
        for(k=grid->etop[j];k<grid->Nke[j];k++) {
          a0[k] = (alpha*0.5*(d[k]-fabs(d[k])) + 0.5*(1-alpha)*d[k])/(0.5*(2.0*dzm[k]));
          b0[k] = (alpha*0.5*(d[k]+fabs(d[k])-d[k+1]+fabs(d[k+1]))+0.5*(1-alpha)*(d[k]-d[k+1]))/(0.5*(2.0*dzm[k]));
          c0[k] = -(alpha*0.5*(d[k+1]+fabs(d[k+1])) + 0.5*(1-alpha)*d[k+1])/(0.5*(2.0*dzm[k]));			
        }
      }

      // add the vertical momentum advection part for the new general vertical coordinate
      // the vertical momentum is divided into two parts
      // wdudz=d(wu)dz-udwdz
      // the first part is the same as the original method while exchange w into omega
      if(prop->vertcoord!=1 && prop->nonlinear && grid->Nke[j]-grid->etop[j]>1 && prop->thetaM>=0)
      {
        // may be useless for vertical coordinate without z-level
        if(grid->ctop[nc1]>grid->ctop[nc2]) {
          n0=nc2;
          n1=nc1;
          l0=def2;
          l1=def1;
        } else {
          n0=nc1;
          n1=nc2;
          l0=def1;
          l1=def2;
        }

        // Don't do advection on vertical faces without water on both sides.
        for(k=0;k<grid->ctop[n1];k++)
          d[k]=0;
        // omega[top]=0 omega[bot]=0
        for(k=grid->ctop[n1];k<grid->Nke[j];k++)
          d[k] =(l1*vert->omega_old[n0][k]+l0*vert->omega_old[n1][k])/grid->dg[j];
        d[grid->Nke[j]]=0; // Assume w=0 at a corners (even if w is nonzero on one side of the face)
        for(k=grid->etop[j];k<grid->Nke[j];k++) 
        {
	  /*
          a0[k] = (d[k]-fabs(d[k]))/(2.0*dzm[k]);
          b0[k] = (d[k]+fabs(d[k])-d[k+1]+fabs(d[k+1]))/(2.0*dzm[k]);
          c0[k] = -(d[k+1]+fabs(d[k+1]))/(2.0*dzm[k]);                        
	  */
	  // This is the original method which allowed first-order upwind (alpha=1) or central differencing (alpha=0)
	  alpha=0;
          a0[k] = (alpha*0.5*(d[k]-fabs(d[k])) + 0.5*(1-alpha)*d[k])/(0.5*(2.0*dzm[k]));
          b0[k] = (alpha*0.5*(d[k]+fabs(d[k])-d[k+1]+fabs(d[k+1]))+0.5*(1-alpha)*(d[k]-d[k+1]))/(0.5*(2.0*dzm[k]));
          c0[k] = -(alpha*0.5*(d[k+1]+fabs(d[k+1])) + 0.5*(1-alpha)*d[k+1])/(0.5*(2.0*dzm[k]));			
        }
      }

      // add on explicit diffusion to RHS (utmp)
      if(grid->Nke[j]-grid->etop[j]>1) { // more than one vertical layer on edge
        // Explicit part of the viscous term over wetted parts of edge
        // for the interior cells
        // for the new vertical coordinate, the bottom layer should be the layer zfb>bufferheight since interior dzf can be zero.
        for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++){
            //phys->utmp[j][k]+=dt*(1-theta)*(a[k]*phys->u[j][k-1]-(a[k]+b[k])*phys->u[j][k]+b[k]*phys->u[j][k+1]);
            phys->utmp[j][k]+=dt*(fac2)*(a[k]*phys->u[j][k-1]-(a[k]+b[k])*phys->u[j][k]+b[k]*phys->u[j][k+1]) + 
              dt*(fac3)*(a[k]*phys->u_old[j][k-1]-(a[k]+b[k])*phys->u_old[j][k]+b[k]*phys->u_old[j][k+1]);
        }

        // Top cell
        // account for no slip conditions which are assumed if CdT = -1 
        if(phys->CdT[j] == -1){ // no slip on top
          // phys->utmp[j][grid->etop[j]]+=
          //   dt*(1-theta)*(a[grid->etop[j]]*-phys->u[j][grid->etop[j]]-
          //     (a[grid->etop[j]]+b[grid->etop[j]])*phys->u[j][grid->etop[j]]+
          //     b[grid->etop[j]]*phys->u[j][grid->etop[j]+1]);
          phys->utmp[j][grid->etop[j]]+=
            dt*(fac2)*(a[grid->etop[j]]*-phys->u[j][grid->etop[j]]-
              (a[grid->etop[j]]+b[grid->etop[j]])*phys->u[j][grid->etop[j]]+
              b[grid->etop[j]]*phys->u[j][grid->etop[j]+1]) + dt*(fac3)*(a[grid->etop[j]]*-phys->u_old[j][grid->etop[j]]-
              (a[grid->etop[j]]+b[grid->etop[j]])*phys->u_old[j][grid->etop[j]]+
              b[grid->etop[j]]*phys->u_old[j][grid->etop[j]+1]);
        }
        else{ // standard drag law code
          // phys->utmp[j][grid->etop[j]]+=dt*(1-theta)*(-(b[grid->etop[j]]+2.0*phys->CdT[j]*
          //       fabs(phys->utmp[j][grid->etop[j]])/
          //       (2.0*dzm[grid->etop[j]]))*
          //     phys->u[j][grid->etop[j]]
          //     +b[grid->etop[j]]*phys->u[j][grid->etop[j]+1]);
          phys->utmp[j][grid->etop[j]]+=dt*(fac2)*(-(b[grid->etop[j]]+2.0*phys->CdT[j]*
                fabs(phys->utmp[j][grid->etop[j]])/
                (2.0*dzm[grid->etop[j]]))*
              phys->u[j][grid->etop[j]]+b[grid->etop[j]]*phys->u[j][grid->etop[j]+1]) + dt*(fac3)*(-(b[grid->etop[j]]+2.0*phys->CdT[j]*
                fabs(phys->u_old[j][grid->etop[j]])/
                (2.0*dzm[grid->etop[j]]))*
              phys->u_old[j][grid->etop[j]]
              +b[grid->etop[j]]*phys->u_old[j][grid->etop[j]+1]); //not sure about the utmp->uold here 
        }

        // Bottom cell
        // account for no slip conditions which are assumed if CdB = -1
        if(phys->CdB[j] == -1){ // no slip on bottom
         // some sort of strange error here... in previous code, now fixed 
         // phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(phys->CdB[j]+phys->CdT[j])/
         //(2.0*dzm[grid->etop[j]])*
        //fabs(phys->u[j][grid->etop[j]])*phys->u[j][grid->etop[j]];
           phys->utmp[j][grid->Nke[j]-1]+=dt*(1-theta)*(a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2]-
               (a[grid->Nke[j]-1]+b[grid->Nke[j]-1])*phys->u[j][grid->Nke[j]-1]+
               b[grid->Nke[j]-1]*-phys->u[j][grid->Nke[j]-1]);
          //phys->utmp[j][grid->Nke[j]-1]+=erf((prop->rtime)/2830)*(dt*(fac2)*(a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2]-
	  //   (a[grid->Nke[j]-1]+b[grid->Nke[j]-1])*phys->u[j][grid->Nke[j]-1]+
	  //   b[grid->Nke[j]-1]*-phys->u[j][grid->Nke[j]-1]) + dt*(fac3)*(a[grid->Nke[j]-1]*phys->u_old[j][grid->Nke[j]-2]-
	  //   (a[grid->Nke[j]-1]+b[grid->Nke[j]-1])*phys->u_old[j][grid->Nke[j]-1]+
	  //								  b[grid->Nke[j]-1]*-phys->u_old[j][grid->Nke[j]-1]));

        } else{ 
          // standard drag law code
          if(prop->vertcoord==1)
            if(!prop->subgrid)
              // phys->utmp[j][grid->Nke[j]-1]+=dt*(1-theta)*(
              //     a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2] -
              //     (a[grid->Nke[j]-1] 
              //      + 2.0*phys->CdB[j]*fabs(phys->utmp[j][grid->Nke[j]-1])/
              //      (2.0*dzm[grid->Nke[j]-1]))*
              //     phys->u[j][grid->Nke[j]-1]);
              phys->utmp[j][grid->Nke[j]-1]+=dt*(fac2)*(
                  a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2] -
                  (a[grid->Nke[j]-1] 
                   + 2.0*phys->CdB[j]*fabs(phys->utmp[j][grid->Nke[j]-1])/
                   (2.0*dzm[grid->Nke[j]-1]))*
                  phys->u[j][grid->Nke[j]-1]) + dt*(fac3)*(
                  a[grid->Nke[j]-1]*phys->u_old[j][grid->Nke[j]-2] -
                  (a[grid->Nke[j]-1] 
                   + 2.0*phys->CdB[j]*fabs(phys->u_old[j][grid->Nke[j]-1])/
                   (2.0*dzm[grid->Nke[j]-1]))*
                  phys->u_old[j][grid->Nke[j]-1]);
            else 
              // phys->utmp[j][grid->Nke[j]-1]+=dt*(1-theta)*(
              //     a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2] -
              //     (a[grid->Nke[j]-1] 
              //      + phys->CdB[j]*fabs(phys->utmp[j][grid->Nke[j]-1])/
              //     subgrid->dzboteff[j])*phys->u[j][grid->Nke[j]-1]);
              phys->utmp[j][grid->Nke[j]-1]+=dt*(fac2)*(
                  a[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2] -
                  (a[grid->Nke[j]-1] 
                   + phys->CdB[j]*fabs(phys->utmp[j][grid->Nke[j]-1])/
                  subgrid->dzboteff[j])*phys->u[j][grid->Nke[j]-1]) + dt*(fac3)*(
                  a[grid->Nke[j]-1]*phys->u_old[j][grid->Nke[j]-2] -
                  (a[grid->Nke[j]-1] 
                   + phys->CdB[j]*fabs(phys->u_old[j][grid->Nke[j]-1])/
                  subgrid->dzboteff[j])*phys->u_old[j][grid->Nke[j]-1]);
           else
           {
              if(!prop->subgrid)
                // apply drag forcing at the layer where zfb>buffer height
                // phys->utmp[j][Nkeb]+=dt*(1-theta)*(
                //   a[Nkeb]*phys->u[j][Nkeb-1] -
                //   (a[Nkeb]+2.0*phys->CdB[j]*fabs(phys->utmp[j][Nkeb])/
                //    (2.0*dzm[Nkeb]))*phys->u[j][Nkeb]);
                phys->utmp[j][Nkeb]+=dt*(fac2)*(
                  a[Nkeb]*phys->u[j][Nkeb-1] -
                  (a[Nkeb]+2.0*phys->CdB[j]*fabs(phys->utmp[j][Nkeb])/
                   (2.0*dzm[Nkeb]))*phys->u[j][Nkeb]) + dt*(fac3)*(
                  a[Nkeb]*phys->u_old[j][Nkeb-1] -
                  (a[Nkeb]+2.0*phys->CdB[j]*fabs(phys->u_old[j][Nkeb])/
                   (2.0*dzm[Nkeb]))*phys->u_old[j][Nkeb]);
              else
                // phys->utmp[j][Nkeb]+=dt*(1-theta)*(
                //   a[Nkeb]*phys->u[j][Nkeb-1] -
                //   (a[Nkeb]+phys->CdB[j]*fabs(phys->utmp[j][Nkeb])/
                //    subgrid->dzboteff[j])*phys->u[j][Nkeb]);     
                phys->utmp[j][Nkeb]+=dt*(fac2)*(
                  a[Nkeb]*phys->u[j][Nkeb-1] -
                  (a[Nkeb]+phys->CdB[j]*fabs(phys->utmp[j][Nkeb])/
                   subgrid->dzboteff[j])*phys->u[j][Nkeb]) + dt*(fac3)*(
                  a[Nkeb]*phys->u_old[j][Nkeb-1] -
                  (a[Nkeb]+phys->CdB[j]*fabs(phys->u_old[j][Nkeb])/
                   subgrid->dzboteff[j])*phys->u_old[j][Nkeb]);              
             // for other cell give a 100 drag coefficient
             for(k=Nkeb+1;k<grid->Nke[j];k++)
                // phys->utmp[j][k]+=-dt*(1-theta)*2.0*100*fabs(phys->utmp[j][k])/
                //    (2.0*dzm[k])*phys->u[j][k];
                phys->utmp[j][k]+=-dt*(fac2)*2.0*100*fabs(phys->utmp[j][k])/
                   (2.0*dzm[k])*phys->u[j][k] -dt*(fac3)*2.0*100*fabs(phys->u_old[j][k])/
                   (2.0*dzm[k])*phys->u_old[j][k];
           }
        }
      } else {  // one layer for edge
        // drag on bottom boundary
        if(phys->CdB[j] == -1){ // no slip on bottom
          // phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(
          //     2.0*(2.0*(prop->nu + c[k]))*phys->u[j][grid->etop[j]]/
          //     ((2.0*dzm[grid->etop[j]])*
          //      (2.0*dzm[grid->etop[j]])));
          phys->utmp[j][grid->etop[j]]-=2.0*dt*(fac2)*(
              2.0*(2.0*(prop->nu + c[k]))*phys->u[j][grid->etop[j]]/
              ((2.0*dzm[grid->etop[j]])*
               (2.0*dzm[grid->etop[j]]))) + 2.0*dt*(fac3)*(
              2.0*(2.0*(prop->nu + c[k]))*phys->u_old[j][grid->etop[j]]/
              ((2.0*dzm[grid->etop[j]])*
               (2.0*dzm[grid->etop[j]])));
        }
        else{ // standard drag law formation on bottom
          // need to change for the new vertical coordinate if there is only one layer
          // same as the original model
          // subgrid part
          if(!prop->subgrid)
            // phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(phys->CdB[j])/
            //   (2.0*dzm[grid->etop[j]])*
            //   fabs(phys->utmp[j][grid->etop[j]])*phys->u[j][grid->etop[j]];
            phys->utmp[j][grid->etop[j]]-=2.0*dt*(fac2)*(phys->CdB[j])/
              (2.0*dzm[grid->etop[j]])*
              fabs(phys->utmp[j][grid->etop[j]])*phys->u[j][grid->etop[j]] + 2.0*dt*(fac3)*(phys->CdB[j])/
              (2.0*dzm[grid->etop[j]])*
              fabs(phys->u_old[j][grid->etop[j]])*phys->u_old[j][grid->etop[j]];
          else
            // phys->utmp[j][grid->etop[j]]-=dt*(1-theta)*(phys->CdB[j])/
            //   subgrid->dzboteff[j]*fabs(phys->utmp[j][grid->etop[j]])*phys->u[j][grid->etop[j]];  
            phys->utmp[j][grid->etop[j]]-=dt*(fac2)*(phys->CdB[j])/
              subgrid->dzboteff[j]*fabs(phys->utmp[j][grid->etop[j]])*phys->u[j][grid->etop[j]]+dt*(fac3)*(phys->CdB[j])/
              subgrid->dzboteff[j]*fabs(phys->u_old[j][grid->etop[j]])*phys->u_old[j][grid->etop[j]];        
        }

        // drag on top boundary
        if(phys->CdT[j] == -1){ // no slip on top
          // phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(
          //     2.0*(2.0*(prop->nu + c[k]))*phys->u[j][grid->etop[j]]/
          //     ((2.0*dzm[grid->etop[j]])*
          //      (2.0*dzm[grid->etop[j]])));
          phys->utmp[j][grid->etop[j]]-=2.0*dt*(fac2)*(
              2.0*(2.0*(prop->nu + c[k]))*phys->u[j][grid->etop[j]]/
              ((2.0*dzm[grid->etop[j]])*
               (2.0*dzm[grid->etop[j]]))) + 2.0*dt*(fac3)*(
              2.0*(2.0*(prop->nu + c[k]))*phys->u_old[j][grid->etop[j]]/
              ((2.0*dzm[grid->etop[j]])*
               (2.0*dzm[grid->etop[j]])));
        }
        else{ // standard drag law formulation on top
          if(!prop->subgrid)
            // phys->utmp[j][grid->etop[j]]-=2.0*dt*(1-theta)*(phys->CdT[j])/
            //   (2.0*dzm[grid->etop[j]])*
            //   fabs(phys->utmp[j][grid->etop[j]])*phys->u[j][grid->etop[j]];
            phys->utmp[j][grid->etop[j]]-=2.0*dt*(fac2)*(phys->CdT[j])/
              (2.0*dzm[grid->etop[j]])*
              fabs(phys->utmp[j][grid->etop[j]])*phys->u[j][grid->etop[j]] + 2.0*dt*(fac3)*(phys->CdT[j])/
              (2.0*dzm[grid->etop[j]])*
              fabs(phys->u_old[j][grid->etop[j]])*phys->u_old[j][grid->etop[j]];
          else
            // phys->utmp[j][grid->etop[j]]-=dt*(1-theta)*(phys->CdT[j])/
            //   subgrid->dzboteff[j]*fabs(phys->utmp[j][grid->etop[j]])*phys->u[j][grid->etop[j]];     
            phys->utmp[j][grid->etop[j]]-=dt*(fac2)*(phys->CdT[j])/
              subgrid->dzboteff[j]*fabs(phys->utmp[j][grid->etop[j]])*phys->u[j][grid->etop[j]] + dt*(fac3)*(phys->CdT[j])/
              subgrid->dzboteff[j]*fabs(phys->u_old[j][grid->etop[j]])*phys->u_old[j][grid->etop[j]];          
        }
      }

      // add marsh explicit term
      if(!prop->subgrid){
        if(prop->marshmodel)
          MarshExplicitTerm(grid,phys,prop,j,theta,dt,myproc);
      }else{
        if(!subgrid->dragpara && prop->marshmodel)
          MarshExplicitTerm(grid,phys,prop,j,theta,dt,myproc);
      }

      // add on explicit vertical momentum advection only if there is more than one vertical layer edge.
      if(prop->vertcoord==1 && prop->nonlinear && prop->thetaM>=0 && grid->Nke[j]-grid->etop[j]>1) {
        for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++)
          phys->utmp[j][k]-=prop->dt*(1-prop->thetaM)*(a0[k]*phys->u[j][k-1]+b0[k]*phys->u[j][k]+c0[k]*phys->u[j][k+1]);
	
      	// Top boundary
        phys->utmp[j][grid->etop[j]]-=prop->dt*(1-prop->thetaM)*((a0[grid->etop[j]]+b0[grid->etop[j]])*phys->u[j][grid->etop[j]]+c0[grid->etop[j]]*phys->u[j][grid->etop[j]+1]);
	
      	// Bottom boundary
        phys->utmp[j][grid->Nke[j]-1]-=prop->dt*(1-prop->thetaM)*(a0[grid->Nke[j]-1]*phys->u[j][grid->Nke[j]-2]+(b0[grid->Nke[j]-1]+c0[grid->Nke[j]-1])*phys->u[j][grid->Nke[j]-1]);
      }

      // add new vertical momentum advection for the general vertical coordinate
      // here phys->u still store u^n 
      if(prop->vertcoord!=1 && prop->nonlinear && prop->thetaM>=0 && grid->Nke[j]-grid->etop[j]>1)
      {
        // conservative form part
        for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++)
          phys->utmp[j][k]-=prop->dt*(a0[k]*(fac2*phys->u_old[j][k-1]+fac3*phys->u_old2[j][k-1])+
                b0[k]*(fac2*phys->u_old[j][k]+fac3*phys->u_old2[j][k])+c0[k]*(fac2*phys->u_old[j][k+1]+fac3*phys->u_old2[j][k+1]));
        // Top boundary
        phys->utmp[j][grid->etop[j]]-=prop->dt*((a0[grid->etop[j]]+b0[grid->etop[j]])*
                (fac2*phys->u_old[j][grid->etop[j]]+fac3*phys->u_old2[j][grid->etop[j]])+
                c0[grid->etop[j]]*(fac2*phys->u_old[j][grid->etop[j]+1]+fac3*phys->u_old2[j][grid->etop[j]+1]));
        
        // Bottom boundary
        phys->utmp[j][grid->Nke[j]-1]-=prop->dt*(a0[grid->Nke[j]-1]*
                (fac2*phys->u_old[j][grid->Nke[j]-2]+fac3*phys->u_old2[j][grid->Nke[j]-2])+
                (b0[grid->Nke[j]-1]+c0[grid->Nke[j]-1])*(fac2*phys->u_old[j][grid->Nke[j]-1]+fac3*phys->u_old2[j][grid->Nke[j]-1]));
      
        // second part udomegadz
        if(prop->wetdry)
          for(k=grid->etop[j];k<grid->Nke[j];k++)
            phys->utmp[j][k]+=prop->dt*(fac2*phys->u[j][k]+fac3*phys->u_old2[j][k])*
            (def2*(vert->omega_old[nc1][k]-vert->omega_old[nc1][k+1])+def1*(vert->omega_old[nc2][k]-vert->omega_old[nc2][k+1]))/
            grid->dg[j]/(0.5*(2.0*dzm[k]));
      }

      // Now set up the coefficients for the tridiagonal inversion for the
      // implicit part.  These are given from the arrays above in the discrete operator
      // d^2U/dz^2 = -theta dt a_k U_{k-1} + (1+theta dt (a_k+b_k)) U_k - theta dt b_k U_{k+1}
      // = RHS of utmp

      // Right hand side U** is given by d[k] here.
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
        e1[k]=1.0;
        d[k]=phys->utmp[j][k];
      }


      if(grid->Nke[j]-grid->etop[j]>1) { // for more than one vertical layer
        // Top cells
        c[grid->etop[j]]=-fac1*dt*b[grid->etop[j]];
        // account for no slip conditions which are assumed if CdT = -1 
        if(phys->CdT[j] == -1){ // no slip
          b[grid->etop[j]]=1.0+fac1*dt*(a[grid->etop[j]]+a[grid->etop[j]+1]+b[grid->etop[j]]);
        } else { // standard drag law
          b[grid->etop[j]]=1.0+fac1*dt*(b[grid->etop[j]]+
              2.0*phys->CdT[j]*fabs(phys->utmp[j][grid->etop[j]])/
              (2.0*dzm[grid->etop[j]]));
        }
        a[grid->etop[j]]=0;     // set a_1=0 (not used in tridiag solve)

        // Bottom cell
        c[grid->Nke[j]-1]=0;   // set c_N=0 (not used in tridiag solve)
        // account for no slip conditions which are assumed if CdB = -1  
        if(phys->CdB[j] == -1){ // no slip
          b[grid->Nke[j]-1]=1.0+fac1*dt*(a[grid->Nke[j]-1]+b[grid->Nke[j]-1]+b[grid->Nke[j]-2]);
        } else { 
          // standard drag law
          if(prop->vertcoord==1)
            if(!prop->subgrid)
              b[grid->Nke[j]-1]=1.0+fac1*dt*(a[grid->Nke[j]-1]+
                  2.0*phys->CdB[j]*fabs(phys->utmp[j][grid->Nke[j]-1])/
                  (2.0*dzm[grid->Nke[j]-1]));
            else
              b[grid->Nke[j]-1]=1.0+fac1*dt*(a[grid->Nke[j]-1]+
                  phys->CdB[j]*fabs(phys->utmp[j][grid->Nke[j]-1])/
                  subgrid->dzboteff[j]);       
          else
          {
            if(!prop->subgrid)
              b[Nkeb]=1.0+fac1*dt*(a[Nkeb]+
                2.0*phys->CdB[j]*fabs(phys->utmp[j][Nkeb])/
                (2.0*dzm[Nkeb]));
            else
              b[Nkeb]=1.0+fac1*dt*(a[Nkeb]+
                phys->CdB[j]*fabs(phys->utmp[j][Nkeb])/
                subgrid->dzboteff[j]);    
            // for the layer below the effective layer (zc>bufferheight) give 100 drag coefficient
            for(k=Nkeb+1;k<grid->Nke[j];k++){
              b[k]=1.0+fac1*dt*2.0*100*fabs(phys->utmp[j][k])/
                (2.0*dzm[k]);
              a[k]=0;
              c[k]=0;
            }
          }
        }
        if(prop->vertcoord==1)
        {
          a[grid->Nke[j]-1]=-fac1*dt*a[grid->Nke[j]-1];
          // Interior cells
          for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++) {
            c[k]=-fac1*dt*b[k];
            b[k]=1.0+fac1*dt*(a[k]+b[k]);
            a[k]=-fac1*dt*a[k];
          }
        } else {
          // defined bottom cell
          a[Nkeb]=-fac1*dt*a[Nkeb];
          c[Nkeb]=0;
          // interior cell
          for(k=grid->etop[j]+1;k<Nkeb;k++) {
            c[k]=-fac1*dt*b[k];
            b[k]=1.0+fac1*dt*(a[k]+b[k]);
            a[k]=-fac1*dt*a[k];
          }            
        }
      } else {

        // for a single vertical layer
        b[grid->etop[j]] = 1.0;

        // account for no slip conditions which are assumed if CdB = -1  
        if(phys->CdB[j] == -1){ // no slip
          b[grid->etop[j]]+=4.0*fac1*dt*2.0*(prop->nu+c[k])/
            ((2.0*dzm[grid->etop[j]])*
             (2.0*dzm[grid->etop[j]]));
        }
        else{
          if(!prop->subgrid)
            b[grid->etop[j]]+=2.0*fac1*dt*fabs(phys->utmp[j][grid->etop[j]])/
              (2.0*dzm[grid->etop[j]])*
              (phys->CdB[j]);
          else
            b[grid->etop[j]]+=fac1*dt*fabs(phys->utmp[j][grid->etop[j]])/
              subgrid->dzboteff[j]*(phys->CdB[j]); 
        }

        // account for no slip conditions which are assumed if CdT = -1 
        if(phys->CdT[j] == -1){
          b[grid->etop[j]]+=4.0*fac1*dt*2.0*(prop->nu+c[k])/
            ((2.0*dzm[grid->etop[j]])*
             (2.0*dzm[grid->etop[j]]));
        }
        else{
          if(!prop->subgrid)
            b[grid->etop[j]]+=2.0*fac1*dt*fabs(phys->utmp[j][grid->etop[j]])/
              (2.0*dzm[grid->etop[j]])*
              (phys->CdT[j]);
          else
            b[grid->etop[j]]+=fac1*dt*fabs(phys->utmp[j][grid->etop[j]])/
              subgrid->dzboteff[j]*(phys->CdT[j]);
        }
      }	  

      // Now add implicit term from marsh friction
      if(!prop->subgrid)
      {
        if(prop->marshmodel)
          for(jv=marsh->marshtop[j];jv<grid->Nke[j];jv++)
            b[jv]+=MarshImplicitTerm(grid,phys,prop,j,jv,theta,dt,myproc);
      } else {
        if(!subgrid->dragpara && prop->marshmodel)
        {
          for(jv=marsh->marshtop[j];jv<grid->Nke[j];jv++)
            b[jv]+=MarshImplicitTerm(grid,phys,prop,j,jv,theta,dt,myproc);
        }
      }

      // Now add on implicit terms for vertical momentum advection, only if there is more than one layer
      if(prop->vertcoord==1 && prop->nonlinear && prop->thetaM>=0 && grid->Nke[j]-grid->etop[j]>1) {
        for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++) {
          a[k]+=prop->dt*prop->thetaM*a0[k];
          b[k]+=prop->dt*prop->thetaM*b0[k];
          c[k]+=prop->dt*prop->thetaM*c0[k];
        }
        // Top boundary
        b[grid->etop[j]]+=prop->dt*prop->thetaM*(a0[grid->etop[j]]+b0[grid->etop[j]]);
        c[grid->etop[j]]+=prop->dt*prop->thetaM*c0[grid->etop[j]];
        // Bottom boundary 
        a[grid->Nke[j]-1]+=prop->dt*prop->thetaM*a0[grid->Nke[j]-1];
        b[grid->Nke[j]-1]+=prop->dt*prop->thetaM*(b0[grid->Nke[j]-1]+c0[grid->Nke[j]-1]);
      }

      // add vertical momentum advection in the new general vertical coordinate
      if(prop->vertcoord!=1 && prop->nonlinear && grid->Nke[j]-grid->etop[j]>1 && prop->thetaM>=0) {
        // conservative part d(\omega u)dz
        for(k=grid->etop[j]+1;k<grid->Nke[j]-1;k++) {
          a[k]+=prop->dt*fac1*a0[k];
          b[k]+=prop->dt*fac1*b0[k];
          c[k]+=prop->dt*fac1*c0[k];
        }
        // Top boundary
        b[grid->etop[j]]+=prop->dt*fac1*(a0[grid->etop[j]]+b0[grid->etop[j]]);
        c[grid->etop[j]]+=prop->dt*fac1*c0[grid->etop[j]];
        // Bottom boundary 
        a[grid->Nke[j]-1]+=prop->dt*fac1*a0[grid->Nke[j]-1];
        b[grid->Nke[j]-1]+=prop->dt*fac1*(b0[grid->Nke[j]-1]+c0[grid->Nke[j]-1]);
        // second part -udwdz
        if(prop->wetdry)
          for(k=grid->etop[j];k<grid->Nke[j];k++)
            b[k]-=prop->dt*fac1*
            (def2*(vert->omega_old[nc1][k]-vert->omega_old[nc1][k+1])+def1*(vert->omega_old[nc2][k]-vert->omega_old[nc2][k+1]))/grid->dg[j]/dzm[k];
      }

      // implicit method for u/JdJdt term
      // fully implicit
      if(prop->vertcoord!=1 && prop->nonlinear && !prop->wetdry)
        if(vert->dJdtmeth==0)
        {
          def1 = grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
          def2 = grid->def[nc2*grid->maxfaces+grid->gradf[2*j+1]];
          dgf = def1+def2;
          for(k=grid->etop[j];k<grid->Nke[j];k++) 
          {
            phys->utmp[j][k]-=(0*phys->u_old[j][k]+0*phys->u_old2[j][k])*
             (def2/dgf*(1-grid->dzzold[nc1][k]/grid->dzz[nc1][k])+def1/dgf*(1-grid->dzzold[nc2][k]/grid->dzz[nc2][k]));
           
            b[k]+=1*(def2/dgf*(1-grid->dzzold[nc1][k]/grid->dzz[nc1][k])+def1/dgf*(1-grid->dzzold[nc2][k]/grid->dzz[nc2][k])); 
          }           
        }  



      for(k=grid->etop[j];k<grid->Nke[j];k++) {

        if(dzm[k]==0) {
          printf("Exiting because j %d dzz[%d][%d]=%f or dzz[%d][%d]=%f dv1 %e dv2 %e nk1 %d nk2 %d nke %d\n",j,
              nc1,k,grid->dzz[nc1][k],nc2,k,grid->dzz[nc2][k],grid->dv[nc1],grid->dv[nc2],grid->Nk[nc1],grid->Nk[nc2],grid->Nke[j]);
          exit(0);
        }

        if(a[k]!=a[k]) printf("a[%d] problems, dzz[%d][%d]=%f\n",k,j,k,grid->dzz[j][k]);
        
        if(b[k]!=b[k] || b[k]==0){
          if(prop->subgrid)
            printf("proc %d n %d ne %d b[%d] problems, b=%f dzf %e nke %d etop %d Nk %d %d dv %e %e hmin %e %e Cd %e\n",myproc,prop->n,j, k,b[k],
              grid->dzf[j][k],grid->Nke[j],grid->etop[j],grid->Nk[grid->grad[2*j]],grid->Nk[grid->grad[2*j+1]],grid->dv[grid->grad[2*j]],grid->dv[grid->grad[2*j+1]],
              subgrid->hmin[grid->grad[2*j]],subgrid->hmin[grid->grad[2*j+1]],phys->CdB[j]);
          else
            printf("proc %d n %d ne %d b[%d] problems, b=%f dzf %e nke %d etop %d Nk %d %d dv %e %e Cd %e\n",myproc,prop->n,j, k,b[k],
              grid->dzf[j][k],grid->Nke[j],grid->etop[j],grid->Nk[grid->grad[2*j]],grid->Nk[grid->grad[2*j+1]],grid->dv[grid->grad[2*j]],grid->dv[grid->grad[2*j+1]]
              ,phys->CdB[j]);                
        }

        if(c[k]!=c[k]) printf("c[%d] problems\n",k);
      }

      // Now utmp will have U*** in it, which is given by A^{-1}U**, and E will have
      // A^{-1}e1, where e1 = [1,1,1,1,1,...,1]^T 
      // Store the tridiagonals so they can be used twice (TriSolve alters the values
      // of the elements in the diagonals!!!
      for(k=0;k<grid->Nke[j];k++) {
        a0[k]=a[k];
        b0[k]=b[k];
        c0[k]=c[k];
      }

      if(grid->Nke[j]-grid->etop[j]>1) { // more than one layer (z level)
        TriSolve(&(a[grid->etop[j]]),&(b[grid->etop[j]]),&(c[grid->etop[j]]),
            &(d[grid->etop[j]]),&(phys->utmp[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);
        TriSolve(&(a0[grid->etop[j]]),&(b0[grid->etop[j]]),&(c0[grid->etop[j]]),
            &(e1[grid->etop[j]]),&(E[j][grid->etop[j]]),grid->Nke[j]-grid->etop[j]);	
      } else {  // one layer (z level)
        phys->utmp[j][grid->etop[j]]/=b[grid->etop[j]];
        E[j][grid->etop[j]]=1.0/b[grid->etop[j]];
      }

      // Now vertically integrate E to create the vertically integrated flux-face
      // values that comprise the coefficients of the free-surface solver.  This
      // will create the D vector, where D=DZ^T E (which should be given by the
      // depth when there is no viscosity.
      phys->D[j]=0;
      for(k=grid->etop[j];k<grid->Nke[j];k++) 
      {  
        phys->D[j]+=E[j][k]*grid->dzf[j][k];
      } 
    }
  }

  j=2;
  /*
  for(k=0;k<3;k++) {
    //    printf("k=%d, dzf=%.2e, utmp=%.2e, (a,b,c)=%.2e,%.2e,%.2e\n",
    //	   k,grid->dzf[j][k],phys->utmp[j][k],a0[k],b0[k],c0[k]);
  }
  */
  theta=theta0;

  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  for(j=0;j<grid->Ne;j++) 
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      if(phys->utmp[j][k]!=phys->utmp[j][k]) {
        if(prop->subgrid)
          printf("n %d Error in function Predictor at j=%d k=%d Nke %d etop %d (U***=nan) cd %e nc1 %d nc2 %d V %e %e \n",prop->n,j,k,grid->Nke[j],grid->etop[j],phys->CdB[j],grid->grad[2*j],grid->grad[2*j+1],\
                  subgrid->Veff[grid->grad[2*j]],subgrid->Veff[grid->grad[2*j+1]]);
        else
          printf("n %d Error in function Predictor at j=%d k=%d Nke %d etop %d (U***=nan) cd %e nc1 %d nc2 %d dzf=%.2e mark=%d\n",prop->n,j,k,grid->Nke[j],grid->etop[j],phys->CdB[j],grid->grad[2*j],grid->grad[2*j+1],grid->dzf[j][k],grid->mark[j]);
        exit(1);
      }

  // So far we have U*** and D.  Now we need to create h* in htmp.   This
  // will comprise the source term for the free-surface solver.  Before we
  // do this we need to set the new velocity at the open boundary faces and
  // place them into utmp.  
  BoundaryVelocities(grid,phys,prop,myproc,comm);
  OpenBoundaryFluxes(NULL,phys->utmp,NULL,grid,phys,prop);

  for(j=0;j<grid->Ne;j++) 
    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      if(phys->utmp[j][k]!=phys->utmp[j][k]) {
        printf("n %d Error in function Predictor at j=%d k=%d (U***=nan) cd %e nc1 %d nc2 %d\n",prop->n,j,k,phys->CdB[j],grid->grad[2*j],grid->grad[2*j+1]);
        exit(1);
      }

  // for computational cells
  // now phys->u=u^n phys->u_old2=u^n-1 phys->utmp=u*
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    sum = 0;

    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      nc1=grid->grad[2*ne];
      nc2=grid->grad[2*ne+1];
      if(nc1==-1)
        nc1=nc2;
      if(nc2==-1)
        nc2=nc1;
      for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
      {  
        sum+=(fac2*phys->u_old[ne][k]+fac1*phys->utmp[ne][k]+fac3*phys->u_old2[ne][k])*
          grid->df[ne]*normal*grid->dzf[ne][k];  
      }
    }
    if(prop->subgrid)
    {
      phys->htmp[i]=subgrid->Veff[i]-dt*sum;
      StoreSubgridOldAceffandVeff(grid, myproc);
      subgrid->rhs[i]=phys->htmp[i];
    }
    else {
      phys->htmp[i]=grid->Ac[i]*phys->h[i]-dt*sum
	+grid->Ac[i]*(phys->zB[i]-phys->zBold[i]);
      if(phys->htmp[i]!=phys->htmp[i]){
	printf("n %d something wrong on the source term of h at cell %d, sum=%f\n",prop->n,i,sum);
	exit(1);
      }
    }
  }

  //printf("made it to about halfway \n");
 

  for(i=0;i<grid->Nc;i++){
    // store the old h as h_n-1 in the next time step
    phys->h_old[i]=phys->h[i];
    phys->dhdt[i]=phys->h[i];
  }
  
  sum=0;
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    sum+=phys->htmp[i];
  }

  // test whether the total mass of the model is positive
  if(sum<0 && prop->subgrid)
  {
    printf("something wrong in b>=0 n=%d sum=%e\n",prop->n,sum);
      MPI_Finalize();
      exit(EXIT_FAILURE);
  }

  // culvert model part
  if(prop->culvertmodel)
    CulvertInitIteration(grid,phys, prop, 1,myproc);

  // Now we have the required components for the CG solver for the free-surface:
  //
  // h^{n+1} - g*(theta*dt)^2/Ac * Sum_{faces} D_{face} dh^{n+1}/dn df N = htmp
  //
  // L(h) = b
  //
  // L(h) = h + 1/Ac * Sum_{faces} D_{face} dh^{n+1}/dn N
  // b = htmp
  //
  // As the initial guess let h^{n+1} = h^n, so just leave it as it is to
  // begin the solver.

  if(!prop->culvertmodel)
  {
    nf=0;
    min=INFTY;
    while(1){
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];
        if(prop->subgrid)
          subgrid->residual[i]=-subgrid->Veff[i]+subgrid->Aceff[i]*phys->h[i];
      }

      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];
        if(prop->subgrid){
          phys->htmp[i]=subgrid->rhs[i];
          phys->htmp[i]-=(subgrid->Veff[i]-phys->h[i]*subgrid->Aceff[i]); 
        }
        
        if(nf==0 && prop->subgrid){
          if(subgrid->rhs[i]!=subgrid->rhs[i])
            printf("%d right hand side wrong=%f\n",i,subgrid->rhs[i]);
          if(subgrid->Veff[i]!=subgrid->Veff[i])
            printf("%d Veff wrong=%f\n",i,subgrid->Veff[i]);
          if(subgrid->Aceff[i]!=subgrid->Aceff[i])
            printf("%d Aceff wrong=%f\n",i,subgrid->Aceff[i]);
          if(phys->htmp[i]!=phys->htmp[i])
            printf("%d htmp wrong=%f\n",i,phys->htmp[i]);
        }
      }
      
      // store h for each iteration
      if(prop->subgrid)
        for(i=0;i<grid->Nc;i++)
          subgrid->hiter[i]=phys->h[i];

      //      printf("HERE\n");
      //      BiCGSolve(grid,phys,prop,myproc,numprocs,comm);
      CGSolve(grid,phys,prop,myproc,numprocs,comm); 
      
      // for original suntans	
      if(!prop->subgrid)
        break;
      
      // subgrid part
      if(prop->subgrid)
        UpdateSubgridVeff(grid, phys, prop, myproc);

      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];

        if(prop->subgrid) {
          subgrid->residual[i]+=subgrid->Veff[i]-subgrid->Aceff[i]*phys->h[i];
        }
      }     

      ISendRecvCellData2D(subgrid->residual,grid,myproc,comm);
      sum=InnerProduct(subgrid->residual,subgrid->residual,grid,myproc,numprocs,comm);

      if(nf==0){
        if(sum>1)
          sum0=sum;
        else
          sum0=1;
      }

      if(prop->subgrid)
        UpdateSubgridAceff(grid, phys, prop, myproc);

      if(sqrt(sum)<subgrid->eps)
        break;

      if(sqrt(sum/sum0)<subgrid->eps)
        break;

      nf++;
      if(min>sqrt(sum/sum0)){
        if(prop->subgrid)
          for(i=0;i<grid->Nc;i++)
            subgrid->hiter_min[i]=subgrid->hiter[i];
        min=sqrt(sum/sum0);        
      }

      if(nf>10) { 
        if(fabs(sqrt(sum/sum0)-min)<0.001)
          break;

        /*if(nf==11)
        {
          for(i=0;i<grid->Nc;i++)
            phys->h[i]=subgrid->hiter_min[i];
          ISendRecvCellData2D(phys->h,grid,myproc,comm);
          UpdateSubgridVeff(grid, phys, prop, myproc);
          UpdateSubgridAceff(grid, phys, prop, myproc);
        }*/

        if(nf>50)
        {
           printf("n %d nf %d something maybe wrong for convergence at time step min %e sum %e sum0 %e r %e\n",prop->n,nf,min,sum,sum0,sqrt(sum/sum0));
           printf("iteration for subgrid is more than 50 times. stop program\n");
           exit(1);
        }
      }
    }
  } else {

    // outer loop
    UpdateCulvertQcoef(grid,prop,0, myproc);
    nf=0;
    min2=INFTY;
    while(1)
    {
      for(i=0;i<grid->Nc;i++){
        culvert->pressure3[i]=phys->h[i];           
      }    
      min=INFTY;
      nf1=0;

      while(1){
        // inner loop
        if(prop->subgrid)
        {
          UpdateSubgridVeff(grid, phys, prop, myproc);
          UpdateSubgridAceff(grid, phys, prop, myproc);
        }
        
        // initialize each iteration 
        CulvertInitIteration(grid,phys,prop,-1,myproc);
        // calculate the source term for free surface eqns. based on casulli's method
        CulvertIterationSource(grid,phys,prop,theta,dt,myproc);      
        // CG solver for free surface
        CGSolve(grid,phys,prop,myproc,numprocs,comm);
 
        if(prop->subgrid)
          UpdateSubgridVeff(grid, phys, prop, myproc);

        // calculate residual
        CheckCulvertCondition(grid,phys,prop,myproc);

        if(prop->subgrid)
          UpdateSubgridAceff(grid, phys, prop, myproc);

        // calculate norm error
        ISendRecvCellData2D(culvert->condition,grid,myproc,comm);

        // the total residual
        culvert->sum=InnerProduct(culvert->condition,culvert->condition,grid,myproc,numprocs,comm);

        // infinity norm
        if(nf1==0){
          if(culvert->sum>1)
            sum0=culvert->sum;
          else
            sum0=1;
        }

        if(sqrt(culvert->sum)<culvert->eps)
          break;

        // check culvert condition
        if(sqrt(culvert->sum/sum0)<culvert->eps)
          break;
      
        if(min>sqrt(culvert->sum/sum0))
          min=sqrt(culvert->sum/sum0);
      
        if(nf1>10)
        {
          printf("proc %d 1 n %d nf %d something maybe wrong for convergence at time step min %e sum %e sum0 %e r %e\n",myproc,prop->n,nf1 ,
            min,culvert->sum,sum0,sqrt(culvert->sum/sum0));

          if(fabs(sqrt(culvert->sum/sum0)-min)<1e-3)
            break;

          if(nf1>50){
            printf("iteration for subgrid is more than 50 times. stop program\n");
            exit(1);
          }
        }

        nf1++;
      }
  
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];           
        culvert->condition2[i]=culvert->Qcoef[i]*(culvert->top[i]-phys->h[i]);
        culvert->pressure2[i]=phys->h[i];
      } 

      UpdateCulvertQcoef(grid,prop,0,myproc);

      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];           
        culvert->condition2[i]+=culvert->Qcoef[i]*(-culvert->top[i]+phys->h[i]);
      } 

      // calculate norm error
      ISendRecvCellData2D(culvert->condition2,grid,myproc,comm);

      // the total residual
      culvert->sum=InnerProduct(culvert->condition2,culvert->condition2,grid,myproc,numprocs,comm);
      
      if(nf==0){
        if(culvert->sum>1)
          sum=culvert->sum;
        else
          sum=1;
      }
      nf++;

      if(sqrt(culvert->sum)<culvert->eps)
        break;

      // check culvert condition
      if(sqrt(culvert->sum/sum)<culvert->eps)
        break; 

      if(min2>sqrt(culvert->sum/sum)){
        min2=sqrt(culvert->sum/sum);
        for(i=0;i<grid->Nc;i++)
          culvert->pressure4[i]=culvert->pressure3[i];
      }

      if(nf>10) { 
        if(nf>50){
          exit(1);
        }
        //exit(EXIT_FAILURE);
        if(fabs(sqrt(culvert->sum/sum)-min2)<0.001)
          break;

        if(nf==20)
        {
          printf("myproc %d 2 n %d nf %d something maybe wrong for convergence at time step min %e sum %e sum0 %e r %e\n",myproc,prop->n,nf,min2,culvert->sum,sum,sqrt(culvert->sum/sum));
          for(i=0;i<grid->Nc;i++)
          {
            phys->h[i]=culvert->pressure[i];
            culvert->pressure2[i]=phys->h[i];
          }
          ISendRecvCellData2D(phys->h,grid,myproc,comm);
          ISendRecvCellData2D(culvert->pressure2,grid,myproc,comm);
          UpdateCulvertQcoef(grid,prop,0,myproc);
        }
        if(fabs(sqrt(culvert->sum/sum)-min2)<0.001)
          break;        
      }
    }
  }

  // Add back the implicit barotropic term to obtain the 
  // hydrostatic horizontal velocity field.
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->u[j][k]=phys->utmp[j][k]-prop->grav*fac1*dt*E[j][k]*
        (phys->h[nc1]-phys->h[nc2])/grid->dg[j];
    }

    // set dry cells (with zero height) to have zero velocity
    // CHANGE
    if(grid->etop[j]==grid->Nke[j]-1 && grid->dzz[nc1][grid->etop[j]]<=1*DRYCELLHEIGHT &&
        grid->dzz[nc2][grid->etop[j]]<=1*DRYCELLHEIGHT)
      phys->u[j][grid->etop[j]]=0;
  }
  // STOP dzm HERE

  // correct cells drying below DRYCELLHEIGHT above the 
  // bathymetry
  // make sure Culvert height is much bigger than DRYCELLHEIGHT when culvertmodel==1 
  for(i=0;i<grid->Nc;i++){
    // for a cell with 0 dzf at all edges reset to the old phys->h[i]
    // It will not affect the results since dzf is all zero. And any neighbouring cell will have no effects     
    if(prop->subgrid && prop->wetdry)
    {
      flag=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        if(grid->etop[ne]<(grid->Nke[ne]-1) || grid->dzf[ne][grid->Nke[ne]-1]>0)
          flag=1;
      }
      if(!flag)
        phys->h[i]=phys->h_old[i]; 
    }

    if(phys->h[i]<=(-grid->dv[i]+DRYCELLHEIGHT)) {
      phys->hcorr[i]=-grid->dv[i]+DRYCELLHEIGHT-phys->h[i];
      phys->h[i]=-grid->dv[i]+DRYCELLHEIGHT;
      phys->active[i]=0;
      //phys->s[i][grid->Nk[i]-1]=0;
      //phys->T[i][grid->Nk[i]-1]=0;
      //if(prop->computeSediments && prop->n>1+prop->nstart)
        //for(k=0;k<sediments->Nsize;k++)
          //sediments->SediC[k][i][grid->Nk[i]-1]=0;

    } else {
      phys->hcorr[i]=0;
      phys->active[i]=1;
    }

    //if(phys->h[i]<=(-grid->dv[i]+1e-3))
      //phys->active[i]=0;
  }

  int neigh;
  flag=0;
  REAL Ac, htmp, u_im;
  for(i=0;i<grid->Nc;i++) {
    sum=0;
    Ac=grid->Ac[i];
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      for(nf=0;nf<grid->nfaces[i];nf++) {               
	ne = grid->face[i*grid->maxfaces+nf];
	neigh = grid->neigh[i*grid->maxfaces+nf];
	if(neigh==-1)
	  neigh=i;
	
	if(k<grid->Nke[ne]) {
	  u_im = fac1*phys->u[ne][k] + fac2*phys->u_old[ne][k] + fac3*phys->u_old2[ne][k];
	  sum+=prop->dt*u_im*grid->dzf[ne][k]*grid->normal[i*grid->maxfaces+nf]*grid->df[ne]/Ac;
	}
      }
    }

    flag=0;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      if(grid->mark[grid->face[i*grid->maxfaces+nf]]==2 || 
	 grid->mark[grid->face[i*grid->maxfaces+nf]]==3 ||
	 grid->mark[grid->face[i*grid->maxfaces+nf]]==6) {
	flag=1;
	break;
      }
    }
    
    htmp=phys->h_old[i]-sum;
    if(!flag && fabs((phys->h[i]-htmp)/(phys->h_old[i]+grid->dv[i]))>1e-10) {
      //      printf("h at location (xv,yv)=%.2f,%.2f after CGSOLVE differs by %.3e\n",
      //	     grid->xv[i],grid->yv[i],fabs((htmp-phys->h[i])/(phys->h_old[i]+grid->dv[i])));
      break;
    }
  }
  //  if(flag)
  //    printf("Proc %d: After CGSOLVE depth differs from value inferred from continuity by more than 1e-10.\n",myproc);

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    phys->dhdt[i]=(phys->h[i]-phys->dhdt[i])/dt;
  }

  // added culvert part
  if(prop->culvertmodel){
    StoreCulvertPressure(phys->h, grid->Nc, 1, myproc);
  }

  if(prop->subgrid)
  {
    UpdateSubgridVeff(grid, phys, prop, myproc);
    UpdateSubgridAceff(grid, phys, prop, myproc);
  }

  // add back residual in free surface solver to ensure restrict mass conservation
  if(prop->wetdry && prop->subgrid){
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      sum0=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        normal = grid->normal[i*grid->maxfaces+nf];
        for(k=grid->etop[ne];k<grid->Nke[ne];k++){ 
          sum0+=prop->dt*(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k])*grid->df[ne]*normal*grid->dzf[ne][k];
        }
      }
      subgrid->Verr[i]=subgrid->Veff[i]-subgrid->Veffold[i]+sum0;
      subgrid->Veff[i]=subgrid->Veffold[i]-sum0;   
    }
  }
   
  if(prop->subgrid)
    ISendRecvCellData2D(subgrid->Veff,grid,myproc,comm);

  if(prop->subgrid)
    for(i=0;i<grid->Nc;i++){
      if(subgrid->Veff[i]<=(DRYCELLHEIGHT*grid->Ac[i]))
      {
        phys->h[i]=-grid->dv[i]+DRYCELLHEIGHT;
        phys->active[i]=0;
        subgrid->Veff[i]=(DRYCELLHEIGHT*grid->Ac[i]);
        //phys->s[i][grid->Nk[i]-1]=0;
        //if(prop->computeSediments && prop->n>1+prop->nstart)
          //for(k=0;k<sediments->Nsize;k++)
            //sediments->SediC[k][i][grid->Nk[i]-1]=0;
      }
      if(subgrid->Veff[i]<1e-3*grid->Ac[i])
        phys->active[i]=0;
    }
  
  // recalculate new free surface based on the modified volume
  if(prop->subgrid && prop->wetdry){
    UpdateSubgridFreeSurface(grid,phys,prop,myproc);
    ISendRecvCellData2D(phys->h,grid,myproc,comm);  
    if(prop->culvertmodel){
      ISendRecvCellData2D(culvert->pressure,grid,myproc,comm);  
      ISendRecvCellData2D(culvert->pressure2,grid,myproc,comm);  
    }
    UpdateSubgridVeff(grid, phys, prop, myproc);
    UpdateSubgridAceff(grid, phys, prop, myproc);
    UpdateSubgridHeff(grid, phys, prop, myproc);
  }

  // Use the new free surface to add the implicit part of the free-surface
  // pressure gradient to the horizontal momentum.
  //
  // This was removed because it violates the assumption of linearity in that
  // the discretization only knows about grid cells below grid->ctopold[].
  /*
     for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
     j = grid->edgep[jptr];

     nc1 = grid->grad[2*j];
     nc2 = grid->grad[2*j+1];

     if(grid->etop[j]>grid->etopold[j]) 
     for(k=0;k<grid->etop[j];k++)
     phys->u[j][k]=0;
     else 
     for(k=grid->etop[j];k<grid->etopold[j];k++)
     phys->u[j][k]=phys->utmp[j][k]-prop->grav*theta*dt*
     (phys->h[nc1]-phys->h[nc2])/grid->dg[j];
     }
     */

  // Set the flux values at the open boundary (marker=2).  These
  // were set to utmp previously in OpenBoundaryFluxes.
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      phys->u[j][k] = phys->utmp[j][k];
  }

  // Now set the fluxes at the free-surface boundary by assuming dw/dz = 0
  // this is for type 3 boundary conditions
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      if(grid->mark[ne]==3) {
        for(k=grid->etop[ne];k<grid->Nke[ne];k++) {
          phys->u[ne][k] = 0;
          sum=0;
          for(nf1=0;nf1<grid->nfaces[i];nf1++)
            sum+=phys->u[grid->face[i*grid->maxfaces+nf1]][k]*grid->df[grid->face[i*grid->maxfaces+nf1]]*grid->normal[i*grid->maxfaces+nf1];
          phys->u[ne][k]=-sum/grid->df[ne]/grid->normal[i*grid->maxfaces+nf];
        }
      }
    } 
  }

  if(prop->vertcoord!=1 && prop->vertcoord!=5)
    if(vert->modifydzf)
    {
      VerifyFluxHeight(grid,prop,phys,myproc);
      UpdateCellcenteredFreeSurface(grid,prop,phys,myproc);
      ISendRecvCellData2D(phys->h,grid,myproc,comm);  
    }

  // Now update the vertical grid spacing with the new free surface.
  // can comment this out to linearize the free surface 
  if(prop->vertcoord==1)
    UpdateDZ(grid,phys,prop, 0); 
  else {
    /*
    // Need to compute SfHp and SfHm for use in UpdateLayerThickness()
    if(prop->vertcoord==2 || prop->vertcoord==4)
      TvdFluxHeight(grid, phys, prop, vert->dzfmeth,comm, myproc);
    */
    SetFluxHeight(grid,phys,prop,vert->dzfmeth,comm,myproc);
    UpdateLayerThickness(grid, prop, phys, 0,myproc, numprocs, comm);

    ISendRecvCellData3D(grid->dzz,grid,myproc,comm);

    // OBF 012623 Fix dzz so that sum(dzz)=depth
    /*
    for(i=0;i<grid->Nc;i++){
        sum=0;
        for(k=grid->ctop[i];k<grid->Nk[i];k++)
          sum+=grid->dzz[i][k];
	for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	  grid->dzz[i][k]*=(phys->h[i]+grid->dv[i]-phys->zB[i])/sum;
	}
    }
    */

    /* Update top cell height to ensure h+d = sum(dzz)
    for(i=0;i<grid->Nc;i++){
      sum=0;
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	sum+=grid->dzz[i][k];
      
      grid->dzz[i][grid->ctop[i]]=grid->dzz[i][grid->ctop[i]]+(phys->h[i]+grid->dv[i]-phys->zB[i])-sum;
    }
    */

    /* Print an error of CG solver produced depth is not consistent with depth from continuity
    for(i=0;i<grid->Nc;i++){
      flag=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        if(grid->mark[grid->face[i*grid->maxfaces+nf]]==2 || 
	   grid->mark[grid->face[i*grid->maxfaces+nf]]==3) {
          flag=1;
          break;
        }
      }

      if(!flag) {
        sum=0;
        for(k=0;k<grid->Nk[i];k++)
          sum+=grid->dzz[i][k];
        if(fabs(sum-(phys->h[i]+grid->dv[i]-phys->zB[i]))>1e-6)
          printf("Time step %d, i=%d, error. |Sum(dzz)-H|=%.3e\n",
		 prop->n, i,fabs(sum-(phys->h[i]+grid->dv[i]))/(phys->h[i]+grid->dv[i]));
      }
    }
    */
    
    // OBF 042723 Added damping in small layers
    /*
    REAL Cd=0.1;
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	if(grid->dzf[j][k]==0) {
	  printf("DZF is ZERO AFTER!\n");
	  exit(1);
	}
	
	if(grid->dzf[j][k]<0.01*grid->dz[k]) 
	  phys->u[j][k]/=(1.0+prop->dt*Cd*fabs(phys->utmp[j][k])/grid->dzf[j][k]);
      }
    }
    */
    
    // compute the new zc
    //MPI_Barrier(comm);
    ComputeZc(grid,prop,phys,1,myproc);
  }
  
  // update vertical ac for scalar transport
  if(prop->subgrid)
    UpdateSubgridVerticalAceff(grid, phys, prop, 0, myproc);
}

// For debugging - never used to solve for free surface but is a 2d version of the 3d
// version that is used to solve for the nonhydrostatic pressure when the cross terms are
// present.
static void BiCGSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {
  int i, iptr, n, niters;
  REAL *x, *r, *p, *z, *r0, *v, *t, *s, eps, eps0, epsW = 10E-10, sum;
  REAL rho0 = 1, rho, alpha = 1, omg = 1, beta, tmp;
  int m0=0, n0=0;
  
  z = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  v = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  s = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  t = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  r0 = (REAL *)SunMalloc(grid->Nc*sizeof(REAL), "BiCGSolveN");
  /*
  x = wave->N[m0][n0];
  r = wave->Ntmp[m0][n0];
  p = wave->Nold[m0][n0];
  */
  x = phys->h;
  r = phys->hold;
  p = phys->htmp;

  for(i = 0; i< grid->Nc; i++)
    r[i]=0;
  
    niters = prop->maxiters;
  //niters = 1000;

  for(i = 0; i< grid->Nc; i++){
    z[i] = 0;
    v[i] = 0;
    s[i] = 0;
    t[i] = 0;
    r0[i] = 0;

  }

  HCoefficients(phys->hcoef,phys->hfcoef,grid,phys,prop);  
  
  ISendRecvCellData2D(x, grid, myproc, comm);
  OperatorH(x,z,phys->hcoef,phys->hfcoef,grid,phys,prop);  
  //  OperatorN(m0, n0, x, z, grid, prop);

  for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
    i = grid->cellp[iptr];
    p[i] = p[i]-z[i]; //Now initial P stores the initial residual.
    r[i] = p[i];  //Initial residual
    r0[i] = r[i]; //Another initial residual for BiCG
  }

  for(iptr=grid->celldist[1]; iptr<grid->celldist[2]; iptr++){
    i = grid->cellp[iptr];
    p[i] = 0;
  }
  eps0=eps=InnerProduct(r,r,grid,myproc,numprocs,comm);
  
  if (eps > epsW){
    if(!prop->resnorm) eps0 = 1;

    for(n = 0; n < niters && eps!=0; n++){ 
      rho = InnerProduct(r0,r,grid,myproc,numprocs,comm);
      beta = rho/rho0*alpha/omg; //When n = 0, rho0, alpha, omg = 1.
 
     
      for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
	i = grid->cellp[iptr];
	p[i] = r[i]+beta*(p[i]-omg*v[i]); 
      }


      ISendRecvCellData2D(p, grid, myproc, comm);
      OperatorH(p,v,phys->hcoef,phys->hfcoef,grid,phys,prop);        
      //      OperatorN(m0, n0, p, v, grid, prop);
      tmp = InnerProduct(r0,v,grid,myproc,numprocs,comm);
      tmp = 1/tmp;   
      alpha = rho*tmp;

    
     
      for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
      	i = grid->cellp[iptr];
      	s[i] = r[i]-alpha*v[i];
      }
      
      eps = InnerProduct(s, s, grid, myproc, numprocs, comm);
      
      if(VERBOSE>3) printf("BiCGSolve Iteration (1): %d, resid=%e, proc=%d\n",n,eps,myproc);
      if(eps<SMALL){
	x[i] += alpha*p[i];
	break;
      }
    
      ISendRecvCellData2D(s, grid, myproc, comm);
      OperatorH(s,t,phys->hcoef,phys->hfcoef,grid,phys,prop);              
      //      OperatorN(m0, n0, s, t, grid, prop);
      
      tmp = InnerProduct(t, t, grid, myproc, numprocs, comm);
      tmp = 1/tmp;
      omg = InnerProduct(t, s, grid, myproc, numprocs, comm)*tmp;
    
      
      for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
      	i = grid->cellp[iptr];
      	x[i] += alpha*p[i] + omg*s[i];
	r[i] = s[i]-omg*t[i];
      }
      rho0 = rho;
      
      eps = InnerProduct(r, r, grid, myproc, numprocs, comm);
      
      if(VERBOSE>3) printf("BiCGSolve Iteration (2): %d, resid=%e, proc=%d\n",n,eps,myproc);
      if(eps < epsW)
	break;
    }
  }

  if(myproc==0 && VERBOSE>3) {
    if(eps0 < epsW) {
      printf("Step %d, norm of action density (%d, %d) source at = %e is already small\n", prop->n, m0, n0, eps0);
    } else {
      if(n==niters) printf("Warning... Step %d, action density (%d, %d) iteration not converging after %d steps! RES=%e > %.2e\n",
			   prop->n, m0, n0, n, eps, SMALL);
      else printf("Step %d, BiCGSolve action density (%d, %d) converged after %d iterations, rsdl = %e < %e\n",
		  prop->n, m0, n0, n, eps, epsW);
    }
  }
  ISendRecvCellData2D(x, grid, myproc, comm);

  SunFree(z, grid->Nc*sizeof(REAL), "BiCGSolveN");
  SunFree(v, grid->Nc*sizeof(REAL), "BiCGSolveN");
  SunFree(s, grid->Nc*sizeof(REAL), "BiCGSolveN");
  SunFree(t, grid->Nc*sizeof(REAL), "BiCGSolveN");
  SunFree(r0, grid->Nc*sizeof(REAL), "BiCGSolveN");
}

/*
 * Function: CGSolve
 * Usage: CGSolve(grid,phys,prop,myproc,numprocs,comm);
 * ----------------------------------------------------
 * Solve the free surface equation using the conjugate gradient algorithm.
 *
 * The source term upon entry is in phys->htmp, which is placed into p, and
 * the free surface upon entry is in phys->h, which is placed into x.
 *
 */
static void CGSolve(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {

  int i, iptr, n, niters;
  REAL *x, *r, *rtmp, *p, *z, mu, nu, eps, eps0, alpha, alpha0;

  x = phys->h;
  r = phys->hold;
  rtmp = phys->htmp2;
  z = phys->htmp3;
  p = phys->htmp;

  niters = prop->maxiters;

  // Create the coefficients for the operator
  if(!prop->culvertmodel){
    HCoefficients(phys->hcoef,phys->hfcoef,grid,phys,prop);
  }else{
    CulvertHCoefficients(phys->hcoef,phys->hfcoef,grid,phys,prop,myproc);
  }

  // For the boundary term (marker of type 3):
  // 1) Need to set x to zero in the interior points, but
  //    leave it as is for the boundary points.
  // 2) Then set z=Ax and substract b = b-z so that
  //    the new problem is Ax=b with the boundary values
  //    on the right hand side acting as forcing terms.
  // 3) After b=b-z for the interior points, then need to
  //    set b=0 for the boundary points.

  /* Fix to account for boundary cells (type 3) in h */
  // 1) x=0 interior cells 
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    x[i]=0;
  }
  ISendRecvCellData2D(x,grid,myproc,comm);
  OperatorH(x,z,phys->hcoef,phys->hfcoef,grid,phys,prop);

  // 2) b = b-z
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    p[i] = p[i] - z[i];    
    r[i] = p[i];
    x[i] = 0;
  }    
  // 3) b=0 for the boundary cells
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) { 
    i = grid->cellp[iptr]; 

    p[i] = 0; 
  }     

  // continue with CG as expected now that boundaries are handled
  if(prop->hprecond==1) {
    HPreconditioner(r,rtmp,grid,phys,prop);
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      p[i] = rtmp[i];
    }
    alpha = alpha0 = InnerProduct(r,rtmp,grid,myproc,numprocs,comm);
  } else {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      p[i] = r[i];
    }
    alpha = alpha0 = InnerProduct(r,r,grid,myproc,numprocs,comm);
  }
  if(!prop->resnorm) alpha0 = 1;

  if(prop->hprecond==1)
    eps=eps0=InnerProduct(r,r,grid,myproc,numprocs,comm);
  else
    eps=eps0=alpha0;

  // Iterate until residual is less than prop->epsilon
  for(n=0;n<niters && eps!=0 && alpha!=0;n++) {

    ISendRecvCellData2D(p,grid,myproc,comm);
    OperatorH(p,z,phys->hcoef,phys->hfcoef,grid,phys,prop);

    mu = 1/alpha;
    nu = alpha/InnerProduct(p,z,grid,myproc,numprocs,comm);

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      x[i] += nu*p[i];
      r[i] -= nu*z[i];
    }
    if(prop->hprecond==1) {
      HPreconditioner(r,rtmp,grid,phys,prop);
      alpha = InnerProduct(r,rtmp,grid,myproc,numprocs,comm);
      mu*=alpha;
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];

        p[i] = rtmp[i] + mu*p[i];
      }
    } else {
      alpha = InnerProduct(r,r,grid,myproc,numprocs,comm);
      mu*=alpha;
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];

        p[i] = r[i] + mu*p[i];
      }
    }

    if(prop->hprecond==1)
      eps=InnerProduct(r,r,grid,myproc,numprocs,comm);
    else
      eps=alpha;

    if(VERBOSE>3 && myproc==0) printf("CGSolve free-surface Iteration: %d, resid=%e\n",n,sqrt(eps/eps0));
    if(sqrt(eps/eps0)<prop->epsilon) 
      break;
  }
  if(myproc==0 && VERBOSE>2){
    if(eps==0){
      printf("Warning...Time step %d, norm of free-surface source is 0.\n",prop->n);
    } else {
      if(n==niters)  printf("Warning... Time step %d, Free-surface iteration not converging after %d steps! RES=%e > %.2e\n",
          prop->n,n,sqrt(eps/eps0),prop->qepsilon);
      else printf("Time step %d, CGSolve free-surface converged after %d iterations, res=%e < %.2e\n",
          prop->n,n,sqrt(eps/eps0),prop->epsilon);
    }
  }
  // Send the solution to the neighboring processors
  ISendRecvCellData2D(x,grid,myproc,comm);
}


/*
 * Function: HPreconditioner
 * Usage: HPreconditioner(r,rtmp,grid,phys,prop);
 * ----------------------------------------------
 * Multiply the vector x by the inverse of the preconditioner M with
 * xc = M^{-1} x
 *
 */
static void HPreconditioner(REAL *x, REAL *y, gridT *grid, physT *phys, propT *prop) {
  int i, iptr;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    y[i]=x[i]/phys->hcoef[i];
  }
}

/*
 * Function: HCoefficients
 * Usage: HCoefficients(coef,fcoef,grid,phys,prop);
 * --------------------------------------------------
 * Compute coefficients for the free-surface solver.  fcoef stores
 * coefficients at the flux faces while coef stores coefficients
 * at the cell center.  If L is the linear operator on x, then
 *
 * L(x(i)) = coef(i)*x(i) - sum(m=1:3) fcoef(ne)*x(neigh)
 * coef(i) = (Ac(i) + sum(m=1:3) tmp*D(ne)*df(ne)/dg(ne))
 * fcoef(ne) = tmp*D(ne)*df(ne)/dg(ne)
 *
 * where tmp = prop->grav*(theta*dt)^2
 *
 */
static void HCoefficients(REAL *coef, REAL *fcoef, gridT *grid, physT *phys, propT *prop) {

  int i, j, iptr, jptr, ne, nf, check;
  REAL tmp, h0, boundary_flag,fac;

  fac=prop->imfac1;

  tmp=prop->grav*pow(fac*prop->dt,2);

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    check=1;
    i = grid->cellp[iptr];
    if(!prop->subgrid)
      coef[i] = grid->Ac[i];
    else
      coef[i] = subgrid->Aceff[i];

    for(nf=0;nf<grid->nfaces[i];nf++) 
      if(grid->neigh[i*grid->maxfaces+nf]!=-1) {
        ne = grid->face[i*grid->maxfaces+nf];
        fcoef[i*grid->maxfaces+nf]=tmp*phys->D[ne]*grid->df[ne]/grid->dg[ne];
        coef[i]+=fcoef[i*grid->maxfaces+nf];

        if(fcoef[i*grid->maxfaces+nf]>0)
          check=0;
      } 
 
    // when 0=0 exists make sure hnew=hold
    if(check) //&& prop->subgrid)
    {
      coef[i]=1.0;
      phys->htmp[i]=phys->h[i];
    }
  }

}

/*
 *
 * Function: InnerProduct
 * Usage: InnerProduct(x,y,grid,myproc,numprocs,comm);
 * ---------------------------------------------------
 * Compute the inner product of two one-dimensional arrays x and y.
 * Used for the CG method to solve for the free surface.
 *
 */
static REAL InnerProduct(REAL *x, REAL *y, gridT *grid, int myproc, int numprocs, MPI_Comm comm) {

  int i, iptr;
  REAL sum, mysum=0;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    mysum+=x[i]*y[i];
  }
  MPI_Reduce(&mysum,&(sum),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Bcast(&sum,1,MPI_DOUBLE,0,comm);

  return sum;
}  

/*
 * Usage: OperatorH(x,y,grid,phys,prop);
 * -------------------------------------
 * Given a vector x, computes the left hand side of the free surface 
 * Poisson equation and places it into y with y = L(x), where
 *
 * L(x(i)) = coef(i)*x(i) + sum(m=1:3) fcoef(ne)*x(neigh)
 * coef(i) = (Ac(i) + sum(m=1:3) tmp*D(ne)*df(ne)/dg(ne))
 * fcoef(ne) = tmp*D(ne)*df(ne)/dg(ne)
 *
 * where tmp = prop->grav*(theta*dt)^2
 *
 */
static void OperatorH(REAL *x, REAL *y, REAL *coef, REAL *fcoef, gridT *grid, physT *phys, propT *prop) {

  int i, j, iptr, jptr, ne, nf;
  REAL tmp = prop->grav*pow(prop->theta*prop->dt,2), h0, boundary_flag;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    y[i] = coef[i]*x[i];
    for(nf=0;nf<grid->nfaces[i];nf++) 
      if(grid->neigh[i*grid->maxfaces+nf]!=-1)
        y[i]-=fcoef[i*grid->maxfaces+nf]*x[grid->neigh[i*grid->maxfaces+nf]];
  }

}
