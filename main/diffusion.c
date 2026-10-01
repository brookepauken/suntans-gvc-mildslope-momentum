/*
 * File: diffusion.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Compute eddy-viscosity via Smagorinsky or for central-differencing as dictated by
 * the lax-wendroff scheme.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "diffusion.h"
#include "util.h"

/*
 * Function: LaxWendroff
 * Usage: LaxWendroff(grid,phys,prop,myproc,comm);
 * -----------------------------------------------
 * Compute the numerical diffusion coefficients required to stablize momentum advection
 * when central-differencing is used (nonlinear=2).
 *
 */
void LaxWendroff(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int i, j, k, nc1, nc2;
  REAL numax=-INFTY, numin=INFTY, numinall, numaxall, laxVertical;

  laxVertical=1;
  if(grid->Nkmax==1)
    laxVertical=0;

  if(prop->nonlinear==2) {
    for(i=0;i<grid->Nc;i++) 
      for(k=0;k<grid->Nk[i];k++) {
	phys->nu_lax[i][k]=
	  0.5*(pow(laxVertical*0.5*(phys->w_old[i][k]+phys->w_old[i][k+1]),2)+
	       pow(phys->uc[i][k],2)+pow(phys->vc[i][k],2))*prop->dt;
	if(phys->nu_lax[i][k]<numin) numin=phys->nu_lax[i][k];
	if(phys->nu_lax[i][k]>numax) numax=phys->nu_lax[i][k];
      }

    if(VERBOSE>2) {
      MPI_Reduce(&numin,&numinall,1,MPI_DOUBLE,MPI_MAX,0,comm);
      MPI_Reduce(&numax,&numaxall,1,MPI_DOUBLE,MPI_MAX,0,comm);
      
      if(myproc==0) printf("Lax-Wendroff diffusion coefficients: numin = %.3e, numax = %.3e\n",numin,numax);
    }
  }
}

void Smagorinsky(gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc) {
  int i, k, iptr, nf, ne, normal, neigh;
  REAL *dudx, *dudy, *dvdx, *dvdy, umag, Cs = 0.1, n1, n2, df, S;
  REAL numax=-INFTY, numin=INFTY, numinall, numaxall, fu, u0, delta_u;
  
  dudx = phys->a;
  dudy = phys->b;
  dvdx = phys->c;
  dvdy = phys->d;

  // Need to eliminate nu_lax when u<u0
  u0 = 2.5;
  delta_u = 0.25;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i=grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      dudx[k]=dudy[k]=dvdx[k]=dvdy[k]=0;
    }

    for(nf=0;nf<grid->nfaces[i];nf++) {
      neigh=grid->neigh[i*grid->maxfaces+nf];
      normal=grid->normal[i*grid->maxfaces+nf];
      ne=grid->face[i*grid->maxfaces+nf];
      n1=grid->n1[ne];
      n2=grid->n2[ne];
      df=grid->df[ne];

      if(neigh==-1) {
	neigh=i;
      }
      for(k=grid->etop[ne];k<grid->Nke[ne];k++) {
	dudx[k]+=0.5*(phys->uc[i][k]+phys->uc[neigh][k])*n1*normal*df;
	dudy[k]+=0.5*(phys->uc[i][k]+phys->uc[neigh][k])*n2*normal*df;
	dvdx[k]+=0.5*(phys->vc[i][k]+phys->vc[neigh][k])*n1*normal*df;
	dvdy[k]+=0.5*(phys->vc[i][k]+phys->vc[neigh][k])*n2*normal*df;
      }
    }
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      dudx[k]/=grid->Ac[i];
      dvdx[k]/=grid->Ac[i];
      dudy[k]/=grid->Ac[i];
      dvdy[k]/=grid->Ac[i];

      S = sqrt(2*dudx[k]*dudx[k]+2*dvdy[k]*dvdy[k]+pow(dudy[k]+dvdx[k],2.0));

      // Use nu_lax when |u|>u0
      umag = sqrt(pow(phys->uc[i][k],2.0)+pow(phys->vc[i][k],2.0));
      fu = 1;0.5*(1+tanh((umag-u0)/delta_u));

      phys->nu_lax[i][k]=Cs*Cs*grid->Ac[i]*S*fu;
      
      if(phys->nu_lax[i][k]<numin) numin=phys->nu_lax[i][k];
      if(phys->nu_lax[i][k]>numax) numax=phys->nu_lax[i][k];
    }
  }

  if(VERBOSE>1) {
    MPI_Reduce(&numin,&numinall,1,MPI_DOUBLE,MPI_MIN,0,comm);
    MPI_Reduce(&numax,&numaxall,1,MPI_DOUBLE,MPI_MAX,0,comm);
    
    if(myproc==0 && !(prop->n%100)) {
      printf("Smagorinsky diffusion coefficients: numin = %.3e, numax = %.3e\n",
	     numinall,numaxall);
    }
  }
}
