/*
 * File: scalars.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * ----------------------------------------
 * This file contains the scalar transport function.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "scalars.h"
#include "util.h"
#include "tvd.h"
#include "subgrid.h"
#include "initialization.h"
#include "check.h"

#define SMALL_CONSISTENCY 1e-5

REAL smin_value, smax_value;

// Local function
//static void ReducePentadiag(REAL *a, REAL *b, REAL *c, REAL *d, REAL *e, REAL *r, REAL *u, int N);
static void GaussSeidel(REAL *a, REAL *b, REAL *c, REAL *d, REAL *e, REAL *r, REAL *u, REAL **u_init, int i, REAL tol, int maxiters, int N);

/*
 * Function: UpdateScalars
 * Usage: UpdateScalars(grid,phys,prop,w_im,scalar,Cn,kappa,kappaH,kappa_tv,theta);
 * --------------------------------------------------------------------------------
 * Update the scalar quantity stored in the array denoted by scal using the
 * theta method for vertical advection and vertical diffusion and Adams-Bashforth
 * for horizontal advection and diffusion.
 *
 * Cn must store the AB terms from time step n-1 for this scalar
 * kappa denotes the vertical scalar diffusivity
 * kappaH denotes the horizontal scalar diffusivity
 * kappa_tv denotes the vertical turbulent scalar diffusivity
 *
 */
void UpdateScalars(gridT *grid, physT *phys, propT *prop, REAL **w_im, REAL **SfH_tm1, REAL **SfH_tm2, REAL **SfHv_t, REAL **SfHv_tm1,  REAL **scal, REAL **scal_old, REAL **boundary_scal, int BCtop, int BCbot, REAL **Cn,
		   REAL kappa, REAL kappaH, REAL **kappa_tv, REAL theta,
		   REAL **src1, REAL **src2, REAL *Ftop, REAL *Fbot, int alpha_top, int alpha_bot,
		   MPI_Comm comm, int myproc, int checkflag, int timestepping, int TVDscheme) 
{
  int i, iptr, j, jptr, ib, k, nf, ktop;
  int Nc=grid->Nc, normal, nc1, nc2, ne;
  REAL df, dg, Ac, dt=prop->dt, fab, *a, *b, *c, *d, *e, *r, *ap, *am, *bd, dznew, mass, *sp, *temp;
  REAL smin, smax, div_local, div_local_max=0, div_da, div_da_max=0, alpha,sum,sum1,sum2;
  int k1, k2, kmin, imin, kmax, imax, mincount, maxcount, allmincount, allmaxcount, flag;
  int  div_local_count, div_da_count;
  REAL fab1,fab2,fab3,fac1,fac2,fac3; // implicit and explicit scheme factors
  prop->TVD = TVDscheme;
  REAL *pent_a, *pent_b, *pent_c, *pent_d, *pent_e; //coefficients of pentdiagonal for QUICK
  REAL temp_k, z_cutoff, num_bot_layers, num_mid_layers, bot_cutoff, sumz, cent, span, kappaH_k;
  REAL **sum_neighs;
  //current

  // These are used mostly debugging to turn on/off vertical and horizontal TVD.
  prop->horiTVD = 1;
  prop->vertTVD = 1;

  ap = phys->ap;
  am = phys->am;
  bd = phys->bp;
  temp = phys->bm;
  a = phys->a;
  b = phys->b;
  c = phys->c;
  d = phys->d;
  e = phys->e;

  r = phys->pent_source; //source term 

  pent_a = phys->pent_a;
  pent_b = phys->pent_b;
  pent_c = phys->pent_c;
  pent_d = phys->pent_d;
  pent_e = phys->pent_e;

  sum_neighs = phys->sum_neighs;

  ktop=0;

  //printf("about to enter fab loop \n");
  // Never use AB2 //?
  if(1) {
    fab=1;
    for(i=0;i<grid->Nc;i++)
      for(k=0;k<grid->Nk[i];k++)
        Cn[i][k]=0;
  } else
      fab=1.5;
  
  
    if(!prop->readOldVelocity){
      if(prop->n==1) {
      fab1=1;
      fab2=fab3=0;
      if(prop->TVD>=5){
        //at first timestep, calculate transverse curvature term for SHARP/QUICK for vertical face interopolation 
        SumNeighborScalars(grid, phys, scal, sum_neighs,comm, myproc);
      }
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
  
    if(prop->n==1 || prop->n==prop->nstart+1) {
      if(prop->TVD>=5){
        //at first timestep, calculate transverse curvature term for SHARP/QUICK for vertical face interopolation 
        SumNeighborScalars(grid, phys, scal, sum_neighs,comm, myproc);
      }
    }
  }



  // store the old value
  // here stmp=scal^n scal_old=scal^n-1
  for(i=0;i<Nc;i++){
    for(k=0;k<grid->Nk[i];k++) { 
      phys->stmp[i][k]=scal[i][k];
    }
    if(prop->n==prop->nstart+1 && prop->TVD>=5){ //these are both w_im, put to local version?
      //at first time step, interpolate intialized scalar onto vertical faces. Save in SfHv_t (Scalar face half vertical at time t)
        GetScalarOnFace(SfHv_t, phys->wp, phys->wm, w_im, grid->dzz,phys->stmp,i,grid->Nk[i],ktop,prop->dt,prop->TVD, BCtop, BCbot, sum_neighs);
        if(prop->readOldVelocity){
          GetScalarOnFace(SfHv_tm1, phys->wp, phys->wm, w_im, grid->dzz,scal_old,i,grid->Nk[i],ktop,prop->dt,prop->TVD, BCtop, BCbot, sum_neighs);
        } else {
        for(k=0;k<grid->Nk[i]+1;k++){
          SfHv_tm1[i][k] = SfHv_t[i][k];
        }
        }
    }
  }

  //printf("through scalar init \n");

  //theta method
  //fac1=prop->thetaS;
  //fac2=1-prop->thetaS;
  //fac3=0;

  //implicit
  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;
  
  // Add on boundary fluxes, using stmp2 as the temporary storage
  // variable

  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      phys->stmp2[i][k]=0; 
  }

  if(boundary_scal) {
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr];
      ib = grid->grad[2*j];
      // Set the value of stmp2 adjacent to the boundary to the value of the boundary.
      // This will be used to add the boundary flux when stmp2 is used again below.
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++)
        phys->stmp2[ib][k]=boundary_scal[jptr-grid->edgedist[2]][k];
    }
  }

  // Compute the scalar on the vertical faces (for horiz. advection)
  // store in SfHp (scalar on face +) anf SfH- (scalar on face -)
  if(prop->TVD && prop->horiTVD){
    HorizontalFaceScalars(grid,phys,prop,scal,phys->u,boundary_scal,prop->TVD,3,comm,myproc);
  }

  //printf("through HorizontalFaceScalars \n");

  //central loop. loop over each cell with pointer i
  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr]; //cell index
    Ac = grid->Ac[i]; //cell area

    if(grid->ctop[i]>=grid->ctopold[i]) { 
      ktop=grid->ctop[i]; 
      dznew=grid->dzz[i][ktop];
    } else {
      ktop=grid->ctopold[i];
      dznew=0;
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++) 
        dznew+=grid->dzz[i][k];      
    }

    // These are the advective components of the tridiagonal
    // at the new time step.
    if(!(prop->TVD && prop->vertTVD))
      for(k=0;k<grid->Nk[i]+1;k++) 
      {
        //add subgrid part
        if(prop->subgrid)
          alpha=subgrid->Acveff[i][k]/Ac;
        else
          alpha=1.0;

        ap[k] = alpha*0.5*(w_im[i][k]+fabs(w_im[i][k]));
        am[k] = alpha*0.5*(w_im[i][k]-fabs(w_im[i][k])); 
      }
    else if (prop->TVD >= 5) {  //for QUICK and SHARP, need a pentadiagonal instead of a tridaigonal
      GetPentDiag(pent_a, pent_b, pent_c, pent_d, pent_e, phys->wp,phys->wm,
		  w_im,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD,BCtop,BCbot);

      for(k=0;k<grid->Nk[i];k++)
	{ //set of pentadiagonal coefficients for implicit part of vertical advection 
	  a[k] = fac1*prop->dt*pent_a[k];
	  b[k] = fac1*prop->dt*pent_b[k];
	  c[k] = fac1*prop->dt*pent_c[k] + grid->dzz[i][k];
	  d[k] = fac1*prop->dt*pent_d[k];
	  e[k] = fac1*prop->dt*pent_e[k];
	}

      c[0] += -grid->dzz[i][0] + dznew; //change dzz to be new dz at top 


    } else  // Compute the ap/am for TVD schemes
    { //calculate coefficients on scalar for vertical advection
      GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
	      w_im,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD,BCtop,BCbot);
      // Vertical TVD is always upwind prop->TVD);
      for(k=0;k<grid->Nk[i]+1;k++) 
      {
        //add subgrid part (do not think this is currently in use)
        if(prop->subgrid)
          alpha=subgrid->Acveff[i][k]/Ac;
        else
          alpha=1.0;

        ap[k]= ap[k]*alpha;
        am[k]= am[k]*alpha;
      }
    
    
      for(k=ktop+1;k<grid->Nk[i];k++) 
	{
	  // add subgrid part
	  if(prop->subgrid)
	    alpha=subgrid->Acceff[i][k]/Ac;
	  else
	    alpha=1.0;
	  b[k-ktop]=fac1*dt*am[k];
	  c[k-ktop]=alpha*grid->dzz[i][k]+fac1*dt*(ap[k]-am[k+1]);
	  d[k-ktop]=-fac1*dt*ap[k+1];

	  a[k-ktop] = 0;
	  e[k-ktop] = 0;
	}
      //top boundary condition
	b[0]=0;
	d[0]=-fac1*dt*ap[ktop+1];
	if(BCtop==1){ //no flux
	  //saying phi_k-1 = phi_k, so add on b to c
	  //need to define am, ap for ktop
	  c[0]=dznew+fac1*dt*(ap[ktop]-am[ktop+1]+am[ktop]);
	} else if(BCtop==2){ //no slip
	  //saying phi_k-1 = -phi_k, so subtract b from c                                            
    c[0]=dznew+fac1*dt*(ap[ktop]-am[ktop+1]-am[ktop]);
	}
      
	if(BCbot==1){
	  // Bottom cell no-flux boundary condition for advection
	  c[(grid->Nk[i]-1)-ktop]+=d[(grid->Nk[i]-1)-ktop];
	} else if(BCbot==2){
	  // Bottom cell no-slip boundary condition for advection
	  c[(grid->Nk[i]-1)-ktop]-=d[(grid->Nk[i]-1)-ktop];
	}
     }

    // Implicit vertical diffusion terms
    if(kappa_tv) {
      for(k=ktop+1;k<grid->Nk[i];k++)
	{
	  if(prop->subgrid)
	    alpha=subgrid->Acveff[i][k]/Ac;
	  else
	    alpha=1.0;
	  
	  bd[k]=alpha*(2.0*kappa+kappa_tv[i][k-1]+kappa_tv[i][k])/
	    (grid->dzz[i][k-1]+grid->dzz[i][k]);
	  
	}
    } else {
      for(k=ktop+1;k<grid->Nk[i];k++)
	{
	  if(prop->subgrid)
	    alpha=subgrid->Acveff[i][k]/Ac;
	  else
	    alpha=1.0;
	  
	  bd[k]=2.0*alpha*kappa/(grid->dzz[i][k-1]+grid->dzz[i][k]);
	}
    }

    for(k=ktop+1;k<grid->Nk[i]-1;k++) 
    {
      b[k-ktop]-=fac1*dt*bd[k];
      c[k-ktop]+=fac1*dt*(bd[k]+bd[k+1]);
      d[k-ktop]-=fac1*dt*bd[k+1];
    }

    if(src1)
      for(k=ktop;k<grid->Nk[i];k++)
      {
        if(prop->subgrid)
          alpha=subgrid->Acceff[i][k]/Ac;
        else
          alpha=1.0;
        c[k-ktop]+=alpha*src1[i][k]*fac1*dt*grid->dzz[i][k];
      }

    // Diffusive fluxes only when more than 1 layer
    if(ktop<grid->Nk[i]-1) {
      // Top cell diffusion
      c[0]+=fac1*dt*(bd[ktop+1]+2*alpha_top*bd[ktop+1]);
      d[0]-=fac1*dt*bd[ktop+1];

      // Bottom cell diffusion
      b[(grid->Nk[i]-1)-ktop]-=fac1*dt*bd[grid->Nk[i]-1];
      c[(grid->Nk[i]-1)-ktop]+=fac1*dt*(bd[grid->Nk[i]-1]+2*alpha_bot*bd[grid->Nk[i]-1]);
    }
  
    // Explicit part into source term d[] 
    for(k=ktop+1;k<grid->Nk[i];k++) 
    {
      if(prop->subgrid)
        alpha=subgrid->Acceffold[i][k]/Ac;
      else 
        alpha=1.0;

      r[k-ktop]=alpha*grid->dzzold[i][k]*phys->stmp[i][k];
    }

    if(src1)
    {
      for(k=ktop+1;k<grid->Nk[i];k++) 
      {
        if(prop->subgrid)
          alpha=subgrid->Acceff[i][k]/Ac;
        else 
          alpha=1.0;
        r[k-ktop]-=alpha*src1[i][k]*dt*grid->dzz[i][k]*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k]);
      }
    }
    r[0]=0;
    if(grid->ctopold[i]<=grid->ctop[i]) {
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
      {
        if(prop->subgrid)
          alpha=subgrid->Acceffold[i][k]/Ac;
        else
          alpha=1.0;
        r[0]+=alpha*grid->dzzold[i][k]*phys->stmp[i][k];
      }

      if(src1)
        for(k=grid->ctopold[i];k<=grid->ctop[i];k++)
        {
          if(prop->subgrid)
            alpha=subgrid->Acceff[i][k]/Ac;
          else
            alpha=1.0;

          r[0]-=alpha*src1[i][k]*dt*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])*grid->dzz[i][k];
        }
    } else {
      if(prop->subgrid)
        alpha=subgrid->Acceffold[i][ktop]/Ac;
      else
        alpha=1.0;      
      //same but top
      r[0]=alpha*grid->dzzold[i][ktop]*phys->stmp[i][ktop];
      if(src1){
        if(prop->subgrid)
          alpha=subgrid->Acceff[i][ktop]/Ac;
        else
          alpha=1.0;          
        r[0]-=alpha*grid->dzz[i][ktop]*src1[i][ktop]*dt*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k]);
      }
    }


    if(!(prop->TVD>=5) || prop->vertTVD==0){
    // Explicit advection and diffusion for TVD schemes
      for(k=ktop+1;k<grid->Nk[i]-1;k++){
	//uses coefficients from current time step to interpolate old values onto top and bottom faces
	r[k-ktop]-=dt*(am[k]*(fac2*phys->stmp[i][k-1]+fac3*scal_old[i][k-1])+
		       (ap[k]-am[k+1])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])-
		       ap[k+1]*(fac2*phys->stmp[i][k+1]+fac3*scal_old[i][k+1]))-
	  dt*(bd[k]*(fac2*phys->stmp[i][k-1]+fac3*scal_old[i][k-1])-
	      (bd[k]+bd[k+1])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])+
	       bd[k+1]*(fac2*phys->stmp[i][k+1]+fac3*scal_old[i][k+1]));
      }
      
    } else{ //for QUICK and SHARP
    //add previous timestep n terms for advection and all diffusion to source term
      for(k=ktop+1;k<grid->Nk[i]-1;k++){
        r[k-ktop]-=dt*(w_im[i][k]*(fac2*SfHv_t[i][k]+fac3*SfHv_tm1[i][k])-
                       w_im[i][k+1]*(fac2*SfHv_t[i][k+1]+fac3*SfHv_tm1[i][k+1]))-
          dt*(bd[k]*(fac2*phys->stmp[i][k-1]+fac3*scal_old[i][k-1])-
              (bd[k]+bd[k+1])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])+
              bd[k+1]*(fac2*phys->stmp[i][k+1]+fac3*scal_old[i][k+1]));	
      }
    //top boundary conditions
    if(ktop<grid->Nk[i]-1) {
      k=ktop;
      r[0]-=dt*(w_im[i][k]*(fac2*SfHv_t[i][k]+fac3*SfHv_tm1[i][k])-w_im[i][k+1]*(fac2*SfHv_t[i][k+1]+fac3*SfHv_tm1[i][k+1]))-
        dt*(-(2*alpha_top*bd[k+1]+bd[k+1])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])+
            bd[k+1]*(fac2*phys->stmp[i][k+1]+fac3*scal_old[i][k+1]));
    }
    //bottom boundary condition
    k=grid->Nk[i]-1; 	 
    r[k-ktop]-=dt*(w_im[i][k]*(fac2*SfHv_t[i][k]+fac3*SfHv_tm1[i][k])-
		   w_im[i][k+1]*(fac2*SfHv_t[i][k+1]+fac3*SfHv_tm1[i][k+1]))+
      dt*(bd[k]*(fac2*phys->stmp[i][k-1]+fac3*scal_old[i][k-1])-
	  (bd[k]+2*alpha_bot*bd[k])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k]));
    
    }

    //printf("through vert source \n");
    
    // If this code is uncommented need to be careful with kappa_tv which can be NULL
    /*if(!prop->subgrid){
      for(k=ktop+1;k<grid->Nk[i]-1;k++) 
        d[k-ktop]-=dt*(am[k]*(fac2*phys->stmp[i][k-1]+fac3*scal_old[i][k-1])+
            (ap[k]-am[k+1])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])-
            ap[k+1]*(fac2*phys->stmp[i][k+1]+fac3*scal_old[i][k+1]))-
            dt*(bd[k]*(fac2*phys->stmp[i][k-1]+fac3*scal_old[i][k-1])-
              (bd[k]+bd[k+1])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])+
              bd[k+1]*(fac2*phys->stmp[i][k+1]+fac3*scal_old[i][k+1]));
    } else {
      // recalculate bd by acveffold
      for(k=ktop+1;k<grid->Nk[i];k++)
      {
        alpha=subgrid->Acveffold[i][k]/Ac;
        bd[k]=alpha*(2.0*kappa+kappa_tv[i][k-1]+kappa_tv[i][k])/
          (grid->dzz[i][k-1]+grid->dzz[i][k]);
      }
      for(k=ktop+1;k<grid->Nk[i]-1;k++) 
        d[k-ktop]-=(1-theta)*dt*(am[k]*phys->stmp[i][k-1]+
            (ap[k]-am[k+1])*phys->stmp[i][k]-
            ap[k+1]*phys->stmp[i][k+1])-
          (1-theta)*dt*(bd[k]*phys->stmp[i][k-1]
              -(bd[k]+bd[k+1])*phys->stmp[i][k]
              +bd[k+1]*phys->stmp[i][k+1]);
    }*/


    if(ktop<grid->Nk[i]-1) {
      //Flux through bottom of top cell
      k=ktop;
      if(!(prop->TVD>=5) || prop->vertTVD==0){
        r[0]=r[0]-dt*((ap[k]+am[k]-am[k+1])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])-
		  ap[k+1]*(fac2*phys->stmp[i][k+1]+fac3*scal_old[i][k+1]))+
      dt*(-(2*alpha_top*bd[k+1]+bd[k+1])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])+
	  bd[k+1]*(fac2*phys->stmp[i][k+1]+fac3*scal_old[i][k+1]));
      } 
       
      // subgrid ??? 
      if(Ftop) r[0]+=dt*(1-alpha_top+2*alpha_top*bd[k+1])*Ftop[i];
    
      // Through top of bottom cell
      k=grid->Nk[i]-1;
      if(!(prop->TVD>=5) || prop->vertTVD==0){
	      r[k-ktop]-=dt*(am[k]*(fac2*phys->stmp[i][k-1]+fac3*scal_old[i][k-1])+
		       (ap[k]-ap[k+1])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k]))-
          dt*(bd[k]*(fac2*phys->stmp[i][k-1]+fac3*scal_old[i][k-1])-
	      (bd[k]+2*alpha_bot*bd[k])*(fac2*phys->stmp[i][k]+fac3*scal_old[i][k])); 
      }
      // subgrid ???
      if(Fbot) r[k-ktop]+=dt*(-1+alpha_bot+2*alpha_bot*bd[k])*Fbot[i];
    }
    // First add on the source term from the previous time step (not currently in use)
    if(grid->ctop[i]<=grid->ctopold[i]) {
      for(k=grid->ctop[i];k<=grid->ctopold[i];k++){
        r[0]+=(1-fab)*Cn[i][grid->ctopold[i]]/(1+abs(grid->ctop[i]-grid->ctopold[i]));
      }
      for(k=grid->ctopold[i]+1;k<grid->Nk[i];k++) {
        r[k-grid->ctopold[i]]+=(1-fab)*Cn[i][k];
      }
    } else {
      for(k=grid->ctopold[i];k<=grid->ctop[i];k++) {
        r[0]+=(1-fab)*Cn[i][k];
      }
	for(k=grid->ctop[i]+1;k<grid->Nk[i];k++){ 
	  r[k-grid->ctop[i]]+=(1-fab)*Cn[i][k];
	}
    }

    for(k=0;k<grid->ctop[i];k++)
      Cn[i][k]=0;

    if(src2)
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
      { 
        if(prop->subgrid)
          alpha=subgrid->Acceff[i][k]/Ac;
        else 
          alpha=1.0; 

        Cn[i][k-ktop]=alpha*dt*src2[i][k]*grid->dzz[i][k];
      }
    else
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	      Cn[i][k]=0;

    // Now create the source term for the current time step
    for(k=0;k<grid->Nk[i];k++)
      ap[k]=0;

    if(prop->vertcoord==4){
      z_cutoff = 0;
      //should not have these hardcoded forever
      num_bot_layers = 4;
      num_mid_layers = 10 + 1;
      bot_cutoff = -280;

      //should do in outside loop and save for each i? 
      for(k=0;k<grid->Nk[i]-(num_bot_layers+num_mid_layers);k++){
        z_cutoff-=grid->dzz[i][k];
      }
    }


    //loop over the faces of each cell
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      df = grid->df[ne];
      dg = grid->dg[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) 
       {
        nc2=nc1;
        if(boundary_scal && (grid->mark[ne]==2 || grid->mark[ne]==3))
          sp=phys->stmp2[nc1];
        else
          sp=phys->stmp[nc1];
      } else 
        sp=phys->stmp[nc2];

      if(!(prop->TVD && prop->horiTVD)) {
        for(k=0;k<grid->Nke[ne];k++) {
          
          temp[k]=UpWind((prop->imfac2*phys->u_old[ne][k]+prop->imfac1*phys->u[ne][k]+prop->imfac3*phys->u_old2[ne][k]),phys->stmp[nc1][k],sp[k]);
          // new edit
          if(k<grid->ctopold[nc1])
            temp[k]=sp[k];
          if(k<grid->ctopold[nc2])
            temp[k]=phys->stmp[nc1][k];
        }

      } else {
        for(k=0;k<grid->Nke[ne];k++) { //set temp=phi{face}, choosing + or - based on sign of u(n+theta)
          if((prop->imfac2*phys->u_old[ne][k]+prop->imfac1*phys->u[ne][k]+prop->imfac3*phys->u_old2[ne][k])>0)
            temp[k]=phys->SfHp[ne][k];
          else
            temp[k]=phys->SfHm[ne][k];

          //here, calculate phi_ex
        if(timestepping){ //only keep track/update face values if timestepping==1
	        if(prop->n==prop->nstart+1){
            if(prop->readOldVelocity){
              HorizontalFaceScalars_oldstep(grid, phys, prop, phys->stmp, SfH_tm1, boundary_scal, prop->TVD, 
                1, 1,comm, myproc);
              HorizontalFaceScalars_oldstep(grid, phys, prop, scal_old, SfH_tm2, boundary_scal, prop->TVD, 
                2, 1,comm, myproc);
            } else {
              SfH_tm1[ne][k] = 0;
              SfH_tm2[ne][k] = 0;
            }
	        }
          temp_k = fab1*temp[k] + fab2*SfH_tm1[ne][k] + fab3*SfH_tm2[ne][k];

          //save old face values for next time step
          SfH_tm2[ne][k] = SfH_tm1[ne][k];
          SfH_tm1[ne][k] = temp[k];
	        temp[k] = temp_k;
        }

	  // OBF 042623 Central differencing
	  //	  temp[k]=0.5*(phys->stmp[nc1][k]+sp[k]);	  
	  }
      }
      

      // this part don't need to be changed
      // this is first term on RHS in new suntans paper eq 61. Horizontal advection summation
      sumz=0;
      for(k=0;k<grid->Nke[ne];k++){
        ap[k]+= dt*df*normal/Ac*(prop->imfac1*phys->u[ne][k]+prop->imfac2*phys->u_old[ne][k]+prop->imfac3*phys->u_old2[ne][k])
          *temp[k]*grid->dzf[ne][k];
        //try simplest version of horizontal diffusion 
        if(prop->vertcoord==4){
          sumz-=.5*grid->dzf[ne][k];

          // cent = .5*(.5*(z_cutoff(nc1)+z_cutoff(nc2)) + bot_cutoff); //put middle of tanh in middle of transition layers
          // span = (1/8)*(.5*(z_cutoff(nc1)+z_cutoff(nc2)) - bot_cutoff); //tanh reaches maxes in about +- 4x(span)

          //should average z cutoff over the edge, but for now we'll test just using (i)
          cent = .5*((z_cutoff) + bot_cutoff); //put middle of tanh in middle of transition layers
          span = (1/8)*((z_cutoff) - bot_cutoff); //tanh reaches maxes in about +- 4x(span)

          kappaH_k = .5*(1+tanh((sumz-cent)/span));
          sumz-=.5*grid->dzf[ne][k];
        } else {
          kappaH_k = kappaH;
        }
        kappaH_k = kappaH;
        Cn[i][k]+= dt*df*normal/Ac*kappaH_k/dg*(phys->stmp[nc1][k]-phys->stmp[nc2][k])*grid->dzf[ne][k];
      }

      //here we want to add in horizontal diffusion, which is a similar expression to above. 

    }

    // for(k=ktop+1;k<grid->Nk[i];k++) //subtract horizontal advection from Cn (source term, set to 0 at start of func) 
    //   Cn[i][k-ktop]-=ap[k];

    // for(k=0;k<=ktop;k++) 
    //   Cn[i][0]-=ap[k];

    // Add on the source from the current time step to the rhs.
    for(k=0;k<grid->Nk[i]-ktop;k++){ //fab currently=1, so d+=Cn. Recall d already includes explicit vertical advection and h(n)phi(n)
      r[k]+=fab*Cn[i][k] - ap[k];
      ap[k]=Cn[i][k]; //changed from a negative
    }
    //printf("through this loop \n");
    
    // Add on the volume correction if h was < -d
    /*
       if(grid->ctop[i]==grid->Nk[i]-1)
       d[grid->Nk[i]-ktop-1]+=phys->hcorr[i]*phys->stmp[i][grid->ctop[i]];
       */

    for(k=ktop;k<grid->Nk[i];k++)
      ap[k]=Cn[i][k-ktop];
    for(k=0;k<=ktop;k++)
      Cn[i][k]=0;
    for(k=ktop+1;k<grid->Nk[i];k++)
      Cn[i][k]=ap[k];
    for(k=grid->ctop[i];k<=ktop;k++)
      Cn[i][k]=ap[ktop]/(1+abs(grid->ctop[i]-ktop));

    if(grid->Nk[i]-ktop>1){ //solve the implicit problem, set scal
      if(prop->TVD>=5 && prop->vertTVD==1){
	    ReducePentadiag(a, b, c, d, e, r, &(scal[i][ktop]),grid->Nk[i]-ktop);

	    //options for Gauss Seidel solver, could put into suntans.dat or somwhere more tunable.
	      //REAL tol = 1E-9;
	      //int maxiters = 1000;
	      //GaussSeidel(a, b, c, d, e, r, &(scal[i][ktop]), phys->stmp, i, tol, maxiters, grid->Nk[i]-ktop);
      }else{
	      TriSolve(b,c,d,r,&(scal[i][ktop]),grid->Nk[i]-ktop);
      }
    }else if(prop->n>1) {
      if(c[0]>0 && phys->active[i])
        scal[i][ktop]=r[0]/c[0];
      else 
        scal[i][ktop]=phys->stmp[i][ktop];
    }

    // Don't update advection (of momentum!) in small cells
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      if(grid->dzz[i][k]<prop->BUFFERHEIGHT) {
	      scal[i][k]=phys->stmp[i][k];
      }
    }

    if(grid->Nk[i]-ktop>1 && !phys->active[i]){
      for(k=ktop;k<grid->Nk[i];k++)
        scal[i][k]=phys->stmp[i][k];
    }
    
    // subgrid flux check for each layer of cell
    // this may induce mass loss!!!!
    // the reason that this is necessary is the wet-dry condition for scalar transport
    // is never Courant number=1 due to the varying Volume/flux height ratio
    if(prop->subgrid)
      if(ktop!=grid->Nk[i]-1)
        for(k=ktop;k<grid->Nk[i];k++){
          if(subgrid->fluxn[i][k]>grid->dzz[i][k]*subgrid->Acceff[i][k])
            scal[i][k]=phys->stmp[i][k];
        }

    for(k=0;k<grid->ctop[i];k++)
      scal[i][k]=0;

    for(k=grid->ctop[i];k<grid->ctopold[i];k++) 
      scal[i][k]=scal[i][ktop];
    
    if(prop->TVD>=5){
      for(k=0;k<grid->Nk[i]+1;k++){
	      SfHv_tm1[i][k] = SfHv_t[i][k];
      }
      //we can calculate new scalar on face here if we don't include the transverse term (calling SumNeighborScalars), but otherwise need to do a 
      //GetScalarOnFace(grid, phys, SfHv_t, phys->wp, phys->wm, w_im, grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD, BCtop, BCbot,, comm, myproc);
    }
    // update scal^old
    for(k=0;k<grid->Nk[i];k++) {
      scal_old[i][k]=phys->stmp[i][k];
    }
  }

  CheckDivergence(grid->dzz,grid->dzzold,phys->u,phys->u_old,phys->u_old2,
    w_im,grid,prop,1e-5,myproc);
    
  if(prop->TVD>=5){
    //calculate transverse curvature correction for vertical face interpolation of newly calculate scalar
    SumNeighborScalars(grid, phys, scal, sum_neighs,
		       comm, myproc);
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr]; //cell index 
      //interpolate scalar onto vertical faces for use in next time loop
      GetScalarOnFace(SfHv_t, phys->wp, phys->wm, w_im, grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD, BCtop, BCbot, sum_neighs);
    }  
  } 

  // Code to check divergence change CHECKCONSISTENCY to 1 in suntans.h
  if(0){
    //  if(CHECKCONSISTENCY && checkflag) {

    if(prop->n==1+prop->nstart) {
      smin=INFTY;
      smax=-INFTY;
      for(i=0;i<grid->Nc;i++) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          if(phys->stmp[i][k]>smax) { 
            smax=phys->stmp[i][k]; 
            imax=i; 
            kmax=k; 
          }
          if(phys->stmp[i][k]<smin) { 
            smin=phys->stmp[i][k]; 
            imin=i; 
            kmin=k; 
          }
        }
      }
      MPI_Reduce(&smin,&smin_value,1,MPI_DOUBLE,MPI_MIN,0,comm);
      MPI_Reduce(&smax,&smax_value,1,MPI_DOUBLE,MPI_MAX,0,comm);
      MPI_Bcast(&smin_value,1,MPI_DOUBLE,0,comm);
      MPI_Bcast(&smax_value,1,MPI_DOUBLE,0,comm);

      if(myproc==0)
        printf("Minimum scalar: %.2f, maximum: %.2f\n",smin_value,smax_value);
    }      

    div_local_count=0;
    div_da_count=0;
    //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];

      flag=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        if(grid->mark[grid->face[i*grid->maxfaces+nf]]==2 || 
            grid->mark[grid->face[i*grid->maxfaces+nf]]==3) {
          flag=1;
          break;
        }
      }

      if(!flag) {
        div_da=0;

        for(k=0;k<grid->Nk[i];k++) {

          div_local=0;
          for(nf=0;nf<grid->nfaces[i];nf++) {
            ne=grid->face[i*grid->maxfaces+nf];
            div_local+=(prop->imfac1*phys->u[ne][k]+prop->imfac2*phys->u_old[ne][k]+prop->imfac3*phys->u_old2[ne][k])
              *grid->dzf[ne][k]*grid->normal[i*grid->maxfaces+nf]*grid->df[ne];
          }
	  div_local+=grid->Ac[i]*(grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt;	  	  
	  div_local+=grid->Ac[i]*(w_im[i][k]-w_im[i][k+1]);
          div_da+=div_local;
	  
          if(k>=grid->ctop[i]) {
            if(fabs(div_local)>SMALL_CONSISTENCY) {// && grid->dzz[imin][0]>DRYCELLHEIGHT) {
	      div_local_max = Max(div_local_max,fabs(div_local));
	      div_local_count++;
	      //	      exit(1);
	      //	      printf("Step: %d, proc: %d, locally-divergent at %d, %d, div=%e\n",
	      //		     prop->n,myproc,i,k,div_local);
	    }
          }
        }
        if(fabs(div_da)>SMALL_CONSISTENCY) {// && phys->h[i]+grid->dv[i]>DRYCELLHEIGHT) {
	  div_da_max = Max(div_da_max,fabs(div_da));
	  div_da_count++;
	  //          printf("i %d  Step: %d, proc: %d, Depth-Ave divergent at i=%d, div=%e\n",
	  //              i,prop->n,myproc,i,div_da);
      }
    }
    }

    if(div_local_count>0)
      printf("Proc %d, Max local divergence = %.3e, count exceeding %.3e = %d\n",
	     myproc,div_local_max,SMALL_CONSISTENCY,div_local_count);
    if(div_da_count>0)
      printf("Proc: %d, Max depth-averaged divergence = %.3e, count exceeding %.3e = %d\n",
	     myproc,div_da_max,SMALL_CONSISTENCY,div_da_count);

    mincount=0;
    maxcount=0;
    smin=INFTY;
    smax=-INFTY;
    //for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];

      flag=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        if(grid->mark[grid->face[i*grid->maxfaces+nf]]==2 || grid->mark[grid->face[i*grid->maxfaces+nf]]==3) 
        {
          flag=1;
          break;
        }
      }

      if(!flag) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          if(scal[i][k]>smax) { 
            smax=scal[i][k]; 
            imax=i; 
            kmax=k; 
          }
          if(scal[i][k]<smin) { 
            smin=scal[i][k]; 
            imin=i; 
            kmin=k; 
          }

          if(scal[i][k]>smax_value+SMALL_CONSISTENCY && grid->dzz[i][k]>DRYCELLHEIGHT)
            maxcount++;
          if(scal[i][k]<smin_value-SMALL_CONSISTENCY && grid->dzz[i][k]>DRYCELLHEIGHT)
            mincount++;
        }
      }
    }
    MPI_Reduce(&mincount,&allmincount,1,MPI_INT,MPI_SUM,0,comm);
    MPI_Reduce(&maxcount,&allmaxcount,1,MPI_INT,MPI_SUM,0,comm);

    if(mincount!=0 || maxcount!=0) 
      printf("Not CWC, step: %d, proc: %d, smin = %e at i=%d,H=%e, smax = %e at i=%d,H=%e\n",
          prop->n,myproc,
          smin,imin,phys->h[imin]+grid->dv[imin],
          smax,imax,phys->h[imax]+grid->dv[imax]);

    if(myproc==0 && (allmincount !=0 || allmaxcount !=0))
      printf("Total number of CWC violations (all procs): s<s_min: %d, s>s_max: %d\n",
          allmincount,allmaxcount);
  }

}

/* Moved to util.c
void ReducePentadiag(REAL *a, REAL *b, REAL *c, REAL *d, REAL *e, REAL *r, REAL *u, int N)
{
  // basically standard tridiagonal solver as at 
  // http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  int k;
  REAL pivot;


  //Remove entries below diagonal one column at a time using diagonal as pivot (diagonal is never 0, but others could be)
  for(k=0;k<N-2;k++){
    pivot=b[k+1]/c[k];
    //b[k+1]-=pivot*c[k];
    c[k+1]-=pivot*d[k];
    d[k+1]-=pivot*e[k];
    r[k+1]-=pivot*r[k];


    pivot=a[k+2]/c[k];
    //a[k+2]-=pivot*c[k];
    b[k+2]-=pivot*d[k];
    c[k+2]-=pivot*e[k];
    r[k+2]-=pivot*r[k];
  }

  k=N-2;
  pivot=b[k+1]/c[k];
  //b[k+1]-=pivot*c[k];
  c[k+1]-=pivot*d[k];
  d[k+1]-=pivot*e[k];
  r[k+1]-=pivot*r[k];


  d[N-1] = 0; //there is no N-1 entry so reset that to 0 again 

  //now solve
  for(k=0;k<N;k++){
    d[k]=d[k]/c[k];
    e[k] = e[k]/c[k];
    r[k]=r[k]/c[k];
    c[k]=1;
  }
  d[N-1] = 0; //there is no N-1 entry so reset that to 0 again 
  
  u[N-1]=r[N-1];
  u[N-2]=r[N-2]-d[N-2]*u[N-1];

  for(k=N-2;k>=0;k--){
    u[k] = r[k] - d[k]*u[k+1] - e[k]*u[k+2]; 
  }


}
*/
void GaussSeidel(double *a, double *b, double *c, double *d, double *e, double *r, double *u, double **u_init, int i, double tol, int maxiters, int N)
{
  int k, iters=0;
  REAL errork = 0;
  REAL checktol = 100;

  for(k=0;k<N;k++){
    u[k] = u_init[i][k];
  }

  while(checktol>tol && iters<maxiters){
    iters = iters+1;
        
    for(k=0;k<N;k++){
      u[k] = 1/c[k]*(r[k] -(a[k]*u[k-2]+b[k]*u[k-1]+d[k]*u[k+1]+e[k]*u[k+2]));
    }

    checktol = 0;
    for(k=0;k<N;k++){
      errork = (r[k] -(a[k]*u[k-2]+b[k]*u[k-1]+c[k]*u[k]+d[k]*u[k+1]+e[k]*u[k+2]));
      checktol += errork*errork;
    }
    checktol = sqrt(checktol);
    if(iters==maxiters){
      printf("did not converge after maxiters, current error %f \n", checktol);
      exit(0);
    }  
  }

}
