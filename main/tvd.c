/*
 * File: tvd.c
 * Author: Zhonghua Zhang and Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains the TVD scalar-computing functions.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "tvd.h"
#include "util.h"
#include "sendrecv.h"
#include <math.h>

// To prevent denominators from going to zero.
#define EPS 1e-12

// Local function
static REAL Psi(REAL r, int TVD);
static REAL GetCF(REAL phi_c);
/*
 * Function: HorizontalFaceScalars
 * Usage: HorizontalFaceScalars(grid, phys, boundary_scal);
 * ---------------------------------------------------------------------------
 * Calculate the horizontal face scalars with upwind/TVD schemes.  
 * SfHp[Ne][Nk] & SfHm[Ne][Nk] are used to store the scalar facial values.
 * where S--scalar, f--face, H--horizontal, p--plus, m--minus;
 *       Ne--the number of horizontal edges, Nk--the number of vertical layers. 
 *
 */
void HorizontalFaceScalars(gridT *grid, physT *phys, propT *prop, REAL **scal, REAL **U, REAL **boundary_scal, int TVD, 
			   int dimen, MPI_Comm comm, int myproc) 
{
  int i, iptr, j, k, m, mf, jptr, ib, nc1, nc2, ne, neigh, normal;
  REAL u_nptheta, Qminus, r, si, sm,fac1,fac2,fac3, rprime, phi_c, CF;
  REAL **sumQC, **sumQ;

  // set new implicit scheme here
  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  // Pointer Variables for TVD scheme
  //Casulli method flux weighted
  sumQ = phys->gradSx;
  sumQC = phys->gradSy;
  
  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      sumQ[i][k]=0;
      sumQC[i][k]=0;

      for(mf=0;mf<grid->nfaces[i];mf++) {
	ne = grid->face[i*grid->maxfaces+mf];
	neigh = grid->neigh[i*grid->maxfaces+mf];
	normal = grid->normal[i*grid->maxfaces+mf];

	u_nptheta = normal*(fac1*U[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k]);
	if(dimen==2)
	  Qminus = 0.5*grid->df[ne]*fabs(u_nptheta-fabs(u_nptheta));
	else
	Qminus = 0.5*grid->dzf[ne][k]*grid->df[ne]*fabs(u_nptheta-fabs(u_nptheta));

	if(neigh!=-1)
	  sumQC[i][k]+=Qminus*(scal[i][k]-scal[neigh][k]);
	sumQ[i][k]+=Qminus;
      }
    }
  }
  ISendRecvCellData3D(sumQ,grid,myproc,comm);  
  ISendRecvCellData3D(sumQC,grid,myproc,comm);  
 
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    for(k=0;k<grid->etop[j];k++) 
      phys->SfHp[j][k] = phys->SfHm[j][k] = 0;
      
    for(k=grid->etop[j];k<grid->Nke[j];k++) {

      u_nptheta = fac1*U[j][k]+fac2*phys->u_old[j][k]+fac3*phys->u_old2[j][k];

      if(u_nptheta>0) {
	i=nc2;
	m=nc1;

	si=scal[nc2][k];
	sm=scal[nc1][k];
      } else {
	i=nc1;
	m=nc2;
	
	si=scal[nc1][k];
	sm=scal[nc2][k];
      }

      if(sumQ[i][k]!=0 && sm!=si)
	r=sumQC[i][k]/(sumQ[i][k]*(sm-si));
      else 
	r=0;

  if(TVD>=5) { //phi_c of 1.6 will give CF=1/8
    if(nc1== -1 || nc2== -1){
      phi_c = 1.6;
    } else if(sumQ[i][k]==0 || sumQ[m][k]==0){
      phi_c = 1.6;
    }else if((sumQC[i][k]/sumQ[i][k])==-(sumQC[m][k]/sumQ[m][k])){
      phi_c = 1.6;
    }else if(fabs((sumQC[i][k]/sumQ[i][k])+(sumQC[m][k]/sumQ[m][k]))<=EPS){
      phi_c = 1.6;
    }else if(sumQC[i][k]==0){
      phi_c = 1.6;
          } else{
      phi_c = (sumQC[i][k]/sumQ[i][k])/((sumQC[i][k]/sumQ[i][k])+(sumQC[m][k]/sumQ[m][k]));
    }

    if(isnan(phi_c)) {
      phi_c = 1.6; 
    }
    
    if(sumQ[i][k]==0){
      rprime=0;
    } else{
      rprime=sumQC[i][k]/(sumQ[i][k]);	
    }

    if(TVD==6){
      CF = GetCF(phi_c);
    } else{
      CF=1/8;
    }

    //only works for BCS where everything is no flux at top and bottom, need to update to include other BCS
    if(k==grid->etop[j]){
      //phys->SfHp[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc2][k]-2*scal[nc2][k]+scal[nc2][k+1]);
      //phys->SfHm[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc1][k]-2*scal[nc1][k]+scal[nc1][k+1]);

      //without transverse
      phys->SfHp[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
      phys->SfHm[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
    }else if(k==(grid->Nke[j]-1)){
      //phys->SfHp[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc2][k-1]-2*scal[nc2][k]+scal[nc2][k]);
      //phys->SfHm[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc1][k-1]-2*scal[nc1][k]+scal[nc1][k]);

      //without transverse
      phys->SfHp[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
      phys->SfHm[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
    }else{
      //Update step                                                                                                                                             
      //phys->SfHp[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc2][k-1]-2*scal[nc2][k]+scal[nc2][k+1]);
      //phys->SfHm[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc1][k-1]-2*scal[nc1][k]+scal[nc1][k+1]);

      //without transverse
      phys->SfHp[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
      phys->SfHm[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
    }
	
	//Update step (no transverse term)
	// phys->SfHp[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
  // phys->SfHm[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;

    } else{ //normal TVD 
      phys->SfHp[j][k] = scal[nc2][k]+0.5*Psi(r,TVD)*(scal[nc1][k]-scal[nc2][k]);
      phys->SfHm[j][k] = scal[nc1][k]-0.5*Psi(r,TVD)*(scal[nc1][k]-scal[nc2][k]);
    }
    }
  }

  // Type 2/3 boundary specifies flux at faces
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[4];jptr++) {
    j = grid->edgep[jptr];
    printf("any edges are this \n");
    
    ib = grid->grad[2*j];
    
    for(k=0;k<grid->etop[j];k++)
      phys->SfHp[j][k] = phys->SfHm[j][k] = 0;
    
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->SfHp[j][k] = boundary_scal[jptr-grid->edgedist[2]][k];  // Coming in if u>0
      phys->SfHm[j][k] = scal[ib][k];                               // Going out if u<0
    }
  }

}

void HorizontalFaceScalars_oldstep(gridT *grid, physT *phys, propT *prop, REAL **scal, REAL **faceval, REAL **boundary_scal, int TVD, 
			   int timestep, int dimen, MPI_Comm comm, int myproc) 
{
  int i, iptr, j, k, m, mf, jptr, ib, nc1, nc2, ne, neigh, normal;
  REAL u_nptheta, Qminus, r, si, sm,fac1,fac2,fac3, rprime, phi_c, CF, dzf;
  REAL **sumQC, **sumQ, **dzz_old2, **dzf_old, **dzf_old2;

  // Pointer Variables for TVD scheme
  //Casulli method 
  sumQ = phys->gradSx;
  sumQC = phys->gradSy;

  fac1 = prop->imfac1;
  fac2 = prop->imfac2;
  fac3 = prop->imfac3;

  if(dimen!=2){ //need old dzfs 

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
    if(prop->vertcoord==4){
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

    //as test
    //wait is this the oscillation problem?
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      if(nc1==-1)
	nc1=nc2;
      if(nc2==-1)
	nc2=nc1;

      for(k=0;k<grid->Nkmax;k++) {
        dzf_old[j][k]=0.5*(grid->dzzold[nc1][k]+grid->dzzold[nc2][k]);
        dzf_old2[j][k]=0.5*(dzz_old2[nc1][k]+dzz_old2[nc2][k]);
      }
    }
  }
  
  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      sumQ[i][k]=0;
      sumQC[i][k]=0;

      for(mf=0;mf<grid->nfaces[i];mf++) {
	ne = grid->face[i*grid->maxfaces+mf];
	neigh = grid->neigh[i*grid->maxfaces+mf];
	normal = grid->normal[i*grid->maxfaces+mf];

	u_nptheta = normal*(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k]);
  if(timestep==1){
    u_nptheta=phys->u_old[ne][k];
  }else if(timestep==2){
    u_nptheta=phys->u_old2[ne][k];
  }

  if(dimen==2)
	  Qminus = 0.5*grid->df[ne]*fabs(u_nptheta-fabs(u_nptheta));
	else{
    if(timestep==1)
      dzf = dzf_old[ne][k];
    else if(timestep==2)
      dzf = dzf_old2[ne][k];
    
	  //Qminus = 0.5*grid->dzf[ne][k]*grid->df[ne]*fabs(u_nptheta-fabs(u_nptheta));
    Qminus = 0.5*dzf*grid->df[ne]*fabs(u_nptheta-fabs(u_nptheta));
    //Qminus = 0.5*grid->df[ne]*fabs(u_nptheta-fabs(u_nptheta));

  }
	//Qminus = 0.5*grid->dzf[ne][k]*grid->df[ne]*fabs(u_nptheta-fabs(u_nptheta));

	if(neigh!=-1)
	  sumQC[i][k]+=Qminus*(scal[i][k]-scal[neigh][k]);
	sumQ[i][k]+=Qminus;
      }
    }
  }
  ISendRecvCellData3D(sumQ,grid,myproc,comm);  
  ISendRecvCellData3D(sumQC,grid,myproc,comm);  
 
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];

    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    for(k=0;k<grid->etop[j];k++) 
      faceval[j][k] =  0;
      
    for(k=grid->etop[j];k<grid->Nke[j];k++) {

      if(timestep==1){
        u_nptheta=phys->u_old[j][k];
      }else if(timestep==2){
        u_nptheta=phys->u_old2[j][k];
      }

      if(u_nptheta>0) {
	i=nc2;
	m=nc1;

	si=scal[nc2][k];
	sm=scal[nc1][k];
      } else {
	i=nc1;
	m=nc2;
	
	si=scal[nc1][k];
	sm=scal[nc2][k];
      }

      if(sumQ[i][k]!=0 && sm!=si)
	r=sumQC[i][k]/(sumQ[i][k]*(sm-si));
      else 
	r=0;

  if(TVD>=5) { //phi_c of 1.6 will give CF=1/8
    if(nc1== -1 || nc2== -1){
      phi_c = 1.6;
    } else if(sumQ[i][k]==0 || sumQ[m][k]==0){
      phi_c = 1.6;
    }else if((sumQC[i][k]/sumQ[i][k])==-(sumQC[m][k]/sumQ[m][k])){
      phi_c = 1.6;
    }else if(fabs((sumQC[i][k]/sumQ[i][k])+(sumQC[m][k]/sumQ[m][k]))<=EPS){
      phi_c = 1.6;
    }else if(sumQC[i][k]==0){
      phi_c = 1.6;
          } else{
      phi_c = (sumQC[i][k]/sumQ[i][k])/((sumQC[i][k]/sumQ[i][k])+(sumQC[m][k]/sumQ[m][k]));
    }

    if(isnan(phi_c)) {
      phi_c = 1.6; 
    }
    
    if(sumQ[i][k]==0){
      rprime=0;
    } else{
      rprime=sumQC[i][k]/(sumQ[i][k]);	
    }

    if(TVD==6){
      CF = GetCF(phi_c);
    } else{
      CF=1/8;
    }

    //only works for BCS where everything is no flux at top and bottom, need to update to include other BCS
    //ie try for 0 condition on w?
    if(u_nptheta>0){
      if(k==grid->etop[j]){
        //faceval[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc2][k]-2*scal[nc2][k]+scal[nc2][k+1]);
        faceval[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
      }else if(k==(grid->Nke[j]-1)){
        //faceval[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc2][k-1]-2*scal[nc2][k]+scal[nc2][k]);
        faceval[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
      }else{
        //Update step                                                                                                                                             
        //faceval[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc2][k-1]-2*scal[nc2][k]+scal[nc2][k+1]);
        faceval[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;

      }
    }else{
      if(k==grid->etop[j]){
        //faceval[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc1][k]-2*scal[nc1][k]+scal[nc1][k+1]);
        faceval[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;

      }else if(k==(grid->Nke[j]-1)){
        //faceval[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc1][k-1]-2*scal[nc1][k]+scal[nc1][k]);
        faceval[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
      }else{
        //Update step                                                                                                                                             
        //faceval[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime + 1/24*(scal[nc1][k-1]-2*scal[nc1][k]+scal[nc1][k+1]);
        faceval[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;

      }
    }
	
	//Update step (no transverse term)
  // if(u_nptheta>0){
  //   faceval[j][k]=scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
  // }else{
  //   faceval[j][k]=scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
  // }
	// phys->SfHp[j][k] = scal[nc2][k]+(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;
  // phys->SfHm[j][k] = scal[nc1][k]-(0.5-CF)*(scal[nc1][k]-scal[nc2][k])+CF*rprime;

    } else{ //normal TVD 
      if(u_nptheta>0){
        faceval[j][k] = scal[nc2][k]+0.5*Psi(r,TVD)*(scal[nc1][k]-scal[nc2][k]);
      } else{
        faceval[j][k] = scal[nc1][k]-0.5*Psi(r,TVD)*(scal[nc1][k]-scal[nc2][k]);
      }
    }
    }
  }

  // Type 2/3 boundary specifies flux at faces
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[4];jptr++) {
    j = grid->edgep[jptr];
    
    ib = grid->grad[2*j];
    
    for(k=0;k<grid->etop[j];k++)
      phys->SfHp[j][k] = phys->SfHm[j][k] = 0;
    
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->SfHp[j][k] = boundary_scal[jptr-grid->edgedist[2]][k];  // Coming in if u>0
      phys->SfHm[j][k] = scal[ib][k];                               // Going out if u<0
    }
  }

  if(dimen!=2){
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
}

/*
 * Function: Psi
 * Usage: Psi(r, TVD);
 * ---------------------------------------------------------------------------
 * Calculate the filter function Psi for TVD schemes 
 * TVD=1  first-order upwind,    TVD=2  Lax-Wendroff
 * TVD=3  Superbee,              TVD=4  Van Leer  
 */
static REAL Psi(REAL r, int TVD){
  switch(TVD) {
  case 1:
    return 0;
    break;
  case 2:
    return 1;
    break;
  case 3:
    return ( Max(0, Max(Min(2*(r),1), Min(r,2))));
    break;
  case 4:
    if(r>0)
      return Min(r,1);
    else
      return 0;
    break;
  default:
    return 0;
    break;
  }
}

/*
 * Function: GetApAm
 * Usage:  GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
 *	           phys->w,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
 * ------------------------------------------------------------------------------
 * Calculate the fluxes for vertical advection using the TVD schemes.
 *
 */
void GetApAm(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
	     REAL **w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot) {
  int k;

  // Implicit vertical advection terms
  for(k=0;k<Nk+1;k++) {
    wp[k] = 0.5*(w[i][k]+fabs(w[i][k]));
    wm[k] = 0.5*(w[i][k]-fabs(w[i][k]));
  }

  //no longer needed
  // Courant number: C=w*dt/dz
  // for(k=1;k<Nk;k++) {
  //   Cp[k] = 2 * wp[k]*dt/(dzz[i][k]+ dzz[i][k-1] );
  //   Cm[k] = 2 * wm[k]*dt/(dzz[i][k]+ dzz[i][k-1] );
  // }
  // k=Nk;
  // Cp[k] = wp[k]*dt/dzz[i][k-1];
  // Cm[k] = wm[k]*dt/dzz[i][k-1];


  // Compute the upwind gradient ratios r
  for(k=ktop+2;k<Nk-1;k++) {
    rp[k-ktop] = (scal[i][k]-scal[i][k+1]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);
    rm[k-ktop] = (scal[i][k-2]-scal[i][k-1]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);
  }

  if(BCtop == 2){ //zero at boundary
    k=ktop; // (scal[-1] = -scal[0]), (scal[-2] = -3*scal[0]) for scal at top face = 0 
    rp[0] = (scal[i][ktop]-scal[i][ktop+1]+EPS) / (-2*scal[i][ktop]+EPS);
    rm[0] = 1;

    k=ktop+1; 
    rp[k-ktop] = (scal[i][k]-scal[i][k+1]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);
    rm[k-ktop] = (-4*scal[i][k-1]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);
  } else{ //no flux throuh boundary
     rp[0]= (scal[i][ktop]-scal[i][ktop+1]+EPS) / (EPS);
     rm[0]= 1; //could be 0, but for consistence with Nk argument... 
  
     //this was the BC in this function originally
     rp[1]= (scal[i][ktop+1]-scal[i][ktop+2]+EPS) / (scal[i][ktop]-scal[i][ktop+1]+EPS);
     rm[1]= EPS / (scal[i][ktop]-scal[i][ktop+1]+EPS);
  }
  
  if(BCbot == 2){ //zero at boundary
    k=Nk-1; // (scal[Nk] = -scal[Nk-1]), (scal[Nk+1] = -3*scal[Nk-1])
    rp[k-ktop] = (scal[i][k]+scal[i][k]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);
    rm[k-ktop] = (scal[i][k-2]-scal[i][k-1]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);

    k=Nk;
    rp[k-ktop] = 1;
    rm[k-ktop] = (scal[i][k-2]-scal[i][k-1]+EPS) / (scal[i][k-1]+scal[i][k+1]+EPS);
  } else { //no flux by default 
    k=Nk-1; // (scal[Nk] == scal[Nk-1] = scal[Nk+1])
    //this was the BC in this function originally
    rp[k-ktop] = (EPS) / (scal[i][k-1]-scal[i][k]+EPS);
    rm[k-ktop] = (scal[i][k-2]-scal[i][k-1]+EPS) / (scal[i][k-1]-scal[i][k]+EPS);

    k=Nk; //added in bottom term
    rp[k-ktop] = 1; //or 0
    rm[k-ktop] = (scal[i][k-2]-scal[i][k-1]+EPS) / (EPS);

  }

  //changed from starting at ktop+1 for use in the top BC. Does not make the assumption that w[]=0
  //also removed the courant number terms for consistency with Casulli and the horizontal
  for(k=ktop;k<Nk;k++) {
    am[k]= 0.5*wp[k]*Psi(rp[k-ktop], TVD)
      + wm[k]*(1-0.5*Psi(rm[k-ktop], TVD));
    ap[k]= wp[k]*(1-0.5*Psi(rp[k-ktop],TVD))
      + 0.5*wm[k]*Psi(rm[k-ktop],TVD);
  }
  //no flux through bottom
  k=Nk;
  am[k]=0;
  ap[k]=0;
}

/*
 * Function: GetApAm
 * Usage:  GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
 *	           phys->w,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
 * ------------------------------------------------------------------------------
 * Calculate the fluxes for vertical advection using the TVD schemes.
 *
 */
void GetApAmNewAdv(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
  REAL *w, REAL **dzz, REAL **scal, int j, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot) {
int k;

// Implicit vertical advection terms
for(k=0;k<Nk+1;k++) {
  wp[k] = 0.5*(w[k]+fabs(w[k]));
  wm[k] = 0.5*(w[k]-fabs(w[k]));
}

//no longer needed
// Courant number: C=w*dt/dz
// for(k=1;k<Nk;k++) {
//   Cp[k] = 2 * wp[k]*dt/(dzz[i][k]+ dzz[i][k-1] );
//   Cm[k] = 2 * wm[k]*dt/(dzz[i][k]+ dzz[i][k-1] );
// }
// k=Nk;
// Cp[k] = wp[k]*dt/dzz[i][k-1];
// Cm[k] = wm[k]*dt/dzz[i][k-1];


// Compute the upwind gradient ratios r
for(k=ktop+2;k<Nk-1;k++) {
rp[k-ktop] = (scal[j][k]-scal[j][k+1]+EPS) / (scal[j][k-1]-scal[j][k]+EPS);
rm[k-ktop] = (scal[j][k-2]-scal[j][k-1]+EPS) / (scal[j][k-1]-scal[j][k]+EPS);
}

if(BCtop == 2){ //zero at boundary
  k=ktop; // (scal[-1] = -scal[0]), (scal[-2] = -3*scal[0]) for scal at top face = 0 
  rp[0] = (scal[j][ktop]-scal[j][ktop+1]+EPS) / (-2*scal[j][ktop]+EPS);
  rm[0] = 1;

  k=ktop+1; 
  rp[k-ktop] = (scal[j][k]-scal[j][k+1]+EPS) / (scal[j][k-1]-scal[j][k]+EPS);
  rm[k-ktop] = (-4*scal[j][k-1]+EPS) / (scal[j][k-1]-scal[j][k]+EPS);
} else{ //no flux throuh boundary
  rp[0]= (scal[j][ktop]-scal[j][ktop+1]+EPS) / (EPS);
  rm[0]= 1; //could be 0, but for consistence with Nk argument... 

  //this was the BC in this function originally
  rp[1]= (scal[j][ktop+1]-scal[j][ktop+2]+EPS) / (scal[j][ktop]-scal[j][ktop+1]+EPS);
  rm[1]= EPS / (scal[j][ktop]-scal[j][ktop+1]+EPS);
}

if(BCbot == 2){ //zero at boundary
  k=Nk-1; // (scal[Nk] = -scal[Nk-1]), (scal[Nk+1] = -3*scal[Nk-1])
  rp[k-ktop] = (scal[j][k]+scal[j][k]+EPS) / (scal[j][k-1]-scal[j][k]+EPS);
  rm[k-ktop] = (scal[j][k-2]-scal[j][k-1]+EPS) / (scal[j][k-1]-scal[j][k]+EPS);

  k=Nk;
  rp[k-ktop] = 1;
  rm[k-ktop] = (scal[j][k-2]-scal[j][k-1]+EPS) / (scal[j][k-1]+scal[j][k+1]+EPS);
} else { //no flux by default 
  k=Nk-1; // (scal[Nk] == scal[Nk-1] = scal[Nk+1])
  //this was the BC in this function originally
  rp[k-ktop] = (EPS) / (scal[j][k-1]-scal[j][k]+EPS);
  rm[k-ktop] = (scal[j][k-2]-scal[j][k-1]+EPS) / (scal[j][k-1]-scal[j][k]+EPS);

  k=Nk; //added in bottom term
  rp[k-ktop] = 1; //or 0
  rm[k-ktop] = (scal[j][k-2]-scal[j][k-1]+EPS) / (EPS);

}

//changed from starting at ktop+1 for use in the top BC. Does not make the assumption that w[]=0
//also removed the courant number terms for consistency with Casulli and the horizontal
  for(k=ktop;k<Nk+1;k++) {
    am[k]= 0.5*wp[k]*Psi(rp[k-ktop], TVD)
    + wm[k]*(1-0.5*Psi(rm[k-ktop], TVD));
    ap[k]= wp[k]*(1-0.5*Psi(rp[k-ktop],TVD))
    + 0.5*wm[k]*Psi(rm[k-ktop],TVD);
  }
}


/*
 * Function: GetApAmVert
 * Usage:  GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
 *	           phys->w,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
 * ------------------------------------------------------------------------------
 * Calculate the fluxes for vertical advection of vertical velocity using the TVD schemes in new advection code.
 *
 */
void GetApAmVert(REAL *ap, REAL *am, REAL *wp, REAL *wm, REAL *Cp, REAL *Cm, REAL *rp, REAL *rm,
  REAL *w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot) {
int k;

// Implicit vertical advection terms
for(k=0;k<Nk;k++) { //cell cented for Weff
wp[k] = 0.5*(w[k]+fabs(w[k]));
wm[k] = 0.5*(w[k]-fabs(w[k]));
}

//no longer needed
// Courant number: C=w*dt/dz
// for(k=1;k<Nk;k++) {
//   Cp[k] = 2 * wp[k]*dt/(dzz[i][k]+ dzz[i][k-1] );
//   Cm[k] = 2 * wm[k]*dt/(dzz[i][k]+ dzz[i][k-1] );
// }
// k=Nk;
// Cp[k] = wp[k]*dt/dzz[i][k-1];
// Cm[k] = wm[k]*dt/dzz[i][k-1];

//not super sure about BCS

// Compute the upwind gradient ratios r
//r's are also cell centered now, but scal is on edges. There is a value for w[0]
for(k=ktop+1;k<Nk-2;k++) {
  rp[k-ktop] = (scal[i][k+1]-scal[i][k+2]+EPS) / (scal[i][k]-scal[i][k+1]+EPS);
  rm[k-ktop] = (scal[i][k-1]-scal[i][k]+EPS) / (scal[i][k]-scal[i][k+1]+EPS);
}

if(BCtop == 2){ //zero at boundary
  k=ktop; // (scal[-1] = -scal[0]), (scal[-2] = -3*scal[0]) for scal at top face = 0 
  rp[0] = (scal[i][k+1]-scal[i][k+2]+EPS) / (scal[i][k]-scal[i][k+1]+EPS);
  rm[0] = (-scal[i][k]-scal[i][k]+EPS) / (scal[i][k]-scal[i][k+1]+EPS); //I think this is right, scal[-1] became -scal[0]
} else{ //no flux throuh boundary
  k=ktop;
  rp[0] = (scal[i][k+1]-scal[i][k+2]+EPS) / (scal[i][k]-scal[i][k+1]+EPS);
  rm[0] = EPS; //(scal[-1] = scal[0])
}

 //no flux by default. 
k=Nk-2; //(scal[Nk] = -scal[Nk-1]), (scal[Nk+1] = -3*scal[Nk-1]) for scal at top face = 0 
rp[k-ktop] = (scal[i][k+1]+EPS) / (scal[i][k]-scal[i][k+1]+EPS);
rm[k-ktop] = (scal[i][k-1]-scal[i][k]+EPS) / (scal[i][k]-scal[i][k+1]+EPS);


//changed from starting at ktop+1 for use in the top BC. Does not make the assumption that w[]=0
//also removed the courant number terms for consistency with Casulli and the horizontal
for(k=ktop;k<Nk-1;k++) {
  am[k]= 0.5*wp[k]*Psi(rp[k-ktop], TVD)
  + wm[k]*(1-0.5*Psi(rm[k-ktop], TVD));
  ap[k]= wp[k]*(1-0.5*Psi(rp[k-ktop],TVD))
  + 0.5*wm[k]*Psi(rm[k-ktop],TVD);
}
//note r[Nk-1] never gets used 
k=Nk-1;
am[k]=0.5*(wp[k]+wm[k]);

}

/*
 * Function: GetPentDiag
 * Usage:  GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
 *	           phys->w,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
 * ------------------------------------------------------------------------------
 * Calculate the coefficients of the pentadiagonal for implicit part of vertical advection for QUICK or SHARP
 *
 */
void GetPentDiagNewAdv(REAL *pent_a, REAL *pent_b, REAL *pent_c, REAL *pent_d, REAL *pent_e, REAL *wp, REAL *wm,
  REAL *w, REAL **dzz, REAL **U, int j, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot) {
int k;

REAL CF, phi_c, phi_c1, CF1, Bu, Bd, B, hkm2, hkm1, hk, hkp1, hkp2;

CF = 1/8;
CF1 = 1/8;


for(k=0;k<Nk+1;k++) {
   wp[k] = 0.5*(w[k]+fabs(w[k]));
   wm[k] = 0.5*(w[k]-fabs(w[k]));
}
//w[k] = w at k-1/2, w[k+1] = w at k-1/2

for(k=0;k<Nk-1;k++) { //loop over k=2 to k=Nk-2, which should be the "normal" rows of the pentdiagonal
 
 if(TVD==6){
   if(wp[k]>0){
     if(U[j][k-1]!=U[j][k+1]){
       phi_c = (U[j][k] - U[j][k+1])/(U[j][k-1]-U[j][k+1]); //hopefully this doesn't change in nonuniform....
     } else{
       phi_c = 1.6; //this will make CF=1/8
     }
   } else if(wm[k] != 0) {
     if(U[j][k]!=U[j][k-2]) {
       phi_c = (U[j][k-1] - U[j][k-2])/(U[j][k]-U[j][k-2]);
     } else{
       phi_c = 1.6;
     }
   } else {
     phi_c = 1.6;
   }

   if(wp[k+1]>0){
     if((U[j][k]!=U[j][k+2])){
       phi_c1 = (U[j][k+1] - U[j][k+2])/(U[j][k]-U[j][k+2]);
     } else{
       phi_c1 = 1.6;
     }
   } else if(wm[k] != 0) {
     if(U[j][k+1]!=U[j][k-1]){
       phi_c1 = (U[j][k] - U[j][k-1])/(U[j][k+1]-U[j][k-1]);
     } else{
       phi_c1 = 1.6;
     }
   } else {
     phi_c1 = 1.6;
   }

   CF = GetCF(phi_c); //k
   CF1 = GetCF(phi_c1); //k+1
 } else {
   CF = 1/8;
   CF1 = 1/8;
 }
 //uniform 
 pent_a[k] = -(wm[k]*(-CF));
 pent_b[k] = -(wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF1));
 pent_c[k] = -(wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF1) - wm[k+1]*(.5+2*CF1));
 pent_d[k] = -(wp[k]*(-CF) - wp[k+1]*(.5+2*CF1) - wm[k+1]*(.5-CF1));
 pent_e[k] = -(-wp[k+1]*(-CF1));
}

//BC=1 (no flux), BC=2 (zero at boundary), BC=3 (no slip at boundary), 

//use QUICK at BCS
CF=1/8;


  if(BCtop == 1 || BCtop == 3) {
  //uniform 
    k = 0;
    pent_b[k]+=pent_a[k];
    pent_c[k]+=pent_b[k];
    pent_a[k]=0;
    pent_b[k]=0;
  // pent_a[k] = 0;
  // pent_b[k] = 0;
  // pent_c[k] = -(wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF)) + (wm[k]*(-CF)) + (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
  // pent_d[k] = -(wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF));
  // pent_e[k] = -(-wp[k+1]*(-CF));

  //printf("a[0]=%f, b[0]=%f, c[0]=%f, d[0]=%f, e[0]=%f \n", pent_a[k], pent_b[k], pent_b[k], pent_d[k], pent_e[k]);

  //if(Nk>1){
    // k=1; //one below top
    // pent_a[k] = 0;
    // pent_b[k] = -(wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF)) + (wm[k]*(-CF));
    // pent_c[k] = -(wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF));
    // pent_d[k] = -(wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF));
    // pent_e[k] = -(-wp[k+1]*(-CF));
  //}
    k=1;
    pent_b[k]+=pent_a[k];
    pent_a[k]=0;

  } else if(BCtop == 2) {
    // k=0; //top, no flux

    //uniform
    k=0; //top    
    pent_c[k]-=3*pent_a[k] + pent_b[k]; 
    pent_a[k]=0;
    pent_b[k]=0;         

    //if(Nk>1){
    k=1; //one below top   
    pent_b[k]-=pent_a[k];
    pent_a[k]=0;                                                                                                                                          
    // pent_a[k] = 0;
    // pent_b[k] = -(wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF)) - (wm[k]*(-CF));
    // pent_c[k] = -(wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF));
    // pent_d[k] = -(wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF));
    // pent_e[k] = -(-wp[k+1]*(-CF));
    //}
  } 

//note BC options are 1=no flux, 2=zero at boundary, currently no option for true no slip (speed at boundary is speed of boundary)

//if(Nk>1){
if(BCbot ==1) {

 //uniform 
 //no flux, so phi_Nk = phi_(Nk+1) = phi_(Nk-1)
 k=Nk-2; //one above bottom
 pent_d[k]+=pent_e[k];
 pent_e[k]=0;

  //  pent_a[k] = -(wm[k]*(-CF));
  //  pent_b[k] = -(wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
  //  pent_c[k] = -(wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF));
  //  pent_d[k] = -(wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF)) + (-wp[k+1]*(-CF));
  //  pent_e[k] = 0;

  k = Nk-1; //bottom
  pent_a[k] = -(wm[k]*(-CF));
  pent_b[k] = -(wp[k]*(.5-CF) + wm[k]*(.5+2*CF));
  pent_c[k] = -(wp[k]*(.5+2*CF) + wm[k]*(.5-CF));
  pent_d[k] = -(wp[k]*(-CF));
  pent_e[k] = 0;

  pent_c[k]+=pent_d[k];
  pent_d[k]=0;
  //  pent_a[k] = -(wm[k]*(-CF));
  //  pent_b[k] = -(wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
  //  pent_c[k] = -(wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF)) + (-wp[k+1]*(-CF)) + (wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF));
  //  pent_d[k] = 0;
  //  pent_e[k] = 0;
} else if(BCbot==2){

 //uniform
 //no slip, so phi_Nk = 2*phi_(Nk+1) = -phi_(Nk-1) 
 k=Nk-2; //one above bottom     
 pent_d[k]-=pent_e[k];
 pent_e[k]=0;                                                                                                                                  
//  pent_a[k] = -(wm[k]*(-CF));
//  pent_b[k] = -(wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
//  pent_c[k] = -(wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF));
//  pent_d[k] = -(wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF)) - (-wp[k+1]*(-CF));
//  pent_e[k] = 0;

 k = Nk-1; //bottom      
 pent_a[k] = -(wm[k]*(-CF));
 pent_b[k] = -(wp[k]*(.5-CF) + wm[k]*(.5+2*CF));
 pent_c[k] = -(wp[k]*(.5+2*CF) + wm[k]*(.5-CF));
 pent_d[k] = -(wp[k]*(-CF));
 pent_e[k] = 0;

 pent_c[k]-=pent_d[k];
 pent_d[k]=0;   
 

}
//}

}


/*
 * Function: GetPentDiag
 * Usage:  GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
 *	           phys->w,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
 * ------------------------------------------------------------------------------
 * Calculate the coefficients of the pentadiagonal for implicit part of vertical advection for QUICK or SHARP
 * W is the vertical velocity being advected (on vertical edges) and w is the effective omega doing the advection (at centers)
 *
 */
void GetPentDiagVert(REAL *pent_a, REAL *pent_b, REAL *pent_c, REAL *pent_d, REAL *pent_e, REAL *wp, REAL *wm,
  REAL *w, REAL **dzz, REAL **W, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot) {
int k;

REAL CF, phi_c, phi_c1, CF1, Bu, Bd, B, hkm2, hkm1, hk, hkp1, hkp2;

CF = 1/8;
CF1 = 1/8;

//w here is effective omega at centers. 0 is top cell, Nk-1 is bottom cell. 
for(k=0;k<Nk;k++) {
   wp[k] = 0.5*(w[k]+fabs(w[k]));
   wm[k] = 0.5*(w[k]-fabs(w[k]));
}
//w[k] = w at k-1/2, w[k+1] = w at k-1/2

for(k=1;k<Nk-1;k++) { //loop over k=2 to k=Nk-2, which should be the "normal" rows of the pentdiagonal
 
  //Just do QUICK, shouldn't need SHARP for momentum anyways 
 /*if(TVD==6){
   if(wp[k]>0){
     if(W[i][k-1]!W[i][k+1]){
       phi_c = (W[i][k] - W[i][k+1])/(W[i][k-1]-W[i][k+1]); //hopefully this doesn't change in nonuniform....
     } else{
       phi_c = 1.6; //this will make CF=1/8
     }
   } else if(wm[k] != 0) {
     if(U[j][k]!=U[j][k-2]) {
       phi_c = (W[i][k-1] - W[j][k-2])/(U[j][k]-U[j][k-2]);
     } else{
       phi_c = 1.6;
     }
   } else {
     phi_c = 1.6;
   }

   if(wp[k+1]>0){
     if((U[j][k]!=U[j][k+2])){
       phi_c1 = (U[j][k+1] - U[j][k+2])/(U[j][k]-U[j][k+2]);
     } else{
       phi_c1 = 1.6;
     }
   } else if(wm[k] != 0) {
     if(U[j][k+1]!=U[j][k-1]){
       phi_c1 = (U[j][k] - U[j][k-1])/(U[j][k+1]-U[j][k-1]);
     } else{
       phi_c1 = 1.6;
     }
   } else {
     phi_c1 = 1.6;
   }

   CF = GetCF(phi_c); //k
   CF1 = GetCF(phi_c1); //k+1
 } else {
   CF = 1/8;
   CF1 = 1/8;
 }
 */ 
 //uniform (remember these are on the FACES)
 //if k=0, then pent[k] is at top FACE and w[k] is top CELL. so for pent_(k-1/2) = pent[k] and w[k-1] = w[k-1]
 pent_a[k] = -(wm[k-1]*(-CF));
 pent_b[k] = -(wp[k-1]*(.5-CF) + wm[k-1]*(.5+2*CF) - wm[k]*(-CF1));
 pent_c[k] = -(wp[k-1]*(.5+2*CF) + wm[k-1]*(.5-CF) - wp[k]*(.5-CF1) - wm[k]*(.5+2*CF1));
 pent_d[k] = -(wp[k-1]*(-CF) - wp[k]*(.5+2*CF1) - wm[k]*(.5-CF1));
 pent_e[k] = -(-wp[k]*(-CF1));
}

//BC=1 (no flux), BC=2 (zero at boundary), BC=3 (no slip at boundary), 

//use QUICK at BCS
CF=1/8;

//deal with very top in main function since it uses different omegas. 

if(BCtop == 1 || BCtop == 3) {

 k=1; //one below top
//  pent_a[k] = 0;
//  pent_b[k] = -(wp[k-1]*(.5-CF) + wm[k-1]*(.5+2*CF) - wm[k]*(-CF)) + (wm[k-1]*(-CF));
//  pent_c[k] = -(wp[k-1]*(.5+2*CF) + wm[k-1]*(.5-CF) - wp[k]*(.5-CF) - wm[k]*(.5+2*CF));
//  pent_d[k] = -(wp[k-1]*(-CF) - wp[k]*(.5+2*CF) - wm[k]*(.5-CF));
//  pent_e[k] = -(-wp[k]*(-CF));

pent_b[k]+=pent_a[k];
pent_a[k]=0;
} else if(BCtop == 2) {
 // k=0; //top, no flux

 //uniform

 k=1; //one below top                                                                                                                                             
//  pent_a[k] = 0;
//  pent_b[k] = -(wp[k-1]*(.5-CF) + wm[k-1]*(.5+2*CF) - wm[k]*(-CF)) - (wm[k-1]*(-CF));
//  pent_c[k] = -(wp[k-1]*(.5+2*CF) + wm[k-1]*(.5-CF) - wp[k]*(.5-CF) - wm[k]*(.5+2*CF));
//  pent_d[k] = -(wp[k-1]*(-CF) - wp[k]*(.5+2*CF) - wm[k]*(.5-CF));
//  pent_e[k] = -(-wp[k]*(-CF));
pent_b[k]-=pent_a[k];
pent_a[k]=0;
 
} 

//note BC options are 1=no flux, 2=zero at boundary, currently no option for true no slip (speed at boundary is speed of boundary)
if(Nk>1){
  k=Nk-1;
  pent_a[k] = -(wm[k-1]*(-CF));
  pent_b[k] = -(wp[k-1]*(.5-CF) + wm[k-1]*(.5+2*CF));
  pent_c[k] = -(wp[k-1]*(.5+2*CF) + wm[k-1]*(.5-CF) - wp[k]*(.5) - wm[k]*(.5));
  pent_d[k] = 0;
  pent_e[k] = 0;

  k=Nk-2;
  pent_d[k]+=pent_e[k];
  pent_e[k]=0;
  
  // if(BCbot ==1) {

  // //uniform 
  // //no flux, so phi_Nk = phi_(Nk+1) = phi_(Nk-1)

  // k=Nk-1;
  // pent_a[k] = -(wm[k-1]*(-CF));
  // pent_b[k] = -(wp[k-1]*(.5-CF) + wm[k-1]*(.5+2*CF));
  // pent_c[k] = -(wp[k-1]*(.5+2*CF) + wm[k-1]*(.5-CF) - wp[k]*(.5) - wm[k]*(.5));
  // pent_d[k] = 0;
  // pent_e[k] = 0;

  // // k = Nk-1; //bottom
  // // pent_a[k] = -(wm[k-1]*(-CF));
  // // pent_b[k] = -(wp[k-1]*(.5-CF) + wm[k-1]*(.5+2*CF) - wm[k]*(-CF));
  // // pent_c[k] = -(wp[k-1]*(.5+2*CF) + wm[k-1]*(.5-CF) - wp[k]*(.5-CF) - wm[k]*(.5+2*CF)) + (-wp[k]*(-CF)) + (wp[k-1]*(-CF) - wp[k]*(.5+2*CF) - wm[k]*(.5-CF));
  // // pent_d[k] = 0;
  // // pent_e[k] = 0;
  // } else if(BCbot==2){

  // //uniform
  // //no slip, so phi_Nk = 2*phi_(Nk+1) = -phi_(Nk-1) 

  // k = Nk-1; //bottom                                                                                                                                               
  // pent_a[k] = -(wm[k-1]*(-CF));
  // pent_b[k] = -(wp[k-1]*(.5-CF) + wm[k-1]*(.5+2*CF) - wm[k]*(-CF));
  // pent_c[k] = -(wp[k-1]*(.5+2*CF) + wm[k-1]*(.5-CF) - wp[k]*(.5-CF) - wm[k]*(.5+2*CF)) - 3*(-wp[k]*(-CF)) - (wp[k-1]*(-CF) - wp[k]*(.5+2*CF) -  wm[k]*(.5-CF));
  // pent_d[k] = 0;
  // pent_e[k] = 0;

  // }
}


}


/*
 * Function: GetPentDiag
 * Usage:  GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
 *	           phys->w,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
 * ------------------------------------------------------------------------------
 * Calculate the coefficients of the pentadiagonal for implicit part of vertical advection for QUICK or SHARP
 *
 */
void GetPentDiag(REAL *pent_a, REAL *pent_b, REAL *pent_c, REAL *pent_d, REAL *pent_e, REAL *wp, REAL *wm,
		 REAL **w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot) {
  int k;

  REAL CF, phi_c, phi_c1, CF1, Bu, Bd, B, hkm2, hkm1, hk, hkp1, hkp2;

  CF = 1/8;
  CF1 = 1/8;


  for(k=0;k<Nk+1;k++) {
      wp[k] = 0.5*(w[i][k]+fabs(w[i][k]));
      wm[k] = 0.5*(w[i][k]-fabs(w[i][k]));
  }
  //w[k] = w at k-1/2, w[k+1] = w at k-1/2

  for(k=2;k<Nk-2;k++) { //loop over k=2 to k=Nk-2, which should be the "normal" rows of the pentdiagonal
    
    if(TVD==6){
      if(wp[k]>0){
	      if(scal[i][k-1]!=scal[i][k+1]){
	        phi_c = (scal[i][k] - scal[i][k+1])/(scal[i][k-1]-scal[i][k+1]); //hopefully this doesn't change in nonuniform....
	      } else{
	        phi_c = 1.6; //this will make CF=1/8
	      }
      } else if(wm[k] != 0) {
        if(scal[i][k]!=scal[i][k-2]) {
          phi_c = (scal[i][k-1] - scal[i][k-2])/(scal[i][k]-scal[i][k-2]);
        } else{
          phi_c = 1.6;
        }
      } else {
        phi_c = 1.6;
      }

      if(wp[k+1]>0){
	      if((scal[i][k]!=scal[i][k+2])){
	        phi_c1 = (scal[i][k+1] - scal[i][k+2])/(scal[i][k]-scal[i][k+2]);
	      } else{
	        phi_c1 = 1.6;
	      }
      } else if(wm[k] != 0) {
	      if(scal[i][k+1]!=scal[i][k-1]){
	        phi_c1 = (scal[i][k] - scal[i][k-1])/(scal[i][k+1]-scal[i][k-1]);
	      } else{
	        phi_c1 = 1.6;
	      }
      } else {
	      phi_c1 = 1.6;
      }

      CF = GetCF(phi_c); //k
      CF1 = GetCF(phi_c1); //k+1
    } else {
      CF = 1/8;
      CF1 = 1/8;
    }
    //uniform 
    pent_a[k] = (wm[k]*(-CF));
    pent_b[k] = (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF1));
    pent_c[k] = (wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF1) - wm[k+1]*(.5+2*CF1));
    pent_d[k] = (wp[k]*(-CF) - wp[k+1]*(.5+2*CF1) - wm[k+1]*(.5-CF1));
    pent_e[k] = (-wp[k+1]*(-CF1));

    //nonuniform 

    //calculate coefficients
    // Bu = CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2);
    // Bd = CF1/4*pow((dzz[i][k]+dzz[i][k+1]), 2);
    // hkm2 = dzz[i][k-2];
    // hkm1 = dzz[i][k-1];
    // hk = dzz[i][k];
    // hkp1 = dzz[i][k+1];
    // hkp2 = dzz[i][k+2];

    // pent_a[k] = wm[k]*(-Bu/hkm1*(1/(.5*(hkm2+hkm1))));
    // pent_b[k] = wp[k]*(.5 - Bu/hk*(1/(.5*(hkm1+hk)))) + wm[k]*(.5 + Bu/hkm1*(1/(.5*(hkm2+hkm1))+1/(.5*(hkm1+hk)))) - wm[k+1]*(-Bd/hk*(1/(.5*(hkm1+hk))));
    // pent_c[k] = wp[k]*(.5 + Bu/hk*(1/(.5*(hk+hkp1)) + 1/(.5*(hkm1+hk)))) + wm[k]*(.5 - Bu/hkm1*(1/(.5*(hkm1+hk)))) - wp[k+1]*(.5 - Bd/hkp1*(1/(.5*(hk+hkp1)))) - wm[k+1]*(.5 + Bd/hk*(1/(.5*(hkm1+hk)) + 1/(.5*hk+hkp1)));
    // pent_d[k] = wp[k]*(-Bu/hk*(1/(.5*(hk+hkp1)))) - wp[k+1]*(.5 + Bd/hkp1*(1/(.5*(hk+hkp1)) + 1/(.5*(hkp1+hkp2)))) - wm[k+1]*(.5 - Bd/hk*(1/(.5*(hk+hkp1))));
    // pent_e[k] = -wp[k+1]*(-Bd/hkp1*(1/(.5*(hkp1+hkp2))));

    //old
    // pent_a[k] = (wp[k]*(-Bu/dzz[i][k-1]*(1/(.5*(dzz[i][k-2]+dzz[i][k-1])))) + wm[k]*0) - (wp[k+1]*0 + wm[k+1]*0);
    // pent_b[k] = (wp[k]*(.5 - Bu/dzz[i][k]*(1/(.5*(dzz[i][k-1]+dzz[i][k])))) + wm[k]*(.5 + Bu/dzz[i][k-1]*(1/(.5*(dzz[i][k-2]+dzz[i][k-1])) + 1/(.5*(dzz[i][k-1]+dzz[i][k]))))) 
    //             - (wp[k+1]*0 + wm[k+1]*(Bd/dzz[i][k]*(1/(.5*(dzz[i][k-1]+dzz[i][k])))));
    // pent_c[k] = (wp[k]*(.5 + Bu/dzz[i][k]*(1/(.5*(dzz[i][k]+dzz[i][k+1])) + 1/(.5*(dzz[i][k-1]+dzz[i][k])))) + wm[k]*(.5 - Bu/dzz[i][k-1]*(1/(.5*(dzz[i][k-1]+dzz[i][k]))))) 
    // - (wp[k+1]*(.5 - Bd/dzz[i][k+1]*(1/(.5*(dzz[i][k]+dzz[i][k+1])))) + wm[k+1]*(.5 + Bd/dzz[i][k]*(1/(.5*(dzz[i][k-1]+dzz[i][k])) + 1/(.5*(dzz[i][k]+dzz[i][k+1])))));
    // pent_d[k] = (wp[k]*(-Bu/dzz[i][k]*(1/(.5*(dzz[i][k]+dzz[i][k+1])))) + wm[k]*0) 
    // - (wp[k+1]*(.5 + Bd/dzz[i][k+1]*(1/(.5*(dzz[i][k]+dzz[i][k+1])) + 1/(.5*(dzz[i][k+1]+dzz[i][k+2])))) + wm[k+1]*(.5 - Bd/dzz[i][k]*(1/(.5*(dzz[i][k]+dzz[i][k+1])))));
    // pent_e[k] = (wp[k]*0 + wm[k]*0) - (wp[k+1]*(-Bd/dzz[i][k+1]*(1/(.5*(dzz[i][k+1]+dzz[i][k+2])))) + wm[k+1]*0);

    }

  //BC=1 (no flux), BC=2 (zero at boundary), BC=3 (no slip at boundary), 

  //use QUICK at BCS
  CF=1/8;

  //BOUNDARY CONDITION SETUP FOR NONUNIFORM GRID

  //set up for top
  /*
  k = 0;
  Bu = CF/4*pow((dzz[i][k]+dzz[i][k]), 2);
  Bd = CF/4*pow((dzz[i][k]+dzz[i][k+1]), 2);

  hkm2 = dzz[i][k];
  hkm1 = dzz[i][k];
  hk = dzz[i][k];
  hkp1 = dzz[i][k+1];
  hkp2 = dzz[i][k+2];

  pent_a[k] = wm[k]*(-Bu/hkm1*(1/(.5*(hkm2+hkm1))));
  pent_b[k] = wp[k]*(.5 - Bu/hk*(1/(.5*(hkm1+hk)))) + wm[k]*(.5 + Bu/hkm1*(1/(.5*(hkm2+hkm1))+1/(.5*(hkm1+hk)))) - wm[k+1]*(-Bd/hk*(1/(.5*(hkm1+hk))));
  pent_c[k] = wp[k]*(.5 + Bu/hk*(1/(.5*(hk+hkp1)) + 1/(.5*(hkm1+hk)))) + wm[k]*(.5 - Bu/hkm1*(1/(.5*(hkm1+hk)))) - wp[k+1]*(.5 - Bd/hkp1*(1/(.5*(hk+hkp1)))) - wm[k+1]*(.5 + Bd/hk*(1/(.5*(hkm1+hk)) + 1/(.5*hk+hkp1)));
  pent_d[k] = wp[k]*(-Bu/hk*(1/(.5*(hk+hkp1)))) - wp[k+1]*(.5 + Bd/hkp1*(1/(.5*(hk+hkp1)) + 1/(.5*(hkp1+hkp2)))) - wm[k+1]*(.5 - Bd/hk*(1/(.5*(hk+hkp1))));
  pent_e[k] = -wp[k+1]*(-Bd/hkp1*(1/(.5*(hkp1+hkp2))));

  //set up for one below top
  k = 1;
  Bu = CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2);
  Bd = CF/4*pow((dzz[i][k]+dzz[i][k+1]), 2);

  hkm2 = dzz[i][k-1];
  hkm1 = dzz[i][k-1];
  hk = dzz[i][k];
  hkp1 = dzz[i][k+1];
  hkp2 = dzz[i][k+2];

  pent_a[k] = wm[k]*(-Bu/hkm1*(1/(.5*(hkm2+hkm1))));
  pent_b[k] = wp[k]*(.5 - Bu/hk*(1/(.5*(hkm1+hk)))) + wm[k]*(.5 + Bu/hkm1*(1/(.5*(hkm2+hkm1))+1/(.5*(hkm1+hk)))) - wm[k+1]*(-Bd/hk*(1/(.5*(hkm1+hk))));
  pent_c[k] = wp[k]*(.5 + Bu/hk*(1/(.5*(hk+hkp1)) + 1/(.5*(hkm1+hk)))) + wm[k]*(.5 - Bu/hkm1*(1/(.5*(hkm1+hk)))) - wp[k+1]*(.5 - Bd/hkp1*(1/(.5*(hk+hkp1)))) - wm[k+1]*(.5 + Bd/hk*(1/(.5*(hkm1+hk)) + 1/(.5*hk+hkp1)));
  pent_d[k] = wp[k]*(-Bu/hk*(1/(.5*(hk+hkp1)))) - wp[k+1]*(.5 + Bd/hkp1*(1/(.5*(hk+hkp1)) + 1/(.5*(hkp1+hkp2)))) - wm[k+1]*(.5 - Bd/hk*(1/(.5*(hk+hkp1))));
  pent_e[k] = -wp[k+1]*(-Bd/hkp1*(1/(.5*(hkp1+hkp2))));

  //set up for one above bottom
  k = Nk-2;
  Bu = CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2);
  Bd = CF/4*pow((dzz[i][k]+dzz[i][k+1]), 2);

  hkm2 = dzz[i][k-2];
  hkm1 = dzz[i][k-1];
  hk = dzz[i][k];
  hkp1 = dzz[i][k+1];
  hkp2 = dzz[i][k+1];

  pent_a[k] = wm[k]*(-Bu/hkm1*(1/(.5*(hkm2+hkm1))));
  pent_b[k] = wp[k]*(.5 - Bu/hk*(1/(.5*(hkm1+hk)))) + wm[k]*(.5 + Bu/hkm1*(1/(.5*(hkm2+hkm1))+1/(.5*(hkm1+hk)))) - wm[k+1]*(-Bd/hk*(1/(.5*(hkm1+hk))));
  pent_c[k] = wp[k]*(.5 + Bu/hk*(1/(.5*(hk+hkp1)) + 1/(.5*(hkm1+hk)))) + wm[k]*(.5 - Bu/hkm1*(1/(.5*(hkm1+hk)))) - wp[k+1]*(.5 - Bd/hkp1*(1/(.5*(hk+hkp1)))) - wm[k+1]*(.5 + Bd/hk*(1/(.5*(hkm1+hk)) + 1/(.5*hk+hkp1)));
  pent_d[k] = wp[k]*(-Bu/hk*(1/(.5*(hk+hkp1)))) - wp[k+1]*(.5 + Bd/hkp1*(1/(.5*(hk+hkp1)) + 1/(.5*(hkp1+hkp2)))) - wm[k+1]*(.5 - Bd/hk*(1/(.5*(hk+hkp1))));
  pent_e[k] = -wp[k+1]*(-Bd/hkp1*(1/(.5*(hkp1+hkp2))));

  //set up for bottom
  k = Nk-1;
  Bu = CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2);
  Bd = CF/4*pow((dzz[i][k]+dzz[i][k]), 2);

  hkm2 = dzz[i][k-2];
  hkm1 = dzz[i][k-1];
  hk = dzz[i][k];
  hkp1 = dzz[i][k];
  hkp2 = dzz[i][k];

  pent_a[k] = wm[k]*(-Bu/hkm1*(1/(.5*(hkm2+hkm1))));
  pent_b[k] = wp[k]*(.5 - Bu/hk*(1/(.5*(hkm1+hk)))) + wm[k]*(.5 + Bu/hkm1*(1/(.5*(hkm2+hkm1))+1/(.5*(hkm1+hk)))) - wm[k+1]*(-Bd/hk*(1/(.5*(hkm1+hk))));
  pent_c[k] = wp[k]*(.5 + Bu/hk*(1/(.5*(hk+hkp1)) + 1/(.5*(hkm1+hk)))) + wm[k]*(.5 - Bu/hkm1*(1/(.5*(hkm1+hk)))) - wp[k+1]*(.5 - Bd/hkp1*(1/(.5*(hk+hkp1)))) - wm[k+1]*(.5 + Bd/hk*(1/(.5*(hkm1+hk)) + 1/(.5*hk+hkp1)));
  pent_d[k] = wp[k]*(-Bu/hk*(1/(.5*(hk+hkp1)))) - wp[k+1]*(.5 + Bd/hkp1*(1/(.5*(hk+hkp1)) + 1/(.5*(hkp1+hkp2)))) - wm[k+1]*(.5 - Bd/hk*(1/(.5*(hk+hkp1))));
  pent_e[k] = -wp[k+1]*(-Bd/hkp1*(1/(.5*(hkp1+hkp2))));
 */

  if(BCtop == 1 || BCtop == 3) {
    // k=0; //top, no flux

    // //add on what would be a and b
    // pent_c[k]+=pent_a[k] + pent_b[k];
    // pent_a[k] = 0;
    // pent_b[k] = 0;

    // k=1; //one below top 

    // //add on what would be 
    // pent_b[k]+=pent_a[k];
    // pent_a[k] = 0;

    //uniform 
    k = 0;
    pent_a[k] = 0;
    pent_b[k] = 0;
    pent_c[k] = (wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF)) + (wm[k]*(-CF)) + (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
    pent_d[k] = (wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF));
    pent_e[k] = (-wp[k+1]*(-CF));

    k=1; //one below top
    pent_a[k] = 0;
    pent_b[k] = (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF)) + (wm[k]*(-CF));
    pent_c[k] = (wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF));
    pent_d[k] = (wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF));
    pent_e[k] = (-wp[k+1]*(-CF));
  } else if(BCtop == 2) {
    // k=0; //top, no flux

    // //subtract what would be a and b
    // pent_c[k]-=(3*pent_a[k] + pent_b[k]);
    // pent_a[k] = 0;
    // pent_b[k] = 0;

    // k=1; //one below top 

    // //add on what would be 
    // pent_b[k]-=pent_a[k];
    // pent_a[k] = 0;

    //uniform
    k=0; //top                                                                                                                                                       
    pent_a[k] = 0;
    pent_b[k] = 0;
    pent_c[k] = (wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF)) - 3*(wm[k]*(-CF)) - (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
    pent_d[k] = (wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF));
    pent_e[k] = (-wp[k+1]*(-CF));

    k=1; //one below top                                                                                                                                             
    pent_a[k] = 0;
    pent_b[k] = (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF)) - (wm[k]*(-CF));
    pent_c[k] = (wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF));
    pent_d[k] = (wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF));
    pent_e[k] = (-wp[k+1]*(-CF));
    
  } 

  //note BC options are 1=no flux, 2=zero at boundary, currently no option for true no slip (speed at boundary is speed of boundary)

  if(BCbot ==1) {
    // k=Nk-2; //one above bottom, no flux

    // //add on what would be e and e
    // pent_c[k]+=pent_d[k] + pent_e[k];
    // pent_d[k] = 0;
    // pent_e[k] = 0;

    // k=Nk-1; //bottom, no flux 

    // //add on what would be e
    // pent_d[k]+=pent_e[k];
    // pent_e[k] = 0;

    //uniform 
    //no flux, so phi_Nk = phi_(Nk+1) = phi_(Nk-1)
    k=Nk-2; //one above bottom
    pent_a[k] = (wm[k]*(-CF));
    pent_b[k] = (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
    pent_c[k] = (wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF));
    pent_d[k] = (wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF)) + (-wp[k+1]*(-CF));
    pent_e[k] = 0;
  
    k = Nk-1; //bottom
    pent_a[k] = (wm[k]*(-CF));
    pent_b[k] = (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
    pent_c[k] = (wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF)) + (-wp[k+1]*(-CF)) + (wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF));
    pent_d[k] = 0;
    pent_e[k] = 0;
  } else if(BCbot==2){
    // k=Nk-2; //one above bottom, no slip

    // //subtract what would be e and d
    // pent_c[k]-=(pent_d[k] + 3*pent_e[k]);
    // pent_d[k] = 0;
    // pent_e[k] = 0;

    // k=Nk-1; //bottom, no flux 

    // //subtract on what would be 
    // pent_d[k]-=pent_e[k];
    // pent_e[k] = 0;


    //uniform
    //no slip, so phi_Nk = 2*phi_(Nk+1) = -phi_(Nk-1) 
    k=Nk-2; //one above bottom                                                                                                                                       
    pent_a[k] = (wm[k]*(-CF));
    pent_b[k] = (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
    pent_c[k] = (wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF));
    pent_d[k] = (wp[k]*(-CF) - wp[k+1]*(.5+2*CF) - wm[k+1]*(.5-CF)) - (-wp[k+1]*(-CF));
    pent_e[k] = 0;

    k = Nk-1; //bottom                                                                                                                                               
    pent_a[k] = (wm[k]*(-CF));
    pent_b[k] = (wp[k]*(.5-CF) + wm[k]*(.5+2*CF) - wm[k+1]*(-CF));
    pent_c[k] = (wp[k]*(.5+2*CF) + wm[k]*(.5-CF) - wp[k+1]*(.5-CF) - wm[k+1]*(.5+2*CF)) - 3*(-wp[k+1]*(-CF)) - (wp[k]*(-CF) - wp[k+1]*(.5+2*CF) -  wm[k+1]*(.5-CF));
    pent_d[k] = 0;
    pent_e[k] = 0;

  }
  
}

/*
 * Function: GetScalarOnFace
 * Usage:  GetApAm(ap,am,phys->wp,phys->wm,phys->Cp,phys->Cm,phys->rp,phys->rm,
 *           phys->w,grid->dzz,scal,i,grid->Nk[i],ktop,prop->dt,prop->TVD);
 * ------------------------------------------------------------------------------
 * Calculate the coefficients of the pentadiagonal for implicit part of vertical advection for QUICK or SHARP
 *
 */
void GetScalarOnFace(REAL **SfHvf, REAL *wp, REAL *wm, REAL **w, REAL **dzz, REAL **scal, int i, int Nk, int ktop, REAL dt, int TVD, int BCtop, int BCbot, REAL **sum_neighs) {
  int k;

  REAL CF;
  REAL phi_c;
  CF=1/8; 

  // Implicit vertical advection terms
  // these all get multiplied by vertical velocity at the current time step, so just set to 1 or 0 depending on direction
  for(k=0;k<Nk+1;k++) {
    if(w[i][k]>0){
      wp[k]=1;
      wm[k]=0;
    } else if(w[i][k]<0){
      wm[k]=1;
      wp[k]=0;
    } else{
      wm[k]=0;
      wp[k]=0;
    }
  }

  //w[k] = w at k-1/2, w[k+1] = w at k-1/2. Same with SfHvf[i][k].
  for(k=2;k<Nk-1;k++) { //loop over k=2 to k=Nk-2, which should be the "normal" rows of the pentdiagonal
    if(TVD==6){
      if(wp[k]>0){
        if(scal[i][k-1]!=scal[i][k+1]){
          phi_c = (scal[i][k] - scal[i][k+1])/(scal[i][k-1]-scal[i][k+1]);
        } else {
          phi_c = 1.6;
        }
      } else if(wm[k] != 0) {
        if(scal[i][k]!=scal[i][k-2]){
          phi_c = (scal[i][k-1] - scal[i][k-2])/(scal[i][k]-scal[i][k-2]);
        } else {
          phi_c = 1.6;
        }
      } else {
        phi_c = 1.6;
      }
      
      CF = GetCF(phi_c);
    } else {
      CF = 1/8;
    }

    //uniform
    //SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-2]) 
    //  + 1/24*(wp[k]*sum_neighs[i][k] + wm[k]*sum_neighs[i][k-1]);

    //without transverse
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-2]);


    //nonuniform (without transverse term)
    //SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
    //             + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k])));
    
    //without transverse term
    //SfHvf[i][k] = .5*(scal[i][k]+scal[i][k+1])-CF*wp[k]*(scal[i][k]-2*scal[i][k+1]+scal[i][k+2])-CF*wm[k]*(scal[i][k+1]-2*scal[i][k]+scal[i][k-1]);
  }

  //Boundary conditions
  //use quick at boundaries 
  CF=1/8;

  //uniform 
  if(BCtop==1){ //no flux
    k=0; //top 
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k])-CF*wp[k]*(scal[i][k]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k]+scal[i][k]);

    k=1; //1 below top 
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-1]);
  } else{ //zero at boundary
    k=0; //top                                                                                                                             
    SfHvf[i][k] = .5*(scal[i][k]-scal[i][k])-CF*wp[k]*(-scal[i][k]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]+2*scal[i][k]-3*scal[i][k]);

    k=1; //1 below top                                                                                                                                              
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]-scal[i][k-1]);
  } 

  if(BCbot==1){ //no flux
    k=Nk-1; //bottom 
    //SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-2]);
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-2]);
    
    k=Nk;
    SfHvf[i][k] = .5*(scal[i][k-1]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k-1]+scal[i][k-1])-CF*wm[k]*(scal[i][k-1]-2*scal[i][k-1]+scal[i][k-2]);
  } else if(BCbot==2){ //zero at boundary
    k=Nk-1; //bottom
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]-scal[i][k])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-2]);
  
    k=Nk;
    SfHvf[i][k] = .5*(-scal[i][k-1]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]+2*scal[i][k-1]-3*scal[i][k-1])-CF*wm[k]*(-scal[i][k-1]-2*scal[i][k-1]+scal[i][k-2]);
  }

  //non uniform so far
  /*
  if(BCtop==1){ //no flux
    k=0; //top
    //SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-2]);
    // SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
    //              + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k])));
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k])-CF/4*pow((dzz[i][k]+dzz[i][k]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k]+dzz[i][k]))*(scal[i][k] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
                  + wm[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k]+dzz[i][k]))*(scal[i][k] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k]))*(scal[i][k] - scal[i][k])));

    k=1; //1 below top 
    //SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-2]);
    // SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
    //              + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k])));
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
                  + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-1]+dzz[i][k-1]))*(scal[i][k-1] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k])));

  } else{ //zero at boundary
    k=0; //top, scal[i][k-1] = -scal[i][k] and scal[i][k-1] = -3*scal[i][k] 
    //SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-2]);                                                                                                                   
    SfHvf[i][k] = .5*(scal[i][k]-scal[i][k])-CF/4*pow((dzz[i][k]+dzz[i][k]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k]+dzz[i][k]))*(-scal[i][k] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
    + wm[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k]+dzz[i][k]))*(-3*scal[i][k] + scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k]))*(-scal[i][k] - scal[i][k])));

    k=1; //1 below top, scal[i][k-2] = -scal[i][k-1]   
    //SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]+scal[i][k-2])                                                                                                                                 
    //SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF*wp[k]*(scal[i][k-1]-2*scal[i][k]+scal[i][k+1])-CF*wm[k]*(scal[i][k]-2*scal[i][k-1]-scal[i][k]);
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
    + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-1]+dzz[i][k-1]))*(-scal[i][k-1] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k])));

  } 

  if(BCbot==1){ //no flux
    k=Nk-1; //one above bottom, scal[i][k+2] = scal[i][k+1]
    // SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
    //              + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k])));
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k]))*(scal[i][k] - scal[i][k])) 
    + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k])));

    k=Nk; //bottom, scal[i][k+1] = scal[i][k+2] = scal[i][k]
    // SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
    //              + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k])));
    SfHvf[i][k] = .5*(scal[i][k-1]+scal[i][k-1])-CF/4*pow((dzz[i][k-1]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-1]+dzz[i][k-1]))*(scal[i][k-1] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k-1]))*(scal[i][k-1] - scal[i][k-1])) 
                  + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k-1]))*(scal[i][k-1] - scal[i][k-1])));; //enforce this 


  } else if(BCbot==2){ //zero at boundary
    k=Nk-1; //one above bottom, scal[i][k+1] = -scal[i][k], scal[i][k+2] = -3*scal[i][k]
    // SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k+1]))*(scal[i][k] - scal[i][k+1])) 
    //              + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k])));
    SfHvf[i][k] = .5*(scal[i][k]+scal[i][k-1])-CF/4*pow((dzz[i][k]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k])*(1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]) -  1/(.5*(dzz[i][k]+dzz[i][k]))*(scal[i][k] + scal[i][k])) 
                  + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k]))*(scal[i][k-1] - scal[i][k]))); 

    k=Nk;  //bottom, scal[i][k+1] = -scal[i][k], scal[i][k+2] = -3*scal[i][k]
    SfHvf[i][k] = .5*(-scal[i][k-1]+scal[i][k-1])-CF/4*pow((dzz[i][k-1]+dzz[i][k-1]), 2)*(wp[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-1]+dzz[i][k-1]))*(scal[i][k-1] + scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k-1]))*(-scal[i][k-1] + 3*scal[i][k-1])) 
                  + wm[k]*(1/dzz[i][k-1])*(1/(.5*(dzz[i][k-2]+dzz[i][k-1]))*(scal[i][k-2] - scal[i][k-1]) -  1/(.5*(dzz[i][k-1]+dzz[i][k-1]))*(scal[i][k-1] + scal[i][k-1])));
    //SfHvf[i][k] = 0;
  }
]*/

}


/*
 * Function: HorizontalFaceU
 * Usage: HorizontalFaceScalars(uc, grid, phys, boundary_scal);
 * ---------------------------------------------------------------------------
 * Calculate the horizontal face values for uc and vc using TVD schemes.  
 * SfHp[Ne][Nk] & SfHm[Ne][Nk] are used to store the facial values.
 */
void HorizontalFaceU(REAL **uc, gridT *grid, physT *phys, propT *prop, int TVD, 
			   MPI_Comm comm, int myproc) 
{
  int i, k, nf, iptr;
  int normal, nc1, nc2, ne;

  REAL df, dg, *sp, Ac, dt=prop->dt;
  REAL *Cp, *Cm, *rp, *rm, **gradSx, **gradSy, **stmp;   

  // For check!
  int iu, ku, nfu, gradflag, faceflag, Cflag, Rflag;

  // Pointer Variables for TVD scheme
  Cp = phys->Cp;
  Cm = phys->Cm;
  rp = phys->rp;
  rm = phys->rm;
  gradSx = phys->gradSx;
  gradSy = phys->gradSy;
  stmp = uc;


  for(i=0; i<grid->Nc; i++) {
    Ac = grid->Ac[i];
   
    // Initialize the gradSx and gradSy
    for(k=0;k<grid->Nk[i];k++){
      gradSx[i][k] = 0;
      gradSy[i][k] = 0;
    }

    // Loop through all faces of the current cell
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      df = grid->df[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) {
	nc2=nc1;
	//if(grid->mark[ne]>1) //(boundary_scal)
	//  sp=phys->stmp2[nc1];
	//else
	sp=stmp[nc1];
      }else 
        sp=stmp[nc2];


      for(k=0;k<grid->Nke[ne];k++) {
	gradSx[i][k]+=1/Ac*0.5*(stmp[nc1][k]+sp[k])*grid->n1[ne]*normal*df; 
	gradSy[i][k]+=1/Ac*0.5*(stmp[nc1][k]+sp[k])*grid->n2[ne]*normal*df;

      }

      for(k=grid->Nke[ne];k<grid->Nk[i];k++) {
	gradSx[i][k]+=1/Ac*stmp[i][k]*grid->n1[ne]*normal*df; 
	gradSy[i][k]+=1/Ac*stmp[i][k]*grid->n2[ne]*normal*df; 	
      }
    } 
  } 
  ISendRecvCellData3D(gradSx,grid,myproc,comm);  
  ISendRecvCellData3D(gradSy,grid,myproc,comm);  
 
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      df = grid->df[ne];
      dg = grid->dg[ne];
      nc1 = grid->grad[2*ne];
      nc2 = grid->grad[2*ne+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) {
	nc2=nc1;
	//if(boundary_scal)
	//  sp=phys->stmp2[nc1];
	//else
	  sp=stmp[nc1];
      } else 
	sp=stmp[nc2];


      for(k=0;k<grid->Nke[ne];k++) {
	Cp[k]= 0.5*(phys->u_old[ne][k]+fabs(phys->u_old[ne][k]))*dt/dg;
	Cm[k]= 0.5*(phys->u_old[ne][k]-fabs(phys->u_old[ne][k]))*dt/dg;
      }

      for(k=0;k<grid->Nke[ne];k++) {	
	rp[k]= 2*(gradSx[nc2][k]*grid->n1[ne]*dg + gradSy[nc2][k]*grid->n2[ne]*dg + EPS )/
	  (stmp[nc1][k]-sp[k]+EPS)-1;
	rm[k]= 2*(gradSx[nc1][k]*grid->n1[ne]*dg + gradSy[nc1][k]*grid->n2[ne]*dg + EPS )/
	  (stmp[nc1][k]-sp[k]+EPS)-1;
      }

      for(k=0;k<grid->Nke[ne];k++) {
	phys->SfHp[ne][k] = sp[k]+0.5*Psi(rp[k],TVD)*(1-Cp[k])*(stmp[nc1][k]-sp[k]);
	phys->SfHm[ne][k] = stmp[nc1][k]-0.5*Psi(rm[k],TVD)*(1+Cm[k])*(stmp[nc1][k]-sp[k]);
      }

      for(k=grid->Nke[ne];k<grid->Nk[nc1];k++) 
	phys->SfHm[ne][k] = stmp[nc1][k];

      for(k=grid->Nke[ne];k<grid->Nk[nc2];k++) 
	phys->SfHp[ne][k] = sp[k];
    } 
  }
}


/*                                                                                                                                                                   
 * Function: GetCF                                                                                                                                         
 * Usage:  GetCF(REAL phi_c);                                                                                          
 * ------------------------------------------------------------------------------                                                                                    
 * Calculate the coefficient for SHARP based on the normalized value of phi_c (see Leonard 1988)
 *                                                          
 *                                                                                                                                                                   
 */
REAL GetCF(REAL phi_c) {
  REAL CF = 1/8; //default 
    
  if(phi_c<=-1 || phi_c>=1.5){
    CF =  0.125;
  } else if(phi_c>-1 && phi_c<=0){
    CF = (0.5+0.125*phi_c)/(1-2*phi_c);
  } else if(phi_c>=1 && phi_c < 1.5) {
    CF = 0.5*(phi_c-1)/(2*phi_c-1);
  } else if(phi_c >= 0.3 && phi_c <= 0.7) {
    CF = 0.125 - 0.2609*(phi_c-1.5) + 0.13613*pow((phi_c-0.5), 2);
  } else {
    CF = (pow(phi_c, 2) - (1+phi_c)*(phi_c - 0.5) - sqrt(phi_c*pow((1-phi_c), 3)))/pow(1-2*phi_c, 2);
  }

  return CF;
  
}

/*
 * Function: SumNeighborScalars
 *
 *
 * Usage:  SumNeighborScalars(gridT *grid, physT *phys, REAL **scal, REAL **sum_neighs,
 *		MPI_Comm comm, int myproc)
 * ------------------------------------------------------------------------------                                                                                    
 * Calculate the sum of (scal[neigh][k]-scal[i][k]) for use in the transverse curvature term 
 * in vertical advection for 2D/3D QUICK and SHARP 
 *  
 */

void SumNeighborScalars(gridT *grid, physT *phys, REAL **scal, REAL **sum_neighs,
			MPI_Comm comm, int myproc) {
  int i, iptr, mf, k, ne, neigh;

  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      sum_neighs[i][k]=0;
      for(mf=0;mf<grid->nfaces[i];mf++) {
	      ne = grid->face[i*grid->maxfaces+mf];
	      neigh = grid->neigh[i*grid->maxfaces+mf];
	      if(neigh!=-1)
	        sum_neighs[i][k]+=scal[neigh][k]-scal[i][k];
      }
    }
  }
  ISendRecvCellData3D(sum_neighs,grid,myproc,comm);  
}
