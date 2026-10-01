/*
 * File: uservertcoordinate.c
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This file include a function to user defined vertical coordinate
 * 
 */

#include "suntans.h"
#include "grid.h"
#include "phys.h"
#include "initialization.h"
#include "boundaries.h"
#include "util.h"
#include "tvd.h"
#include "mympi.h"
#include "scalars.h"
#include "vertcoordinate.h"
#include "physio.h"
#include "subgrid.h"
#include "gridio.h"

REAL ZfromRho(REAL rho);

/*
 * Function: UserDefinedVerticalCoordinate
 * User define vertical coordinate 
 * basically it is a user-defined function to calculate the layer thickness based on 
 * different criterion
 * ----------------------------------------------------
 * the original code has already include 1 z-level, 2 isopycnal, 3 sigma, 4 variational 
 */
void UserDefinedVerticalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k;
  for(i=0;i<grid->Nc;i++)
  { 
    k=grid->ctop[i];
    grid->dzz[i][k]=7.0+phys->h[i];
    for(k=grid->ctop[i]+1;k<10;k++)
      grid->dzz[i][k]=7.0;
    for(k=10;k<50;k++)
      grid->dzz[i][k]=1.5;
    for(k=50;k<grid->Nk[i];k++)
      grid->dzz[i][k]=(grid->dv[i]-130.0)/(grid->Nk[i]-50);
  }
}

/*
 * Function: InitializeVerticalCoordinate
 * to setup the initial condition of dzz for user defined vertical coordinate
 * ----------------------------------------------------
 */
void InitializeVerticalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k;
  for(i=0;i<grid->Nc;i++)
  { 
    k=grid->ctop[i];
    grid->dzz[i][k]=7.0;
    for(k=grid->ctop[i]+1;k<10;k++)
      grid->dzz[i][k]=7.0;
    for(k=10;k<50;k++)
      grid->dzz[i][k]=1.5;
    for(k=50;k<grid->Nk[i];k++)
      grid->dzz[i][k]=(grid->dv[i]-130.0)/(grid->Nk[i]-50);
  }	
}

// Newton method to determine z from rho
REAL ZfromRho(REAL rho) {
  int iters=0, max_iters=100;
  REAL f, f_old, z_new, z, z_old, tol = 1e-3, beta = 1e-3, relax = 0.2;

  z_old = 0;
  z = z_old - 10;
  f_old = beta*ReturnSalinity(0,0,z_old) - rho;
  while(fabs(z-z_old)>tol && iters<max_iters) {
    f = beta*ReturnSalinity(0,0,z) - rho;

    if(f==f_old) break;

    z_new = z - relax*f*(z-z_old)/(f-f_old);

    f_old = f;
    z_old = z;
    z = z_new;

    iters++;
  }
  if(iters==max_iters)
    printf("Warning...iteration to find z in ZfromRho not converging after %d steps\n",iters);
  
  return z;
}

void InitializeIsopycnalCoordinate(gridT *grid, propT *prop, physT *phys, MPI_Comm comm, int myproc)
//void InitializeIsopycnalCoordinate(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i, k, k0, ktop, Nkmax=grid->Nkmax, Nk_mixed, Nk_pycnocline, Nk_lower;  
  REAL L0, a, h1, h2, zeta, eta, L, depth, mydmax, dmax, delta;
  REAL z, *zi, *zi0, *zptr, rhomin, rhomax, beta = 1e-3;

  Nk_mixed = 1;
  Nk_lower = 1;
  Nk_pycnocline = Nkmax-Nk_mixed-Nk_lower;

  mydmax=0;
  for(i=0;i<grid->Nc;i++) {
    if(grid->dv[i]>mydmax) mydmax=grid->dv[i];
  }
  MPI_Reduce(&mydmax,&(dmax),1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Bcast(&dmax,1,MPI_DOUBLE,0,comm);
  //  dmax=mydmax;
  //dmax=3518.385000;
  
  rhomin = beta*ReturnSalinity(0,0,0);
  rhomax = beta*ReturnSalinity(0,0,-dmax);
  REAL drho = (rhomax-rhomin)/((REAL)Nkmax);

  printf("myproc = %d, dmax = %f, drho = %f\n",myproc,dmax,drho);

  zi0 = (REAL *)SunMalloc((Nkmax+1)*sizeof(REAL),"UserVerticalCoordinate");
  zi = (REAL *)SunMalloc((Nkmax+1)*sizeof(REAL),"UserVerticalCoordinate");  
  REAL dz0=dmax/1000;
  REAL rho1, rho2, rho_i = rhomin;

  //  printf("z(%f) = %f\n",rho,ZfromRho(rho));

  /*
  z = 0;
  k=1;
  zi0[0] = 0;
  rho1 = beta*ReturnSalinity(0,0,0);
  rho_i = rho_i + drho;
  while(z>=-dmax) {
    rho2 = beta*ReturnSalinity(0,0,z-dz0);
    //    printf("z=%.2f k=%d rho2=%.2e rho1=%.2e rho_i=%.2e\n",z,k,rho2,rho1,rho_i);
    if(rho2>=rho_i && rho1<rho_i) {
      zi0[k]=z;
      rho_i+=drho;
      if(k==Nkmax)
	break;
      else
	k++;
    }
    rho1=rho2;
    z-=dz0;
  }
  zi0[Nkmax]=-dmax;
  */
  rho_i=0;
  for(k=1;k<Nkmax;k++) {
    rho_i+=drho;
    zi0[k]=ZfromRho(rho_i);
  }
  zi0[Nkmax]=-dmax;

  for(k=0;k<Nkmax;k++) {
    grid->dz[k]=zi0[k]-zi0[k+1];
    printf("proc %d: dz[%d]=%f, rho=%f\n",myproc,k,grid->dz[k],beta*ReturnSalinity(0,0,0.5*(zi0[k]+zi0[k+1])));
  }

  if(myproc==0) 
    WriteVertSpaceData(VERTSPACEFILE,grid,myproc);
  
  a = 1500;
  L0 = 6000;
  L = 200e3;
  for(i=0;i<grid->Nc;i++) {
    //    zeta = -a*exp(-pow((grid->xv[i]-L)/L0,2.0));    
    //    zeta = 1e-6*grid->xv[i];//-a*pow(cosh((grid->xv[i]-L)/L0),-2.0);
    zeta = -a*pow(cosh((grid->xv[i]-L)/L0),-2.0);    
    zi[0]=0;
    for(k=1;k<Nkmax;k++) {
      zi[k]=zi0[k]+zeta;//*sin(PI*zi0[k]/dmax);
    }
    zi[Nkmax]=-dmax;

    for(k=0;k<Nkmax;k++)
      grid->dzz[i][k]=zi[k]-zi[k+1];
  }

  // Adjust dz based on minimum layer height or depth
  for(i=0;i<grid->Nc;i++) {
    z=0;
    for(k=0;k<Nkmax;k++) {
      z-=grid->dzz[i][k];
      if(z<-grid->dv[i]) {
	grid->dzz[i][k]=z+grid->dzz[i][k]+grid->dv[i];

	if(grid->dzz[i][k]<vert->vertdzmin)
	  grid->dzz[i][k]=vert->vertdzmin;
	break;
      }
    }
    for(k0=k+1;k0<Nkmax;k0++)
      grid->dzz[i][k0]=vert->vertdzmin;
  }

  //  ISendRecvCellData3D(grid->dzz,grid,myproc,comm);  
  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      grid->dzzold[i][k]=grid->dzz[i][k];
  }
}

/*
 * Function: InitializeIsopycnalCoordinate
 * User define isopycnal coordinate 
 * define the initial dzz for each cell under isopycnal coordinate
 * ----------------------------------------------------
 */
void InitializeIsopycnalCoordinateOld(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i, k, k0, ktop, Nkmax=grid->Nkmax, Nk_mixed, Nk_pycnocline, Nk_lower;  
  REAL L0, a, h1, h2, zeta, eta, L, depth, dmax, delta;
  REAL z, *zptr;

  Nk_mixed = 1;
  Nk_lower = 1;
  Nk_pycnocline = Nkmax-Nk_mixed-Nk_lower;

  eta = 0;
  delta = 50;
  h1 = 200;
  L = 500e3;
  
  if(Nk_pycnocline==0)
    delta = 0;
  
  dmax=0;
  for(i=0;i<grid->Nc;i++) {
    if(grid->dv[i]>dmax) dmax=grid->dv[i];
  }
      
  for(i=-1;i<grid->Nc;i++) {
    if(i==-1) {
      zptr=grid->dz;
      zeta=0;
      ktop=0;
      depth=dmax;
    } else {
      zptr=grid->dzz[i];

      a = 100;
      L0 = 20000;
      zeta = -a*exp(-pow((grid->xv[i]-L)/L0,2.0));
      
      ktop=0;
      depth=grid->dv[i];
    }
    for(k=ktop;k<Nk_mixed;k++) 
      zptr[k]=(h1+eta-zeta-delta/2)/Nk_mixed;
    for(k=Nk_mixed;k<Nk_mixed+Nk_pycnocline;k++) 
      zptr[k]=delta/Nk_pycnocline;
    for(k=Nk_mixed+Nk_pycnocline;k<Nkmax;k++) 
      zptr[k]=(depth-h1+zeta-delta/2)/Nk_lower;

    zptr[0]=h1-delta/2-zeta;
    zptr[1]=delta;
    zptr[2]=depth-h1-delta/2+zeta;
  }

  // Adjust dz based on minimum layer height or depth
  for(i=0;i<grid->Nc;i++) {
    z=0;
    for(k=0;k<Nkmax;k++) {
      z-=grid->dzz[i][k];
      if(z<-grid->dv[i]) {
	grid->dzz[i][k]=z+grid->dzz[i][k]+grid->dv[i];

	if(grid->dzz[i][k]<vert->vertdzmin)
	  grid->dzz[i][k]=vert->vertdzmin;
	break;
      }
    }
    for(k0=k+1;k0<Nkmax;k0++)
      grid->dzz[i][k0]=vert->vertdzmin;
  }
  
  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      grid->dzzold[i][k]=grid->dzz[i][k];
  }
}

void OldInitializeIsopycnalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i, k, k0, ktop, Nkmax=grid->Nkmax, Nk_mixed, Nk_pycnocline, Nk_lower;  
  REAL L0=1200, a=50, h1=60, delta=60, eta, L=30e3;
  REAL z, *zptr;

  Nk_mixed = 2;
  Nk_pycnocline = 2;
  Nk_lower = Nkmax-Nk_mixed-Nk_pycnocline;

  for(i=-1;i<grid->Nc;i++) {
    if(i==-1) {
      zptr=grid->dz;
      eta=0;
      ktop=0;
    } else {
      zptr=grid->dzz[i];
      //eta=a*cos(PI*grid->xv[i]/10000);
      //      eta=-2*a*pow(1/cosh(grid->xv[i]/2/L0),2)
      //	-2*a*pow(1/cosh((grid->xv[i]-L)/2/L0),2);
      a=100;eta=-a*pow(1/cosh(grid->xv[i]/2/L0),2);

      ktop=grid->ctop[i];
    }

    for(k=ktop;k<Nk_mixed;k++) 
      zptr[k]=(h1-eta-delta/2)/Nk_mixed;
    for(k=Nk_mixed;k<Nk_mixed+Nk_pycnocline;k++) 
      zptr[k]=delta/Nk_pycnocline;
    for(k=Nk_mixed+Nk_pycnocline;k<Nkmax;k++) 
      zptr[k]=(0*grid->dv[i]+600-h1+eta-delta/2)/Nk_lower;
  }
  /*
  for(i=0;i<grid->Nc;i++) {
    k=grid->Nk[i]-1;
    grid->dzz[i][k]=grid->dzz[i][k]+10*cos(2*PI*grid->xv[i]/10000);
    k=grid->Nk[i]-2;
    grid->dzz[i][k]=grid->dzz[i][k]-10*cos(2*PI*grid->xv[i]/10000);
  }
  */
  for(i=0;i<grid->Nc;i++) {
    z=0;
    for(k=0;k<Nkmax;k++) {
      z-=grid->dzz[i][k];
      if(z<=-grid->dv[i]) {
	grid->dzz[i][k]=z+grid->dzz[i][k]+grid->dv[i];
	if(grid->dzz[i][k]<vert->vertdzmin)
	  grid->dzz[i][k]=vert->vertdzmin;
	break;
      }
    }
    for(k0=k+1;k0<Nkmax;k0++)
      grid->dzz[i][k0]=vert->vertdzmin;
  }

  for(i=0;i<grid->Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      grid->dzzold[i][k]=grid->dzz[i][k];
  }
}

/*
 * Function: InitializeVariationalCoordinate
 * Initialize dzz for variational vertical coordinate
 * --------zz--------------------------------------------
 */
void InitializeVariationalCoordinate(gridT *grid, propT *prop, physT *phys,int myproc)
{
  int i,k,sum=0;
  REAL depth, depth0 = 10;

  GetDZ(grid->dz,depth0,depth0,grid->Nkmax,myproc);
  
  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      depth = ReturnDepth(grid->xv[i],grid->yv[i]);
      grid->dzz[i][k]=grid->dz[k];//(depth+ReturnFreeSurface(grid->xv[i],grid->yv[i],0))/grid->Nkmax;
      grid->dzzold[i][k]=grid->dzz[i][k];
    }
  }

  if(myproc==0) 
    WriteVertSpaceData(VERTSPACEFILE,grid,myproc);  
}

/*
 * Function: UserDefinedSigmaCoordinate
 * User define sigma coordinate 
 * basically to define the dsigma for each layer
 * ----------------------------------------------------
 */
void InitializeSigmaCoordinate(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,k, status, sum=0;
  REAL depth = ReturnDepth(100,0);
  REAL r, *dz;
  REAL dsigsum = 0;

  //printf("in initialize sigma coordinate call \n");
  
  //allocate dz array                                                                                                                                       
  dz = (REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"UpdateBottomHeight");

  GetDZ(dz, depth, depth, grid->Nkmax, myproc);
  //printf("made it through dz call");

  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      //vert->dsigma[k]=1.0/grid->Nkmax;
      //printf("dz of k is %f for k=%d \n", dz[k], k);
      vert->dsigma[k]=dz[k]/depth;
      if(i==1){
	dsigsum += dz[k]/depth;
      }
      //printf("dv[i] is %f and h[i] is %f \n", grid->dv[i], phys->h[i]);
      grid->dzz[i][k]=vert->dsigma[k]*(grid->dv[i]+phys->h[i]);
      //grid->dzz[i][k]=vert->dsigma[k]*depth;
      grid->dzzold[i][k]=grid->dzz[i][k];
      // grid->dzzold2[i][k]=grid->dzz[i][k];
      // grid->dzzold3[i][k]=grid->dzz[i][k];
      if(i==0) {
	grid->dz[k]=vert->dsigma[k]*depth;
      }
    }
  }
  
  //printf("sum of dsigma = %f \n", dsigsum);
  
  if(myproc==0) {
    WriteVertSpaceData(VERTSPACEFILE,grid,myproc);}  
  
  SunFree(dz, grid->Nkmax*sizeof(REAL),"UpdateBottomHeight"); //free extra dz variable
}


/*
 * Function: MonitorFunctionForVariationalMethod
 * calculate the value of monitor function for the variational approach
 * to update layer thickness when nonlinear==4
 * ----------------------------------------------------
 * Mii=sqrt(1-alphaM*(drhodz)^2)
 */
void MonitorFunctionForAverageMethod(gridT *grid, propT *prop, physT *phys, int myproc)
{
   int i,k;
   REAL alphaM=0,minM=0.15,max;
   // nonlinear=1 or 5 stable with alpham=320
   // nonlinear=2 stable with alpham=60
   // nonlinear=4 stable with alpham=60
 
   for(i=0;i<grid->Nc;i++)
   {
     max=0;
     vert->Msum[i]=0;
     for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++){
       vert->Mc[i][k]=1000*(phys->rho[i][k-1]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k-1]+grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
       if(fabs(vert->Mc[i][k])>max)
         max=fabs(vert->Mc[i][k]);
     }
     
     // top boundary
     k=grid->ctop[i];
     vert->Mc[i][k]=1000*(phys->rho[i][k]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
     if(fabs(vert->Mc[i][k])>max)
       max=fabs(vert->Mc[i][k]);   
     // bottom boundary
     k=grid->Nk[i]-1;
     vert->Mc[i][k]=1000*(phys->rho[i][k-1]-phys->rho[i][k])/(0.5*grid->dzz[i][k-1]+0.5*grid->dzz[i][k]);
     if(fabs(vert->Mc[i][k])>max)
       max=fabs(vert->Mc[i][k]);   
     if(max<1)
       max=1;
     
     for(k=grid->ctop[i];k<grid->Nk[i];k++){ 
       vert->Mc[i][k]=sqrt(1+alphaM*vert->Mc[i][k]/max*vert->Mc[i][k]/max);
       if(vert->Mc[i][k]<minM)
         vert->Mc[i][k]=minM;     
       vert->Msum[i]+=1/vert->Mc[i][k];
     }
   }
}

/*
 * Function: MonitorFunctionForVariationalMethod
 * calculate the value of monitor function for the variational approach
 * to update layer thickness when nonlinear==4
 * solve the elliptic equation using iteration method
 * ----------------------------------------------------
 * Mii=sqrt(1-alphaM*(drhodz)^2)
 * alpha_H define how much horizontal diffusion 
 * alphaH define how much horizontal density gradient 
 * alphaV define how much vertical density gradient
 */
void MonitorFunctionForVariationalMethod(gridT *grid, propT *prop, physT *phys, int myproc, int numprocs, MPI_Comm comm)
{
  int i,k,j,nf,neigh,ne,kk,nc1,nc2;
  REAL normal,alpha_H=0, alphaH=0, alphaV=10, minM=0.15, maxM=100000,max,tmp;
  REAL max_gradient_v,max_gradient_h=0,max_gradient_h_global,H1,H2,rho1,rho2;
  // initialize everything zero
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      vert->Mc[i][k]=0;
  
  for(j=0;j<grid->Ne;j++)
    for(k=0;k<grid->Nke[j]+1;k++)
      vert->Me_l[j][k]=0;

  // calculate Mc first
  for(i=0;i<grid->Nc;i++)
  {
    // calculate gradient
    max_gradient_v=0;
    for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++){
      vert->Mc[i][k]=RHO0*(phys->rho[i][k-1]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k-1]+grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
      if(fabs(vert->Mc[i][k])>max_gradient_v)
        max_gradient_v=fabs(vert->Mc[i][k]);
    }

    // top boundary
    k=grid->ctop[i];
    vert->Mc[i][k]=RHO0*(phys->rho[i][k]-phys->rho[i][k+1])/(0.5*grid->dzz[i][k]+0.5*grid->dzz[i][k+1]);
    if(fabs(vert->Mc[i][k])>max_gradient_v)
      max_gradient_v=fabs(vert->Mc[i][k]);   

    // bottom boundary
    k=grid->Nk[i]-1;
    vert->Mc[i][k]=RHO0*(phys->rho[i][k-1]-phys->rho[i][k])/(0.5*grid->dzz[i][k-1]+0.5*grid->dzz[i][k]);
    if(fabs(vert->Mc[i][k])>max_gradient_v)
      max_gradient_v=fabs(vert->Mc[i][k]);   
    if(max_gradient_v<1)
      max_gradient_v=1;  
    // calculate monitor function value
    for(k=grid->ctop[i];k<grid->Nk[i];k++){ 
      if(alphaV!=0){
        if(vert->Mc[i][k]/max_gradient_v>(maxM-1)/sqrt(alphaV))
          vert->Mc[i][k]=maxM; 
        else
          vert->Mc[i][k]=sqrt(1+alphaV*vert->Mc[i][k]/max_gradient_v*vert->Mc[i][k]/max_gradient_v);
      } else 
        vert->Mc[i][k]=1;
    }
  }

  // calculate Me_l 
  // calculate gradient
  for(j=0;j<grid->Ne;j++)
  {
    nc1=grid->grad[2*j];
    nc2=grid->grad[2*j+1];
    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;   
    
    // interior layer
    for(k=grid->etop[j]+1;k<grid->Nke[j];k++)
    {
      rho1=grid->dzz[nc1][k-1]/(grid->dzz[nc1][k]+grid->dzz[nc1][k-1])*phys->rho[nc1][k]+
        grid->dzz[nc1][k]/(grid->dzz[nc1][k]+grid->dzz[nc1][k-1])*phys->rho[nc1][k-1];
      rho2=grid->dzz[nc2][k-1]/(grid->dzz[nc2][k]+grid->dzz[nc2][k-1])*phys->rho[nc2][k]+
        grid->dzz[nc2][k]/(grid->dzz[nc2][k]+grid->dzz[nc2][k-1])*phys->rho[nc2][k-1];
      vert->Me_l[j][k]=RHO0*(rho1-rho2)/grid->dg[j];
      if(fabs(vert->Me_l[j][k])>max_gradient_h)
        max_gradient_h=fabs(vert->Me_l[j][k]);
    }

    // top and bottom
    k=grid->etop[j];
    vert->Me_l[j][k]=RHO0*(phys->rho[nc1][k]-phys->rho[nc2][k])/grid->dg[j];
    if(fabs(vert->Me_l[j][k])>max_gradient_h)
      max_gradient_h=fabs(vert->Me_l[j][k]); 
    k=grid->Nke[j];
    vert->Me_l[j][k]=RHO0*(phys->rho[nc1][k-1]-phys->rho[nc2][k-1])/grid->dg[j];
    if(fabs(vert->Me_l[j][k])>max_gradient_h)
      max_gradient_h=fabs(vert->Me_l[j][k]);
  }

  // find max_global and normalize
  MPI_Reduce(&max_gradient_h,&max_gradient_h_global,1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Bcast(&max_gradient_h_global,1,MPI_DOUBLE,0,comm);
  if(max_gradient_h_global<1){
    max_gradient_h_global=1.0;
  }

  for(j=0;j<grid->Ne;j++)
    for(k=grid->etop[j];k<=grid->Nke[j];k++){
      vert->Me_l[j][k]=alpha_H*sqrt(1+alphaH*vert->Me_l[j][k]/max_gradient_h_global*
        vert->Me_l[j][k]/max_gradient_h_global);
    }
}
