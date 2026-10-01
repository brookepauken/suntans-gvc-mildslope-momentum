/*
 * File: vertcoordinate.c
 * Author: Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * This file include all the functions for the new general 
 * vertical coordinate
 * 
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
#include "uservertcoordinate.h"
#include "physio.h"
#include "subgrid.h"

void Generate_dzzolds_varationalvoid(gridT *grid, propT *prop, physT *phys);


/*
 * Function: AllocateVertCoordinate
 * Allocate space for the general vertical coordinate
 * ----------------------------------------------------
 *
 */
void AllocateandInitializeVertCoordinate(gridT *grid, propT *prop, int myproc)
{ 
  int i,j,k;

  // velocity field at layer center for each edge
  vert->uf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->vf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->wf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->qcf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omegaf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zf=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->Nkeb=(int *)SunMalloc(grid->Ne*sizeof(int),"AllocateVertCoordinate");
  vert->zfb=(REAL *)SunMalloc(grid->Ne*sizeof(REAL),"AllocateVertCoordinate");
  vert->modifydzf = MPI_GetValue(DATAFILE,"modifydzf","AllocateVertCoordinate",myproc);
  vert->dJdtmeth = MPI_GetValue(DATAFILE,"dJdtmeth","AllocateVertCoordinate",myproc);
  vert->dzfmeth = MPI_GetValue(DATAFILE,"dzfmeth","AllocateVertCoordinate",myproc);
  vert->thetaT = MPI_GetValue(DATAFILE,"thetaT","AllocateVertCoordinate",myproc);
  vert->vertdzmin = MPI_GetValue(DATAFILE,"vertdzmin","AllocateVertCoordinate",myproc);
  vert->Me_l=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->f_re=(REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"AllocateVertCoordinate");
  vert->CCNpe=(REAL *)SunMalloc(2*grid->Ne*sizeof(REAL),"AllocateVertCoordinate");
  for(j=0;j<grid->Ne;j++)
  {
    vert->Nkeb[j]=0;
    vert->zfb[j]=0.0;
    vert->CCNpe[2*j]=0;
    vert->CCNpe[2*j+1]=0;
    vert->qcf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->uf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->vf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->wf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->f_re[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->omegaf[j]=(REAL *)SunMalloc(grid->Nkc[j]*sizeof(REAL),"AllocateVertCoordinate");
    vert->zf[j]=(REAL *)SunMalloc((grid->Nkc[j]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->Me_l[j]=(REAL *)SunMalloc((grid->Nkc[j]+1)*sizeof(REAL),"AllocateVertCoordinate");
    for(k=0;k<grid->Nkc[j];k++)
    {
      vert->qcf[j][k]=0;
      vert->uf[j][k]=0;
      vert->zf[j][k]=0;
      vert->vf[j][k]=0;
      vert->wf[j][k]=0;
      vert->omegaf[j][k]=0;
      vert->Me_l[j][k]=0;
      vert->f_re[j][k]=0;
    }
    vert->Me_l[j][grid->Nkc[j]]=0;
  }

  // normal vector
  vert->n1=(REAL *)SunMalloc(grid->Nc*grid->maxfaces*sizeof(REAL),"AllocateVertCoordinate");
  vert->n2=(REAL *)SunMalloc(grid->Nc*grid->maxfaces*sizeof(REAL),"AllocateVertCoordinate");
  // velocity field at the top and bottom face of each layer
  vert->ul=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->vl=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega_im=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega_old=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omega_old2=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  // U3=w-udzdx-udzdy
  vert->U3=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->U3_old=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->U3_old2=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->omegac=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zc=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zl=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zcold=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->zlold=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");  
  vert->f_r=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dvdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dudy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dudx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dvdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dwdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dwdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dzdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dzdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dqdy=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dqdx=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->Mc=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->dzztmp=(REAL **)SunMalloc(grid->Nc*sizeof(REAL *),"AllocateVertCoordinate");
  vert->Msum=(REAL *)SunMalloc(grid->Nc*sizeof(REAL),"AllocateVertCoordinate");

  for(i=0;i<grid->Nc;i++)
  {
    vert->ul[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->vl[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega_im[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega_old[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omega_old2[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->U3[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->U3_old[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->U3_old2[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->omegac[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->zc[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->zcold[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->zl[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->zlold[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");    
    vert->dvdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dvdy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->f_r[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dudy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dudx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dwdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dwdy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dzdy[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->dzdx[i]=(REAL *)SunMalloc((grid->Nk[i]+1)*sizeof(REAL),"AllocateVertCoordinate");
    vert->dqdy[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dqdx[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->Mc[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    vert->dzztmp[i]=(REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"AllocateVertCoordinate");
    
    // initialize normal vector
    for(k=0;k<grid->maxfaces;k++)
    {
      vert->n1[i*grid->maxfaces+k]=0;
      vert->n2[i*grid->maxfaces+k]=0;
    }

    for(k=0;k<grid->Nk[i]+1;k++)
    {
      vert->ul[i][k]=0;
      vert->vl[i][k]=0;
      vert->omega[i][k]=0;
      vert->omega_im[i][k]=0;
      vert->omega_old[i][k]=0;
      vert->omega_old2[i][k]=0;
      vert->U3[i][k]=0;
      vert->U3_old[i][k]=0;
      vert->U3_old2[i][k]=0;
    }
    for(k=0;k<grid->Nk[i];k++)
    {
      vert->omegac[i][k]=0;
      vert->zc[i][k]=0;
      vert->zcold[i][k]=0;
      vert->zlold[i][k]=0;      
      vert->zl[i][k]=0;
      vert->f_r[i][k]=0;
      vert->dvdx[i][k]=0;
      vert->dudy[i][k]=0;
      vert->dudx[i][k]=0;
      vert->dvdy[i][k]=0;
      vert->dwdx[i][k]=0;
      vert->dwdy[i][k]=0;
      vert->dzdx[i][k]=0;      
      vert->dzdy[i][k]=0;
      vert->dqdx[i][k]=0;      
      vert->dqdy[i][k]=0;
      vert->Mc[i][k]=0;
      vert->dzztmp[i][k]=0;
    }
    vert->zl[i][grid->Nk[i]]=0;
    vert->Msum[i]=0;
  }

  // nodal
  vert->Ap=(REAL *)SunMalloc(grid->Np*sizeof(REAL),"AllocateVertCoordinate");
  vert->Nkpmin=(int *)SunMalloc(grid->Np*sizeof(int),"AllocateVertCoordinate");
  vert->typep=(int *)SunMalloc(grid->Np*sizeof(int),"AllocateVertCoordinate");
  vert->f_rp=(REAL **)SunMalloc(grid->Np*sizeof(REAL *),"AllocateVertCoordinate");
  for(i=0;i<grid->Np;i++)
  {
    vert->f_rp[i]=(REAL *)SunMalloc(grid->Nkp[i]*sizeof(REAL),"AllocateVertCoordinate"); 
    for(k=0;k<grid->Nkp[i];k++)
      vert->f_rp[i][k]=0.0;
    vert->Ap[i]=0.0;
    vert->Nkpmin[i]=grid->Nkp[i]; 
    vert->typep[i]=1;  
  }

  // temporary array for output 
  vert->tmp = (REAL *)SunMalloc(10*grid->Nc*sizeof(REAL),"AllocateVertCoordinate");
  vert->tmp_nc = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"AllocateVertCoordinate");

  // read vertical coordinate switch
  if(prop->vertcoord==3)
    vert->dsigma=(REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"AllocateVertCoordinate");

    int Mid_Stretch=1;
  if(prop->vertcoord==4 && Mid_Stretch==1)
    vert->dsigma=(REAL *)SunMalloc(grid->Nkmax*sizeof(REAL),"AllocateVertCoordinate");
}

/*
 * Function: InitializeLayerThickness
 * Calculate dz for the first time step based on the vertical coordinate chosen
 * ----------------------------------------------------
 * The switch is vertcoord 
 * 0 for user defined, 1 for z level, 2 for isopycnal,3 for sigma
 * need to modify grid.c to make sure ctop=0 and Nk=Nkmax for all cells 
 */
void InitializeLayerThickness(gridT *grid, propT *prop, physT *phys,MPI_Comm comm,int myproc)
{
  int i,j,k;
  if(prop->vertcoord!=5){
    for(i=0;i<grid->Nc;i++)
    {
      grid->ctopold[i]=grid->ctop[i]=0;
      if(grid->Nk[i]!=grid->Nkmax){
        printf("There is something wrong in setting the Nkmax for each cell, Nk=%d at i=%d (Nkmax=%d)\n",grid->Nk[i],i,grid->Nkmax);
        exit(1);
      }
    }

    for(j=0;j<grid->Ne;j++)
      grid->etopold[j]=grid->etop[j]=0;
  }

  switch(prop->vertcoord)
  {
    case 0: // user defined
      InitializeVerticalCoordinate(grid,prop,phys,myproc);  
      break;    
    case 2:
      InitializeIsopycnalCoordinate(grid,prop,phys,comm,myproc);
      break;
    case 3:
      InitializeSigmaCoordinate(grid, prop, phys, myproc);
      break;
    case 4:
      //InitializeIsopycnalCoordinate(grid,prop,phys,comm,myproc);   
      
      //intializing with isopycnal code for now
      InitializeVariationalCoordinate(grid,prop,phys,comm,myproc);
      break;
  } 

  // set the old value as the initial condition
  for(i=0;i<grid->Nc;i++)
    for(k=grid->ctop[i];k<grid->Nk[i];k++){
      if(!prop->readOldVelocity || prop->vertcoord!=4){
        grid->dzzold[i][k]=grid->dzz[i][k];
       //grid->dzzold2[i][k]=grid->dzz[i][k];
      }
      if(0){
        Generate_dzzolds_varationalvoid(grid, prop, phys);

      }
    }

  ISendRecvCellData3D(grid->dzz, grid, myproc, comm);
  ISendRecvCellData3D(grid->dzzold, grid, myproc, comm);

}


/*
 * Function: UpdateLayerThickness
 * Calculate dz based on the vertical coordinate chosen
 * ----------------------------------------------------
 * The switch is vertcoord 
 * 0 for user defined, 1 for z level, 2 for isopycnal,3 for sigma
 * need to modify grid.c to make sure ctop=0 and Nk=Nkmax for all cells 
 */
void UpdateLayerThickness(gridT *grid, propT *prop, physT *phys, int index, int myproc, int numprocs, MPI_Comm comm)
{
  int i,k,j,nf,ne,dry,max_ind,neigh,num_bot_layers;
  REAL fac1,fac2,fac3,Ac,sum,sum_old,maxdzz,u_im,alpha_f,error,dzf,dzz_iso,dzz_sig,d1,d2,dsigma;

  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  // ctop and etop is always 0 and Nke and Nk is always Nkmax
  if(!index)
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++){
        //here is the problem it seems
          //grid->dzzold2[i][k]=grid->dzzold[i][k];
        grid->dzzold[i][k]=grid->dzz[i][k];
        vert->zcold[i][k]=vert->zc[i][k];
	      vert->zlold[i][k]=vert->zl[i][k];	
      }
    }

  switch(prop->vertcoord)
  {
    case 0:
      UserDefinedVerticalCoordinate(grid,prop,phys,myproc);
      break;
    case 2:
      for(i=0;i<grid->Nc;i++) {
        for(k=grid->ctop[i];k<grid->Nk[i];k++){
          for(nf=0;nf<grid->nfaces[i];nf++) {               
            ne = grid->face[i*grid->maxfaces+nf];
	    
	    u_im = fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k];	    

	    /*
	    if(u_im>=0)
	      dzf=phys->SfHp[ne][k];
	    else
	      dzf=phys->SfHm[ne][k];

	    // Must set grid->dzf to dzf to ensure that Omega vanishes at each face in the layer-averaged
	    // continuity equation
	    grid->dzf[ne][k]=dzf;
	    */
	    grid->dzz[i][k]-=prop->dt*u_im*grid->dzf[ne][k]*
	      grid->normal[i*grid->maxfaces+nf]*grid->df[ne]/grid->Ac[i];
	  }
	}
      }
	
      break;
    case 3:
      for(i=0;i<grid->Nc;i++)
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  grid->dzz[i][k]=vert->dsigma[k]*(phys->h[i]+grid->dv[i]-phys->zB[i]);

      break;
    case 4:
  
    //HYBRID VERTICAL LAYERS   
    // if num_bot_layers, bot_stretch, and num_mid_layers =0, this gives normal variational layers.
    num_bot_layers = prop->botLayers; //40, 4, 30 for most sand wave runs. 40 for 90 layer sandwaves
    REAL num_mid_layers = prop->midLayers; //15, 14, 17??, 26 for ape1,sh2
    REAL cutoffpoint, sum2;
    
    int bot_stretch = 0; //1
    if(num_mid_layers > 0)
      bot_stretch = 1; //currently have bottom stretching on always

    if(bot_stretch){ //calc dsig before update
      //cutoffpoint = 270; //270 for most sand wave runs, 280 for Nk90 config
      //cutoffpoint = 300 - 30.8611; //not quite 270 because of the stretching factor
      cutoffpoint = prop->cuttoffDepth;

      i=1;  
      sum=0;
      for(k=0;k<grid->Nk[i]-(num_bot_layers+num_mid_layers);k++){
	sum+=grid->dzz[i][k];
      }
      for(k=grid->Nk[i]-(num_bot_layers+num_mid_layers);k<(grid->Nk[i]-num_bot_layers);k++){
	vert->dsigma[k]=grid->dzz[i][k]/(cutoffpoint-sum+phys->h[i]); //find dsig for middle layers if n=1
      }

      //and for bottom layers
      for(k=grid->Nk[i]-(num_bot_layers);k<(grid->Nk[i]);k++){
	vert->dsigma[k]=grid->dzz[i][k]/(grid->dv[i]-cutoffpoint-phys->zBold[i]); //find dsig for bottom layers if n=1
      }

    }

      for(i=0;i<grid->Nc;i++) {
        for(k=grid->ctop[i];k<grid->Nk[i]-(num_bot_layers+num_mid_layers);k++) {
          for(nf=0;nf<grid->nfaces[i];nf++) {               
            ne = grid->face[i*grid->maxfaces+nf];

	    u_im = fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k];

	    
	    // if(u_im>=0)
	    //   dzf=phys->SfHp[ne][k];
	    // else
	    //   dzf=phys->SfHm[ne][k];	    

	    // // Can't set grid->dzf to dzf if also trying to ensure sum(dzz^{n+1})=H^{n+1} because this requires
	    // // Depth-averaged continuity equaiton is satisfied, which is based on the original
	    // // dzf called in SetFluxeHeight
	    // grid->dzf[ne][k]=dzf;
	    
	    grid->dzz[i][k]-=prop->dt*u_im*grid->dzf[ne][k]*
	      grid->normal[i*grid->maxfaces+nf]*grid->df[ne]/grid->Ac[i];
          }

	  // Fix the layer heights if smaller than vertdzmin
	  if(grid->dzz[i][k]<vert->vertdzmin)
	    grid->dzz[i][k]=vert->vertdzmin;
	}
      }

  if(bot_stretch){
    //cutoffpoint = 270; //270 for most sand wave runs, 280 for 90 comfig
    //cutoffpoint = 300 - 30.8611; //not quite 270 because of the stretching factor
    cutoffpoint = prop->cuttoffDepth;

    for(i=0;i<grid->Nc;i++){
      sum=0;
      sum2=0;
      for(k=0;k<grid->Nk[i]-(num_bot_layers+num_mid_layers);k++){
        sum+=grid->dzz[i][k];
      }
      for(k=grid->Nk[i]-(num_bot_layers+num_mid_layers);k<(grid->Nk[i]-num_bot_layers);k++){
        grid->dzz[i][k]=vert->dsigma[k]*(cutoffpoint-sum+phys->h[i]);
        sum2+=grid->dzz[i][k];
      }
      alpha_f = (cutoffpoint-sum+phys->h[i])/sum2;
      for(k=grid->Nk[i]-(num_bot_layers+num_mid_layers);k<(grid->Nk[i]-num_bot_layers);k++) {
        grid->dzz[i][k]*=alpha_f;   
      }

      sum2=0;
      alpha_f=0;
      for(k=grid->Nk[i]-(num_bot_layers);k<grid->Nk[i];k++){
        //grid->dzz[i][k]=(1.0/num_bot_layers)*(grid->dv[i]-cutoffpoint-phys->zB[i]);
        grid->dzz[i][k]=vert->dsigma[k]*(grid->dv[i]-cutoffpoint-phys->zB[i]);
        //grid->dzz[i][k]=grid->dzz[i][k]; //stay same? 
        sum2+=grid->dzz[i][k];
      }
      alpha_f = (grid->dv[i]-cutoffpoint-phys->zB[i])/sum2;
      for(k=grid->Nk[i]-(num_bot_layers);k<grid->Nk[i];k++) {
        grid->dzz[i][k]*=alpha_f;   
      }

    }


  } else if(num_bot_layers){ //normal hybrid layers
    //update bottom layers
    cutoffpoint = 300; //270
    for(i=0;i<grid->Nc;i++){
    sum=0;
    for(k=0;k<grid->Nk[i]-(num_bot_layers+num_mid_layers);k++){
      sum+=grid->dzz[i][k];
    }
    for(k=grid->Nk[i]-(num_bot_layers+num_mid_layers);k<(grid->Nk[i]-num_bot_layers);k++){
      grid->dzz[i][k]=(1.0/num_mid_layers)*(cutoffpoint-sum+phys->h[i]);
    }
	  for(k=grid->Nk[i]-(num_bot_layers);k<grid->Nk[i];k++){
	    //grid->dzz[i][k]=(1.0/num_bot_layers)*(grid->dv[i]-cutoffpoint-phys->zB[i]);
	    grid->dzz[i][k]=grid->dzz[i][k]; //stay same? 
	  }
    }

    } else{
      //Normal variational layers 
      //Update to make sure sum of layer heights matches total depth
      for(i=0;i<grid->Nc;i++){
        sum=0;
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          sum+=grid->dzz[i][k];
          }

          alpha_f = (phys->h[i]+grid->dv[i]-phys->zB[i])/sum;
          for(k=grid->ctop[i];k<grid->Nk[i];k++) {
            grid->dzz[i][k]*=alpha_f;
          }
        }
    }
    //through here 
    

    //check grid 
    // for(i=0;i<grid->Nc;i++){
    //   sum=0;
    //   for(k=grid->ctop[i];k<grid->Nk[i];k++) {
    //     sum+=grid->dzz[i][k];
    //     }
    //     //printf("sum of layer heights at cell [%d] is %f at step %d \n", i, sum, prop->n);

    //     alpha_f = ((phys->h[i]+grid->dv[i]-phys->zB[i]) - sum)*1000;
    //     printf("difference in heights (expected total - total)*1000 = %f for i = %d \n", alpha_f, i);

    //   }
      

      
      break;
    case 5:
      UpdateDZ(grid,phys,prop, 0); 
  }
}


/*
 * Function: VariationalVertCoordinate
 * use the variational moving mesh approach (Tang & Tang, 2003) and (Koltakov & Fringer, 2013)
 * ----------------------------------------------------
 * here not solve full elliptic equation but averaging in horizontal equation
 */
void VariationalVertCoordinateAverageMethod(gridT *grid, propT *prop, physT *phys, int myproc)
{
   int i,k,nf,Neigh;
   REAL thetaT=1,thetaS=0.5,sum;
   MonitorFunctionForAverageMethod(grid, prop, phys, myproc);
   // use the monitor function first
   for(i=0;i<grid->Nc;i++)
   {
     sum=0;
     Neigh=0;
     for(nf=0;nf<grid->nfaces[i];nf++)
       if(grid->neigh[i*grid->maxfaces+nf]!=-1)
        Neigh++;
     for(k=0;k<grid->Nk[i];k++){
        vert->dzztmp[i][k]=0;
        // use the monitor function
        grid->dzz[i][k]=(phys->h[i]+grid->dv[i])/vert->Msum[i]/vert->Mc[i][k];
        // time averaging with the old step to slow down the time varying
        grid->dzz[i][k]=thetaT*grid->dzz[i][k]+(1-thetaT)*grid->dzzold[i][k];
        grid->dzz[i][k]=(phys->h[i]+grid->dv[i])/((1-thetaT)*phys->h_old[i]+thetaT*phys->h[i]+grid->dv[i])*grid->dzz[i][k];
        // horizontal averaging
        for(nf=0;nf<grid->nfaces[i];nf++){
          if(grid->neigh[i*grid->maxfaces+nf]!=-1)
            vert->dzztmp[i][k]+=grid->dzz[grid->neigh[i*grid->maxfaces+nf]][k]/Neigh;
        }
        grid->dzz[i][k]=grid->dzz[i][k]*(1-thetaS)+vert->dzztmp[i][k]*thetaS;
        sum+=grid->dzz[i][k];
     }
     for(k=0;k<grid->Nk[i];k++){
       grid->dzz[i][k]=(phys->h[i]+grid->dv[i])/sum*grid->dzz[i][k];

     }
   }
}

/*
 * Function: VariationalVertCoordinate
 * use the variational moving mesh approach (Tang & Tang, 2003) and (Koltakov & Fringer, 2013)
 * ----------------------------------------------------
 * here not solve full elliptic equation but averaging in horizontal equation
 */
void IsoVariational(gridT *grid, propT *prop, physT *phys, int myproc, int numprocs, MPI_Comm comm) {
  int i, iptr, k, nf, ne, dry, max_ind, iter, neigh, maxiters, normal;
  REAL fac1, fac2, fac3, maxdzz, sum, alpha, S, Smax, maxSlope, H_old, H_new, df, Ac, dzf;

  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    dry=0;
    sum=0;
    max_ind=grid->ctop[i];
    maxdzz=grid->dzz[i][grid->ctop[i]];    
    Ac = grid->Ac[i];
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {

      for(nf=0;nf<grid->nfaces[i];nf++) {               
	ne = grid->face[i*grid->maxfaces+nf];
	df = grid->df[ne];
        dzf = grid->dzf[ne][k];	
	normal = grid->normal[i*grid->maxfaces+nf];
	if(k<grid->Nke[ne]) {
	  grid->dzz[i][k]-=
	    prop->dt*(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k])*dzf*df*normal/Ac;
	}
      }
      if(grid->dzz[i][k]<vert->vertdzmin){
	dry++;
	grid->dzz[i][k]=vert->vertdzmin;
      }
      sum+=grid->dzz[i][k];
      if(k>grid->ctop[i]){
	if(grid->dzz[i][k]>maxdzz){
	  max_ind=k;
	  maxdzz=grid->dzz[i][k];
	}
      }
    }
    if(dry) {
      //      if(VERBOSE>1 && myproc==0)
      //	printf("Step %d, correcting %d dry cells at x=%.2f\n",prop->n,dry,grid->xv[i]);

      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	grid->dzz[i][k]*=(grid->dv[i]+phys->h[i])/sum;
      for(k=grid->ctop[i];k<grid->Nk[i];k++)
	if(k==max_ind)
	  grid->dzz[i][k]+=(grid->dv[i]+phys->h[i]-sum);
    }
  }

  ISendRecvCellData3D(grid->dzz,grid,myproc,comm);
  /*
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    vert->zl[i][grid->Nk[i]]=-grid->dv[i];

    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      vert->zl[i][k]=vert->zl[i][k+1]+grid->dzz[i][k];
    }
  }
  ISendRecvCellData3D(grid->dzz,grid,myproc,comm);    
  */
  /*  
  for(iter=0;iter<10;iter++) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      
      
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	phys->a[k]=grid->dzz[i][k]*grid->Ac[i];
	phys->b[k]=grid->Ac[i];
      }
      for(nf=0;nf<grid->nfaces[i];nf++) {               
	ne = grid->neigh[i*grid->maxfaces+nf];
	
	if(ne!=-1) {
	  for(k=grid->ctop[ne];k<grid->Nk[ne];k++) {
	    phys->a[k]+=grid->Ac[ne]*grid->dzz[ne][k];
	    phys->b[k]+=grid->Ac[ne];
	  }
	}
      }
      H_old = H_new = 0;
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {    
	H_old+=grid->dzz[i][k];
	grid->dzz[i][k]=phys->a[k]/phys->b[k];
	H_new+=grid->dzz[i][k];
      }
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	grid->dzz[i][k]*=H_old/H_new;
      }
    }
    ISendRecvCellData3D(grid->dzz,grid,myproc,comm);
  }
  ComputeZc(grid,prop,phys,1,myproc);
  */

  /* Iterate and smooth the zeta levels
  ComputeZc(grid,prop,phys,1,myproc);

  maxiters=10000;
  maxSlope=INFTY;
  Smax = tan(10*PI/180); // Maximum isopycnal slope

  iter=0;
  while(maxSlope>Smax && iter<maxiters) {
    maxSlope=0;
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(nf=0;nf<grid->nfaces[i];nf++) {
	ne = grid->face[i*grid->maxfaces+nf];
	neigh = grid->neigh[i*grid->maxfaces+nf];
	
	if(neigh!=-1) {
	  for(k=grid->ctop[neigh]+1;k<grid->Nk[neigh];k++) {
	    S = (vert->zl[neigh][k]-vert->zl[i][k])/grid->dg[ne];

	    if(fabs(S)>=maxSlope) maxSlope = fabs(S);	    

	    //	    if(phys->rho[i][k-1]==phys->rho[i][k] || S>Smax)
	    if(S>Smax)
	      alpha=100;
	    else
	      alpha=0.0;
	    // alpha ~ 0.5 Ac dg/df
	    //alpha = df/(dg*Ac);
	    vert->zl[i][k]+=alpha*S*grid->df[ne]/grid->Ac[i];
	    //	    vert->zl[i][k]+=alpha*(vert->zl[neigh][k]-vert->zl[i][k]);	      
	  }
	}
      }
    }
    if(iter==maxiters && myproc==0)
      printf("Warning - grid smoothing did not converge after %d iterations.\n",maxiters);
    //    printf("Time step %d, Iteration %d, maxSlope = %.3e, Smax=%.3e\n",prop->n,iter,maxSlope,Smax);

    ISendRecvCellData3D(vert->zl,grid,myproc,comm);    
    iter++;
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      grid->dzz[i][k]=vert->zl[i][k]-vert->zl[i][k+1];
    }
  }
  ComputeZc(grid,prop,phys,1,myproc);  
  ISendRecvCellData3D(grid->dzz,grid,myproc,comm);
  */
}

/*
 * Function: FindBottomLayer
 * Calculate Nkeb and zfb for each edge based on dzf
 * ----------------------------------------------------
 * these two values will be further used in SetDragCoefficients and calculate drag terms
 */
void FindBottomLayer(gridT *grid, propT *prop, physT *phys, int myproc) {
  int j,k,Nkeb,nc1,nc2;
  REAL zfb;
  
  for(j=0;j<grid->Ne;j++) {

    nc1=grid->grad[2*j];
    nc2=grid->grad[2*j+1];
    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;
    
    zfb=0;
    for(k=grid->Nke[j]-1;k>=grid->etop[j];k--) {
      zfb+=0.5*(grid->dzz[nc1][k]+grid->dzz[nc2][k]);

      // find the layer when zfb>buffer height
      if(zfb>=prop->BUFFERHEIGHT) {
	vert->Nkeb[j]=k;
	vert->zfb[j]=zfb;
	break;
      }

      // if not find at the top layer, it means the total dzf is smaller than bufferheight=dry
      // assume the drag layer at top layer
      if(k==grid->etop[j]) {
	vert->Nkeb[j]=k;
	vert->zfb[j]=zfb;
      } 
    }
    //    if(grid->mark[j]==0) printf("zfb = %.3e, Nkeb = %d, buffer=%.3e\n",
    //				vert->zfb[j],vert->Nkeb[j],prop->BUFFERHEIGHT);
  }
}

/*
 * Function: ComputeUl
 * Calculate u v at the top and bottom face of each layer
 * ----------------------------------------------------
 * compute ul vl
 */
void ComputeUl(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,j,k;
  for(i=0;i<grid->Nc;i++)
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
    {
       vert->ul[i][k]=InterpToLayerTopFace(i, k, phys->uc, grid);
       vert->vl[i][k]=InterpToLayerTopFace(i, k, phys->vc, grid);
    }
}

/*
 * Function: ComputeUf
 * Calculate the u v w omegaf at the edge center of each layer and 
 * u v at the top and bottom face of each layer
 * ----------------------------------------------------
 * compute uf vf wf and ul vl
 */
void ComputeUf(gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,j,k,nc1,nc2;
  REAL sum,r1,r2,def1,def2;
  // compute omegac
  for(i=0;i<grid->Nc;i++)
    for(k=grid->ctop[i];k<grid->Nk[i];k++){
      phys->wc[i][k]=(phys->w[i][k]+phys->w[i][k+1])/2;
      vert->omegac[i][k]=(vert->omega[i][k]+vert->omega[i][k+1])/2;
    }

  // uf vf and wf
  for(j=0;j<grid->Ne;j++){
    nc1=grid->grad[2*j];
    nc2=grid->grad[2*j+1];
    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;
    Return_def(&def1,&def2,nc1,nc2,j,grid);	   	            
    //    def1 = grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
    //    def2 = grid->def[nc2*grid->maxfaces+grid->gradf[2*j+1]];
    for(k=grid->etop[j];k<grid->Nke[j];k++)
    {
      vert->uf[j][k]=0.5*(phys->uc[nc1][k]+phys->uc[nc2][k]);
      vert->vf[j][k]=0.5*(phys->vc[nc1][k]+phys->vc[nc2][k]);
      vert->wf[j][k]=0.5*(phys->w[nc1][k]+phys->w[nc2][k]);
      vert->qcf[j][k]=0.5*(phys->qc[nc1][k]+phys->qc[nc2][k]);
      vert->omegaf[j][k]=0.5*(vert->omegac[nc1][k]+vert->omegac[nc2][k]);
    }
    // wf is for dwdx and dwdy which needs to cover the top of each cell layer 
    if(grid->ctop[nc1]!=grid->ctop[nc2])
    {
      if(grid->ctop[nc1]>grid->ctop[nc2])
        nc1=nc2;
      for(k=grid->ctop[nc1];k<grid->etop[j];k++)
        vert->wf[j][k]=0.5*phys->w[nc1][k];
    }
  }
}

/*
 * Function: InterpToLayerTopFace
 * Usage: uface = InterpToLayerTopBotFace(j,k,phys->uc,u,grid);
 * -------------------------------------------------
 * Linear interpolation of a Voronoi-centered value to the top and bottom face of each layer
 * 
 * phi_top = (dz1*u2 + dz2*u1)/(dz1+dz2);
 *
 */
REAL InterpToLayerTopFace(int i, int k, REAL **phi, gridT *grid) {
  int nc1, nc2;
  REAL def1, def2, Dj;
  if(k>grid->ctop[i])
  {
    Dj = grid->dzz[i][k]+grid->dzz[i][k-1];
    def1=grid->dzz[i][k-1];
    def2=grid->dzz[i][k];
  }else{
    Dj = grid->dzz[i][k]/2+grid->dzz[i][k+1]/2;
    def1=0;
    def2=grid->dzz[i][k]/2;    
  }

  if(def2==0 && k<grid->Nk[i]-1 && k!=grid->ctop[i]) {
    // means the interior layer is dry
    return phi[i][k+1];
  } else if(def2==0 && k==grid->Nk[i]-1){
    // means the bottom layer is dry
    return phi[i][k];
  } else if(k==grid->ctop[i]) {
    return phi[i][k]+def2*(phi[i][k]-phi[i][k+1])/Dj;
  } else {
    return (phi[i][k-1]*def2+phi[i][k]*def1)/Dj;
  }
}

/*
 * Function: InterpToLayerTopFace
 * Usage: uface = InterpToLayerTopBotFace(j,k,phys->uc,u,grid);
 * -------------------------------------------------
 * Linear interpolation of a Voronoi-centered value to the top and bottom face of each layer
 * 
 * phi_top = (dz1*u2 + dz2*u1)/(dz1+dz2);
 *
 */
REAL InterpToLayerTopFaceOldStep(int i, int k, REAL **phi, gridT *grid) {
  int nc1, nc2;
  REAL def1, def2, Dj;
  if(k>grid->ctop[i])
  {
    Dj = grid->dzzold[i][k]+grid->dzzold[i][k-1];
    def1=grid->dzzold[i][k-1];
    def2=grid->dzzold[i][k];
  }else{
    Dj = grid->dzzold[i][k]/2+grid->dzzold[i][k+1]/2;
    def1=0;
    def2=grid->dzzold[i][k]/2;    
  }

  if(def2==0 && k<grid->Nk[i]-1 && k!=grid->ctop[i]) {
    // means the interior layer is dry
    return phi[i][k+1];
  } else if(def2==0 && k==grid->Nk[i]-1){
    // means the bottom layer is dry
    return phi[i][k];
  } else if(k==grid->ctop[i]) {
    return phi[i][k]+def2*(phi[i][k]-phi[i][k+1])/Dj;
  } else {
    return (phi[i][k-1]*def2+phi[i][k]*def1)/Dj;
  }
}

/*
 * Function: LayerAveragedContinuity
 * Compute omega for for the top and bottom face of each layer
 * from the layer-averaged continuity equation
 *
 * Note that omega hould be zero for isopycnals but it is not because
 * of the correction of small layer heights.
 * ----------------------------------------------------
 * 
 */
void LayerAveragedContinuity(REAL **omega, gridT *grid, propT *prop, physT *phys, MPI_Comm comm, int myproc)
{
  int i,iptr,k,nf,ne;
  REAL fac1,fac2,fac3,flux,Ac;

  if(prop->n==1){
    fac1=1;
    fac2=0;
    fac3=0;
  }else if(prop->n==2){
    fac1=prop->theta;
    fac2=1.0-prop->theta;
    fac3=0;
  }else{
    fac1=prop->imfac1;
    fac2=prop->imfac2;
    fac3=prop->imfac3;
  }

    fac1=prop->imfac1;
    fac2=prop->imfac2;
    fac3=prop->imfac3;


  //printf("fac1 = %f, fac2 = %f, and fac3 = %f at n=%d \n", fac1, fac2, fac3, prop->n);


  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) 
  {
    i = grid->cellp[iptr];
    
    for(k=0;k<grid->Nk[i];k++){
	if(vert->omega[i][k]!=vert->omega[i][k]){
	  printf("omega nan pre continuity for i=%d, k=%d \n", i, k);
	}
    }
  

    // compute omega from bottom layer using the omega_bot=0
    // worry about the numerical accuracy to ensure omega_top==0
    Ac=grid->Ac[i];
    omega[i][grid->Nk[i]]=0;
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
    {
      flux=0;
      int flag=-1;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
	if(grid->mark[ne]!=0) flag=grid->mark[ne];
	//	if(grid->mark[ne]==1 && fabs(phys->u[ne][k])!=0) { printf("ERROR U NONZERO!\n"); exit(1); }
        if(k<grid->Nke[ne])
          flux+=(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
      }
      REAL div = (grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt+flux;
      
      // if(prop->n<4){
      //   printf("omega_old = %f, and omega_old2 = %f for i=%d, k=%d \n", vert->omega_old[i][k], vert->omega_old2[i][k], i, k);
      // }

      if(!prop->subgrid)
        omega[i][k]=(fac1*omega[i][k+1]+fac2*vert->omega_old[i][k+1]+fac3*vert->omega_old2[i][k+1])-
        (fac2*vert->omega_old[i][k]+fac3*vert->omega_old2[i][k])-flux-
        (grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt;
      else{
        omega[i][k]=(fac1*omega[i][k+1]*subgrid->Acveff[i][k+1]+
          fac2*vert->omega_old[i][k+1]*subgrid->Acveffold[i][k+1]+
          fac3*vert->omega_old2[i][k+1]*subgrid->Acveffold2[i][k+1])-
        (fac2*vert->omega_old[i][k]*subgrid->Acveffold[i][k]+
          fac3*vert->omega_old2[i][k]*subgrid->Acveffold2[i][k])-flux*Ac-
        (subgrid->Acceff[i][k]*grid->dzz[i][k]-subgrid->Acceffold[i][k]*grid->dzzold[i][k])/prop->dt;        
        omega[i][k]/=subgrid->Acveff[i][k];
      }
      omega[i][k]/=fac1;


      //omega[i][k]=0; // OBF 010424 -- THIS IS FOR DEBUGGING!!! Should remove; omega=0 only for isopycnals.
      /*
      if(fabs(div)>1e-10 && flag>-1 && flag!=5)
	printf("Nonzero div = %.3e, Omega=%.3e, dzz=%.3e flag = %d, at i,k=%d,%d\n",div,omega[i][k],grid->dzz[i][k],flag,i,k);
      */
    }

    
    //compute omega_im for updatescalar function
    for(k=grid->ctop[i];k<grid->Nk[i];k++){
      if(!prop->subgrid) 
        vert->omega_im[i][k]=fac1*omega[i][k]+fac2*vert->omega_old[i][k]+fac3*vert->omega_old2[i][k];
      else {
        vert->omega_im[i][k]=(fac2*vert->omega_old[i][k]*subgrid->Acveffold[i][k]+
                fac3*vert->omega_old2[i][k]*subgrid->Acveffold2[i][k]+
                fac1*omega[i][k]*subgrid->Acveff[i][k])/subgrid->Acveff[i][k];
      }
    }
    vert->omega_im[i][grid->Nk[i]]=0;
  }

  //printf("made it through omega_im set on proc %d \n", myproc);

  for(i=0;i<grid->Nc;i++){
    for(k=0;k<grid->Nk[i];k++){
      if(omega[i][k]!=omega[i][k]){
	printf("omega nan post continuity for i=%d, k=%d \n", i, k);
      }
    }
  }

  //printf("omega nan check on proc %d \n", myproc);

  ISendRecvCellData3D(vert->omega_im,grid,myproc,comm);	
  //printf("sent omega_im set on proc %d \n", myproc);

}

/*
 * Function: ComputeVerticalVelocity
 * Compute w for for the top and bottom face of each layer
 * from the layer-averaged continuity equation
 * ----------------------------------------------------
 * 
 */
void ComputeVerticalVelocity(REAL **omega, REAL **zc, gridT *grid, propT *prop, physT *phys, int myproc)
{
  int i,k,nf,ne;
  REAL fac1,fac2,fac3,flux,Ac;
  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  for(i=0;i<grid->Nc;i++)
  {
    Ac=grid->Ac[i];
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--){
      flux=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        if(k<grid->Nke[ne])
          flux+=(fac1*phys->u[ne][k]+fac2*phys->u_old[ne][k]+fac3*phys->u_old2[ne][k])*grid->df[ne]*grid->normal[i*grid->maxfaces+nf]/Ac*grid->dzf[ne][k];
      }
      omega[i][k]=vert->omega[i][k+1]-flux-(grid->dzz[i][k]-grid->dzzold[i][k])/prop->dt;
    }
  }
}

/*
 * Function: ComputeZc
 * Compute the vertical location of each cell center
 * also compute the vertical location of each layer center at edge face
 * ----------------------------------------------------
 */
void ComputeZc(gridT *grid, propT *prop, physT *phys, int index, int myproc)
{
  int i,j,k,nc1,nc2;
  REAL z,def1,def2,Dj;
  for(i=0;i<grid->Nc;i++)
  {
    z=-grid->dv[i]+phys->zB[i];
    vert->zc[i][grid->Nk[i]-1]=z+grid->dzz[i][grid->Nk[i]-1]/2;
    for(k=grid->Nk[i]-2;k>=grid->ctop[i];k--){
      vert->zc[i][k]=vert->zc[i][k+1]+grid->dzz[i][k+1]/2+grid->dzz[i][k]/2;
    }
    if(prop->n==prop->nstart)
      for(k=0;k<grid->Nk[i];k++) 
        vert->zcold[i][k]=vert->zc[i][k];
 
    if(!index){
      z=-grid->dv[i]+phys->zBold[i];
      vert->zcold[i][grid->Nk[i]-1]=z+grid->dzzold[i][grid->Nk[i]-1]/2;
      for(k=grid->Nk[i]-2;k>=grid->ctop[i];k--){
        vert->zcold[i][k]=vert->zcold[i][k+1]+grid->dzzold[i][k+1]/2+grid->dzzold[i][k]/2;
      }      
    }

    vert->zl[i][grid->Nk[i]]=-grid->dv[i]+phys->zB[i];
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      vert->zl[i][k]=vert->zl[i][k+1]+grid->dzz[i][k];
    }
    if(prop->n==prop->nstart)
      for(k=0;k<grid->Nk[i];k++) 
        vert->zlold[i][k]=vert->zl[i][k];

    if(!index){
      vert->zlold[i][grid->Nk[i]]=-grid->dv[i]+phys->zBold[i];
      for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
	vert->zlold[i][k]=vert->zlold[i][k+1]+grid->dzzold[i][k];
      }      
    }

    if(fabs(vert->zl[i][grid->ctop[i]]-phys->h[i])>1e-9*grid->Nkmax){
      //printf("something wrong with the calculation of layerthickness at cell %d %e\n",
        //i,vert->zl[i][grid->ctop[i]]-phys->h[i]);
      //exit(0);
    }
    // OBF 012623    grid->dzz[i][0]+=phys->h[i]-vert->zl[i][grid->ctop[i]];
    // OBF 012623   vert->zl[i][grid->ctop[i]]=phys->h[i];
  }

  for(j=0;j<grid->Ne;j++)
  {
    nc1=grid->grad[2*j];
    nc2=grid->grad[2*j+1];
    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;
    Dj = grid->dg[j];
    Return_def(&def1,&def2,nc1,nc2,j,grid);	   	                
    //    def1 = grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
    //    def2 = grid->def[nc2*grid->maxfaces+grid->gradf[2*j+1]];
    for(k=0;k<=grid->Nke[j];k++){
      vert->zf[j][k]=vert->zl[nc1][k]*def2/Dj+vert->zl[nc2][k]*def1/Dj;
    }
  }
}

/*
 * Function: ComputeOmega
 * Compute omega from the definition omega=w-udzdx-vdzdy-wg and U3=w-udzdx-vdzdy
 * ----------------------------------------------------
 * index=1, w->omega index=0 omega->w index=-1 w->U3
 * index=1 is never used in phys.c
 * omega is always calculated from continuity
 */
void ComputeOmega(gridT *grid, propT *prop, physT *phys, int index, int myproc)
{
  int i,k;
  REAL alpha,alpha_s=0, dzdt; // alpha=0 do not include slope terms in U3
  if(prop->vertcoord==5)
    alpha=0;

  //alpha = 1;
  alpha = 0; 

  if(index==1) {
    // Compute omega = w - u dz/dx - v dz/dy - dz/dt
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
        vert->omega[i][k]=phys->w[i][k]
	  -alpha*vert->ul[i][k]*vert->dzdx[i][k]
	  -alpha*vert->vl[i][k]*vert->dzdy[i][k]
	  -(vert->zl[i][k]-vert->zlold[i][k])/prop->dt;
      }
    }
  } else if(index==0) {
    // Compute w = omega + u dz/dx + v dz/dy + dz/dt
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
        phys->w[i][k]=vert->omega[i][k]
	  +alpha*vert->ul[i][k]*vert->dzdx[i][k]
	  +alpha*vert->vl[i][k]*vert->dzdy[i][k]
	  +(vert->zl[i][k]-vert->zlold[i][k])/prop->dt;
      }
    }
  } else {
    // Compute U3 = w - u dz/dx - v dz/dy
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->Nk[i]-1;k>grid->ctop[i];k--) {
        vert->U3[i][k]=phys->w[i][k]-
	  (alpha*0.5*(grid->dzz[i][k]*phys->uc[i][k-1]+grid->dzz[i][k-1]*phys->uc[i][k])*vert->dzdx[i][k]+
	   alpha*0.5*(grid->dzz[i][k]*phys->vc[i][k-1]+grid->dzz[i][k-1]*phys->vc[i][k])*vert->dzdy[i][k])/
	  (grid->dzz[i][k-1]+grid->dzz[i][k]);
      }
      k=grid->ctop[i];
      vert->U3[i][k]=phys->w[i][k]-
	(alpha*phys->uc[i][k]*vert->dzdx[i][k]+alpha*phys->vc[i][k]*vert->dzdy[i][k]);
	
      k=grid->Nk[i];
      // BDF2 scheme to obtain dzdt from z^{n+1}, z^{n}, and z^{n-1}
      dzdt = (1.5*phys->zB[i] - 2.0*phys->zBold[i] + 0.5*phys->zBold2[i])/prop->dt;
      /* This needs to be done in w predictor
      phys->w[i][k]=dzdt
	+alpha*phys->uc[i][k-1]*vert->dzdx[i][k-1]
      	+alpha*phys->vc[i][k-1]*vert->dzdy[i][k-1];
      */
      vert->U3[i][k]=dzdt;
    }
  }
}

/*
 * Function: ComputeOmega
 * Compute omega from the definition omega=w-udzdx-vdzdy-wg and U3=w-udzdx-vdzdy
 * ----------------------------------------------------
 * index=1, w->omega index=0 omega->w index=-1 w->U3
 * index=1 is never used in phys.c
 * omega is always calculated from continuity
 */
void ComputeOmegaOlds(gridT *grid, propT *prop, physT *phys, MPI_Comm comm,int myproc)
{
  int i,j,k;
  REAL alpha,alpha_s=0, dzdt; // alpha=0 do not include slope terms in U3
  REAL **dzzold2;

  //when we enter, omega = n, omega_old = omega_n-1 and omega_old2 = omega_n-2
  //after it will shift so omega_old = omega_n, omega_old2 = omega_old
  
  //and w_old = w, w_old2 = w_old

  if(prop->vertcoord==5)
    alpha=0;

  //alpha = 1;
  alpha = 0;

  /*start with n-1*/


  //first get ul, vl
  for(i=0;i<grid->Nc;i++)
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
    {
      //for n-1
       vert->ul[i][k]=InterpToLayerTopFaceOldStep(i, k, phys->uold, grid);
       vert->vl[i][k]=InterpToLayerTopFaceOldStep(i, k, phys->vold, grid);

       if(vert->ul[i][k] != vert->ul[i][k]){
        printf("ul is nan for i=%d, k=%d \n", i, k);
       }

       if(vert->vl[i][k] != vert->vl[i][k]){
        printf("vl is nan for i=%d, k=%d \n", i, k);
       }
    }

  //for some reason, can't add dzzold2 variable globally so we will read it in once here. 
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];
  FILE *fid;


  //allocate space for dzzold2
  dzzold2 = (REAL **)malloc(grid->Nc*sizeof(REAL *));
  for(i=0;i<grid->Nc;i++) {
    dzzold2[i] = (REAL *)malloc(grid->Nk[i]*sizeof(REAL));
  } 

  //add in to read in dzz_old2
  if(prop->vertcoord==4){
    MPI_GetFile(filename,DATAFILE,"dzz_t2_init_file","ComputeOmegaOlds",myproc);  
    sprintf(str,"%s.%d",filename,myproc);  
    fid = MPI_FOpen(str,"r","ComputeOmegaOlds",myproc);

    for(i=0;i<grid->Nc;i++) {
      fread(dzzold2[i],sizeof(REAL),grid->Nkmax,fid);
    }
    fclose(fid);
  } else{
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        dzzold2[i][k] = grid->dzzold[i][k];
        if(grid->ctop[i] != 0)
          printf("grid->ctop[i] = %d for i = %d, k=%d \n", grid->ctop[i], i, k);
      }
    }
  }

  ISendRecvCellData3D(dzzold2,grid,myproc,comm);
  
  //next need zl and zlold (copied the relevatn part of ComputeZc)
  for(i=0;i<grid->Nc;i++){
  
    vert->zl[i][grid->Nk[i]]=-grid->dv[i]+phys->zBold[i]; 
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      vert->zl[i][k]=vert->zl[i][k+1]+grid->dzzold[i][k]; //dzzold is n-1 

      if(vert->zl[i][k] != vert->zl[i][k]){
        printf("zl is nan for i=%d, k=%d \n", i, k);
       }
    }



    vert->zlold[i][grid->Nk[i]]=-grid->dv[i]+phys->zBold2[i]; 
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      vert->zlold[i][k]=vert->zlold[i][k+1]+dzzold2[i][k]; //should be dzzold2 but we're getting mem issues
      if(vert->zlold[i][k] != vert->zlold[i][k]){
        printf("zlold is nan for i=%d, k=%d \n", i, k);
       }
    }      
  }


  
//now can do dzdx, dzdy
ComputeCellAveragedHorizontalGradientCell(vert->dzdx, 0, vert->zl, 1, grid, prop, phys, myproc);
ComputeCellAveragedHorizontalGradientCell(vert->dzdy, 1, vert->zl, 1, grid, prop, phys, myproc);



//finally, compute omegaold
// Compute omega = w - u dz/dx - v dz/dy - dz/dt
for(i=0;i<grid->Nc;i++) {
  vert->omega_old[i][grid->Nk[i]]=0;
  for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
    vert->omega_old[i][k]=phys->w_old[i][k]
  -alpha*vert->ul[i][k]*vert->dzdx[i][k]
  -alpha*vert->vl[i][k]*vert->dzdy[i][k]
  -(vert->zl[i][k]-vert->zlold[i][k])/prop->dt; //this will become omega_old2
  }
}

// //and U3_old     
// // Compute U3 = w - u dz/dx - v dz/dy
for(i=0;i<grid->Nc;i++) {
   for(k=grid->Nk[i]-1;k>grid->ctop[i];k--) {
      vert->U3_old[i][k]=phys->w_old[i][k]-
	(alpha*0.5*(grid->dzzold[i][k]*phys->uold[i][k-1]+grid->dzzold[i][k-1]*phys->uold[i][k])*vert->dzdx[i][k]+
	 alpha*0.5*(grid->dzzold[i][k]*phys->vold[i][k-1]+grid->dzzold[i][k-1]*phys->vold[i][k])*vert->dzdy[i][k])/
	(grid->dzzold[i][k-1]+grid->dzzold[i][k]);
   }
   k=grid->ctop[i];
    vert->U3_old[i][k]=phys->w[i][k]-
      (alpha*phys->uold[i][k]*vert->dzdx[i][k]+alpha*phys->vold[i][k]*vert->dzdy[i][k]);

    k=grid->Nk[i];
    vert->U3_old[i][k]=0; //force to 0 here for now
 }


/* Now do current time step n*/

  //first get ul, vl
  for(i=0;i<grid->Nc;i++)
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--)
    {
      //for n-1
        vert->ul[i][k]=InterpToLayerTopFace(i, k, phys->uc, grid);
        vert->vl[i][k]=InterpToLayerTopFace(i, k, phys->vc, grid);
    }

  //next need zl and zlold (copied the relevatn part of ComputeZc)
  for(i=0;i<grid->Nc;i++){
    vert->zl[i][grid->Nk[i]]=-grid->dv[i]+phys->zB[i]; 
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      vert->zl[i][k]=vert->zl[i][k+1]+grid->dzz[i][k];
    }


  vert->zlold[i][grid->Nk[i]]=-grid->dv[i]+phys->zBold[i];
  for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
    vert->zlold[i][k]=vert->zlold[i][k+1]+grid->dzzold[i][k];
  }      
  }
  
  //now can do dzdx, dzdy
  ComputeCellAveragedHorizontalGradientCell(vert->dzdx, 0, vert->zl, 1, grid, prop, phys, myproc);
  ComputeCellAveragedHorizontalGradientCell(vert->dzdy, 1, vert->zl, 1, grid, prop, phys, myproc);

  //finally, compute omegaold
  // Compute omega = w - u dz/dx - v dz/dy - dz/dt
  for(i=0;i<grid->Nc;i++) {
    vert->omega[i][grid->Nk[i]]=0;
    for(k=grid->Nk[i]-1;k>=grid->ctop[i];k--) {
      vert->omega[i][k]=phys->w[i][k]
    -alpha*vert->ul[i][k]*vert->dzdx[i][k]
    -alpha*vert->vl[i][k]*vert->dzdy[i][k]
    -(vert->zl[i][k]-vert->zlold[i][k])/prop->dt; //changed phys->w_old to phys_w
    }
  }
  
  //now put zc's back 
  ComputeZc(grid,prop,phys,0,myproc); //index was 1

  //to compute U3
  ComputeOmega(grid, prop, phys,-1, myproc); //this is U3_n, soon to be U3_old. Now U3_n-1.. 

  for(i=0;i<grid->Nc;i++) {
    free(dzzold2[i]);
    }
  free(dzzold2);


}

/*
 * Function: VerticalCoordinateHorizontalSource
 * Compute all the preparation due to the new general vertical
 * coordinate for the horizontal source of momentum conservation
 * ----------------------------------------------------
 * 
 */
void VertCoordinateHorizontalSource(gridT *grid, physT *phys, propT *prop,
    int myproc, int numprocs, MPI_Comm comm)
{
   int i,k;
   // compute the relative vorticity
   // 1. compute v and u at each edge center
   ComputeUf(grid, prop, phys,myproc);

   // comment out due to use conservative method
   // 2. compute dvdx and dudy
   if(prop->nonlinear)
   {
     ComputeCellAveragedHorizontalGradient(vert->dvdx, 0, vert->vf, grid, prop, phys, myproc);
     ComputeCellAveragedHorizontalGradient(vert->dudy, 1, vert->uf, grid, prop, phys, myproc);
     ComputeCellAveragedHorizontalGradient(vert->dvdy, 1, vert->vf, grid, prop, phys, myproc);
     ComputeCellAveragedHorizontalGradient(vert->dudx, 0, vert->uf, grid, prop, phys, myproc);
     if(prop->nonhydrostatic)
     {
       ComputeCellAveragedHorizontalGradient(vert->dwdy, 1, vert->wf, grid, prop, phys, myproc);
       ComputeCellAveragedHorizontalGradient(vert->dwdx, 0, vert->wf, grid, prop, phys, myproc);
     }
   }
   // 3. compute f_r
   ComputeRelativeVorticity(grid,phys,prop,myproc);
  
   // 4. prepare for nonhydrostatic calculation
   if(prop->nonhydrostatic){
     ComputeCellAveragedHorizontalGradient(vert->dqdy, 1, vert->qcf, grid, prop, phys, myproc);
     ComputeCellAveragedHorizontalGradient(vert->dqdx, 0, vert->qcf, grid, prop, phys, myproc);     
     //ComputeCellAveragedHorizontalGradient(vert->dwdy, 1, vert->wf, grid, prop, phys, myproc);
     //ComputeCellAveragedHorizontalGradient(vert->dwdx, 0, vert->wf, grid, prop, phys, myproc);
     //ComputeCellAveragedHorizontalGradient(vert->dvdy, 1, vert->vf, grid, prop, phys, myproc);
     //ComputeCellAveragedHorizontalGradient(vert->dudx, 0, vert->uf, grid, prop, phys, myproc);
   }
}

/*
 * Function: ComputeCellAveragedGradient
 * Compute the cell averaged gradient for scalars
 * ----------------------------------------------------
 * direction=0 x gradient, direction=1 y gradient
 */
void ComputeCellAveragedHorizontalGradient(REAL **gradient, int direction, REAL **scalar, gridT *grid, propT *prop, 
	physT *phys, int myproc)
{
  int i,k,nf,ne;
  REAL vec;

  for(i=0;i<grid->Nc;i++)
  {
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
    {
      gradient[i][k]=0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
        if(direction==0)
          vec=vert->n1[i*grid->maxfaces+nf];
        else
          vec=vert->n2[i*grid->maxfaces+nf];

        gradient[i][k]+=scalar[ne][k]*grid->df[ne]*vec;
      }

      gradient[i][k]/=grid->Ac[i];
    }
  }
}

// Computes the gradient from cell-centered quantities rather than face-centered, as in
// ComputeCellAveragedHorizontalGradient();
void ComputeCellAveragedHorizontalGradientCell(REAL **gradient, int direction, REAL **scalar,  int w_centered, gridT *grid, propT *prop, physT *phys, int myproc) {
  int i, iptr, k, nf, ne, neigh, nc1, nc2, Nkmax;
  REAL face_value, vec, xc, yc;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i=grid->cellp[iptr];

    if(w_centered)
      Nkmax=grid->Nk[i]+1;
    else
      Nkmax=grid->Nk[i];
    
    for(k=grid->ctop[i];k<Nkmax;k++) {
      gradient[i][k]=0;
      
      for(nf=0;nf<grid->nfaces[i];nf++) {
        ne = grid->face[i*grid->maxfaces+nf];
	nc1 = grid->grad[2*ne];
	nc2 = grid->grad[2*ne+1];

	if(nc1!=-1 && nc2!=-1)
	  face_value = 0.5*(scalar[nc1][k]+scalar[nc2][k]);
	else if(nc1==-1) {
	  if(nc2==-1) { printf("ERROR nc2!\n"); exit(1); }
	  face_value = scalar[nc2][k];
	} else {
	  if(nc1==-1) { printf("ERROR nc1!\n"); exit(1); }	  
	  face_value = scalar[nc1][k];
	}
	
        if(direction==0)
          vec=grid->n1[ne];
        else
          vec=grid->n2[ne];

        gradient[i][k]+=face_value*vec*grid->df[ne]*grid->normal[i*grid->maxfaces+nf];
      }
      gradient[i][k]/=grid->Ac[i];
    }
  }
}

/*
 * Function: ComputeNormalVector
 * Compute the outpointing normalized vector from cell center to each face center
 * ----------------------------------------------------
 * n1 is the x-component and n2 is y-component
 */
/*void ComputeNormalVector(gridT *grid, physT *phys, int myproc)
{ 
  int i,nf,ne;
  REAL D,Dx,Dy;
  for(i=0;i<grid->Nc;i++)
    for(nf=0;nf<grid->nfaces[i];nf++)
    {
      ne=grid->face[i*grid->maxfaces+nf];
      D=(grid->xv[i]-grid->xe[ne])*(grid->xv[i]-grid->xe[ne])+
      (grid->yv[i]-grid->ye[ne])*(grid->yv[i]-grid->ye[ne]);
      D=sqrt(D);
      Dx=grid->xe[ne]-grid->xv[i];
      Dy=grid->ye[ne]-grid->yv[i];
      vert->n1[i*grid->maxfaces+nf]=Dx/D;
      vert->n2[i*grid->maxfaces+nf]=Dy/D;
    }
}*/

/*
 * Function: ComputeNormalVector
 * Compute the outpointing normalized vector from cell center to each face center
 * ----------------------------------------------------
 * n1 is the x-component and n2 is y-component
 */
void ComputeNormalVector(gridT *grid, physT *phys,  int myproc, MPI_Comm comm) {
  int j, nc1, nc2, myperiodic=0, periodic;
  REAL nx1, nx2, ny1, ny2;
  

  //printf("hello from proc %d \n", myproc);

  for(j=0;j<grid->Ne;j++) {
    nc1=grid->grad[2*j];
    nc2=grid->grad[2*j+1];
    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1)
      nc2=nc1;
    
    nx1=grid->xe[j]-grid->xv[nc1];
    nx2=grid->xe[j]-grid->xv[nc2];
    ny1=grid->ye[j]-grid->yv[nc1];
    ny2=grid->ye[j]-grid->yv[nc2];

    // The grid is periodic if the cell centers are on the same side of the face
    if(nx1*nx2+ny1*ny2>0 && nc1!=nc2) {
      //printf("proc %d says periodic \n", myproc);
      myperiodic=1;
      break;
    }
  }

  //printf("before reduce on proc %d \n", myproc);

  MPI_Reduce(&myperiodic,&(periodic),1,MPI_INT,MPI_SUM,0,comm); //removed parenthesis around periodic (after &)
  
  //printf("after reduce on proc %d \n", myproc);
  
  MPI_Bcast(&periodic,1,MPI_INT,0,comm);
  //printf("periodic is %d on proc %d after bcast \n", periodic, myproc);

  //printf("through reduce and bcase \n");

  if(myproc==0 && VERBOSE>2) {
    if(periodic)
      printf("The grid is periodic!\n");
    else
      printf("The grid is not periodic!\n");
  }
  
  if(periodic)
    ComputeNormalVectorPeriodic(grid,phys,myproc);
  else
    ComputeNormalVectorNonPeriodic(grid,phys,myproc);
}
  
void ComputeNormalVectorPeriodic(gridT *grid, physT *phys, int myproc)
{ 
  int i, nf, ne, p1, p2;
  REAL Dx, Dy, D, xc, yc;

  for(i=0;i<grid->Nc;i++) {
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne=grid->face[i*grid->maxfaces+nf];

      p1 = grid->cells[i*grid->maxfaces+nf];
      if(nf==grid->nfaces[i]-1)
	p2 = grid->cells[i*grid->maxfaces];
      else
	p2 = grid->cells[i*grid->maxfaces+nf+1];

      xc = 0.5*(grid->xp[p1]+grid->xp[p2]);
      yc = 0.5*(grid->yp[p1]+grid->yp[p2]);      
      Dx = xc-grid->xv[i];
      Dy = yc-grid->yv[i];
      D = sqrt(Dx*Dx+Dy*Dy);

      vert->n1[i*grid->maxfaces+nf]=Dx/D;
      vert->n2[i*grid->maxfaces+nf]=Dy/D;
    }
    //printf("made it through periodic loop %d of %d \n", i, grid->Nc);
  }
}

void ComputeNormalVectorNonPeriodic(gridT *grid, physT *phys, int myproc)
{ 
  int i,nf,ne, p1,p2;
  REAL D,Dx,Dy,a,b,c,x,y,DD;
  for(i=0;i<grid->Nc;i++)
    for(nf=0;nf<grid->nfaces[i];nf++)
    {
      ne=grid->face[i*grid->maxfaces+nf];
      p1 = grid->edges[NUMEDGECOLUMNS*ne];
      p2 = grid->edges[NUMEDGECOLUMNS*ne+1];
      if(grid->xp[p2]!=grid->xp[p1]){
        a=(grid->yp[p2]-grid->yp[p1])/(grid->xp[p2]-grid->xp[p1]);
        b=-1;
        c=(grid->yp[p1]*grid->xp[p2]-grid->yp[p2]*grid->xp[p1])/(grid->xp[p2]-grid->xp[p1]);
      } else {
        a=1;
        b=0;
        c=-grid->xp[p1];
      }
      // calculate orthogonal distance
      D=fabs(grid->xv[i]*a+b*grid->yv[i]+c)/sqrt(a*a+b*b);
      // calculate x and y intersection point
      x=(b*(b*grid->xv[i]-a*grid->yv[i])-a*c)/(a*a+b*b);
      y=(a*(-b*grid->xv[i]+a*grid->yv[i])-b*c)/(a*a+b*b);
      Dx=x-grid->xv[i];
      Dy=y-grid->yv[i];
      vert->n1[i*grid->maxfaces+nf]=Dx/D;
      vert->n2[i*grid->maxfaces+nf]=Dy/D;
    }
}

/* 
 * Function: OpenVertCoordinateFiles
 * Usage: OpenFiles(myproc);
 * ------------------------------
 * Open all of the files used for i/o to store the file pointers.
 * include dzz, zc, omega for each cell at different time steps
 * 
 */
void OpenVertCoordinateFiles(gridT *grid,int mergeArrays, int myproc)
{
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];
  MPI_GetFile(filename,DATAFILE,"zcFile","OpenFiles",myproc);
  if(mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  vert->zcFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"dzzFile","OpenFiles",myproc);
  if(mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  vert->dzzFID = MPI_FOpen(str,"w","OpenFiles",myproc);

  MPI_GetFile(filename,DATAFILE,"omegaFile","OpenFiles",myproc);
  if(mergeArrays)
    strcpy(str,filename);
  else
    sprintf(str,"%s.%d",filename,myproc);
  vert->omegaFID = MPI_FOpen(str,"w","OpenFiles",myproc);
}

/*
 * Function: OutputVertcoordinate
 * Usage: OutputVertcoordinate(grid,phys,prop,myproc,numprocs,comm);
 * ---------------------------------------------------------------------------
 * Output the data every ntout steps as specified in suntans.dat
 * now include zc,dzz,omega
 *
 */
void OutputVertCoordinate(gridT *grid, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, j, jptr, k, nwritten, arraySize, writeProc,nc1,nc2;
  char str[BUFFERLENGTH], filename[BUFFERLENGTH];

  if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart) { 
    Write3DData(vert->zl,vert->tmp,prop->mergeArrays,vert->zcFID,
    "Error outputting uc-data!\n",grid,numprocs,myproc,comm);
    Write3DData(grid->dzz,vert->tmp,prop->mergeArrays,vert->dzzFID,
    "Error outputting uc-data!\n",grid,numprocs,myproc,comm);
    Write3DData(vert->omega,vert->tmp,prop->mergeArrays,vert->omegaFID,
    "Error outputting uc-data!\n",grid,numprocs,myproc,comm);
  }

  if(prop->n==prop->nsteps+prop->nstart && myproc==0) {
    fclose(vert->zcFID);
    fclose(vert->dzzFID);
    fclose(vert->omegaFID);
  }
}

/*
 * Function: VertCoordinateBasic
 * Setup the basics and initial condition for the vertical coordinate
 * ----------------------------------------------------
 * 
 */
void VertCoordinateBasic(gridT *grid, propT *prop, physT *phys, int myproc, MPI_Comm comm)
{
  int i,k,j;
  // allocate subgrid struture first
  vert=(vertT *)SunMalloc(sizeof(vertT),"VertCoordinateBasic");

  // allocate necessary variable
  AllocateandInitializeVertCoordinate(grid,prop, myproc);

  // output
  OpenVertCoordinateFiles(grid,prop->mergeArrays,myproc);

  // if not z-level setup dzz
  InitializeLayerThickness(grid, prop, phys,comm,myproc);
  
  // compute zc
  ComputeZc(grid,prop,phys,1,myproc);
  
  // compute normal vector
  ComputeNormalVector(grid,phys,myproc,comm);
  
  //  ComputeCellAveragedHorizontalGradient(vert->dzdx, 0, vert->zf, grid, prop, phys, myproc);
  //  ComputeCellAveragedHorizontalGradient(vert->dzdy, 1, vert->zf, grid, prop, phys, myproc);
  ComputeCellAveragedHorizontalGradientCell(vert->dzdx, 0, vert->zl, 1, grid, prop, phys, myproc);
  ComputeCellAveragedHorizontalGradientCell(vert->dzdy, 1, vert->zl, 1, grid, prop, phys, myproc);
  /*
  k = grid->Nk[i]-1;
  for(i=0;i<grid->Nc;i++)
    printf("%f %f %f\n",grid->xv[i],vert->dzdx[i][k],vert->dzdy[i][k]);


  exit(0);
  */
  //printf("exited grads \n");

  //  ComputeCellAveragedHorizontalGradientCell(vert->dzdx, 0, vert->zc, grid, prop, phys, myproc);
  //  ComputeCellAveragedHorizontalGradientCell(vert->dzdy, 1, vert->zc, grid, prop, phys, myproc);  

  // compute nodal information
  // prepare for the calculation of nodal relative vorticity
  ComputeNodalData(grid,myproc);

  //printf("exited nodal  \n");
}

/*
 * Function: VertCoordinateBasicRestart
 * Setup the basics and initial condition for the vertical coordinate
 * for restarting run
 * ----------------------------------------------------
 * 
 */
void VertCoordinateBasicRestart(gridT *grid, propT *prop, physT *phys, int myproc, MPI_Comm comm)
{
  int i,k,j;
  // allocate subgrid struture first
  vert=(vertT *)SunMalloc(sizeof(vertT),"VertCoordinateBasic");

  // allocate necessary variable
  AllocateandInitializeVertCoordinate(grid,prop, myproc);

  // output 
  OpenVertCoordinateFiles(grid,prop->mergeArrays,myproc);

  // compute zc
  ComputeZc(grid,prop,phys,0,myproc);

  // sigma coordinate need to recalculate dsigma for each layer
  if(prop->vertcoord==3)
    ComputeDSigma(grid,phys,myproc);

  // compute normal vector
  ComputeNormalVector(grid,phys,myproc,comm);

  //  ComputeCellAveragedHorizontalGradient(vert->dzdx, 0, vert->zf, grid, prop, phys, myproc);
  //  ComputeCellAveragedHorizontalGradient(vert->dzdy, 1, vert->zf, grid, prop, phys, myproc);
  ComputeCellAveragedHorizontalGradientCell(vert->dzdx, 0, vert->zl, 1, grid, prop, phys, myproc);
  ComputeCellAveragedHorizontalGradientCell(vert->dzdy, 1, vert->zl, 1, grid, prop, phys, myproc);  

  // compute nodal information
  // prepare for the calculation of nodal relative vorticity
  ComputeNodalData(grid,myproc);
}


/*
 * Function: ComputeDSigma
 * based on the dzz to calculate dsigma for restart function
 * ----------------------------------------------------
 * dsigma[k]=dzz[k]/(h+d)
 */
void ComputeDSigma(gridT *grid, physT *phys, int myproc){
  int i=0,k;
  for(k=0;k<grid->Nk[i];k++)
    vert->dsigma[k]=grid->dzz[i][k]/(phys->h[i]+grid->dv[i]+phys->zB[i]);
}

/*
 * Function: StoreVariables
 * Usage: StoreVariables(grid,phys);
 * ---------------------------------
 * Store the old values of s, u, and w into stmp3, u_old, and wtmp2,
 * respectively.
 *
 */
void StoreVertVariables(gridT *grid, physT *phys) {
  int i, j, k, iptr, jptr;

  for(i=0;i<grid->Nc;i++) 
    for(k=0;k<grid->Nk[i];k++) {
      // store omega^n-1 and omega^n
      vert->omega_old2[i][k]=vert->omega_old[i][k];
      vert->U3_old2[i][k]=vert->U3_old[i][k];
      vert->omega_old[i][k]=vert->omega[i][k];
      vert->U3_old[i][k]=vert->U3[i][k];
    }
}

/*
 * Function: UpdateCellCenteredFreeSurface
 * Usage: UpdateCellCenteredFreeSurface(grid,prop,phys,myproc);
 * ---------------------------------
 * update the free surface height after the nonhydrostatic pressure is solved
 *
 */
void UpdateCellcenteredFreeSurface(gridT *grid, propT *prop, physT *phys, int myproc)
{
   int i,iptr,k,nf,ne;
   REAL sum,normal,tmp;
   for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    sum = 0;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];
      for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
      {  
        sum+=(prop->imfac1*phys->u[ne][k]+prop->imfac2*phys->u_old[ne][k]+prop->imfac3*phys->u_old2[ne][k])*
          grid->df[ne]*normal*grid->dzf[ne][k];
      }
    }
    tmp=phys->h_old[i]-sum*prop->dt/grid->Ac[i];
    phys->h[i]=tmp;
  }
}

/*
 * Function: UpdateFluxHeight
 * Usage: UpdateFluxHeight(grid,prop,phys,myproc);
 * ---------------------------------
 * Update the flux height use the upwind the scheme with u*, u^n-1, u^n-2
 * The reason is to ensure a positive layer height for the isopycnal coordinate
 * The flux height set by the SetFluxHeight is calculated by u^n-1 which need to be verified by
 * u^im which is used to calculate the flux into/out of a grid cell.
 * The new flux height will be used to further modified free surface height/ layer thickness
 * only work when vertcoord==2
 */
void SetOldFluxHeight(REAL **dzf_old, REAL **dzf_old2, REAL **dzz_t2, gridT *grid, physT *phys, propT *prop, MPI_Comm comm, int myproc) {
  int j, jptr, k, nc1, nc2;

  // Set flux height at boundaries as the height of interior point
  for(jptr=grid->edgedist[1];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];

    for(k=0;k<grid->Nkmax;k++){
      dzf_old[j][k]=grid->dzz[nc1][k];
      dzf_old2[j][k]=grid->dzz[nc1][k];
    }
  }

  if(prop->vertcoord==2 || prop->vertcoord==4) {
    //TvdFluxHeight(grid, phys, prop, dzfmeth, comm, myproc);  
    int iptr, i, ib;
    int num_layers_bot = prop->botLayers; //40, 4, 30 for most sand wave runs. 40 for Nk=90 sandwaves
  
    // give values for the boundary flux height
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[4];jptr++) {
        j=grid->edgep[jptr];
        ib=grid->grad[2*j];
  
        for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
          phys->boundary_tmp[jptr-grid->edgedist[2]][k]=grid->dzzold[ib][k];
        }
    }
    // calculate tvd face values
    // will put directly into dzf_old2
    HorizontalFaceScalars_oldstep(grid,phys,prop,grid->dzzold,dzf_old,phys->boundary_tmp,vert->dzfmeth,1,2,comm,myproc);


    // now for two steps ago
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[4];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
        phys->boundary_tmp[jptr-grid->edgedist[2]][k]=dzz_t2[ib][k];
      }
  }
    // calculate tvd face values
    // will put directly into dzf_old2
   
    //consider changing this?
    HorizontalFaceScalars_oldstep(grid,phys,prop,dzz_t2,dzf_old2,phys->boundary_tmp,vert->dzfmeth,2,2,comm,myproc);

    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr];
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
    
    for(k=grid->Nkmax-num_layers_bot;k<grid->Nkmax;k++){
      dzf_old[j][k]=0.5*(grid->dzzold[nc1][k]+grid->dzzold[nc2][k]);
      dzf_old2[j][k]=0.5*(dzz_t2[nc1][k]+dzz_t2[nc2][k]);
    }
  }

} else {
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    
    for(k=0;k<grid->Nkmax;k++) {
      dzf_old[j][k]=0.5*(grid->dzzold[nc1][k]+grid->dzzold[nc2][k]);
      dzf_old2[j][k]=0.5*(dzz_t2[nc1][k]+dzz_t2[nc2][k]);
    }
  }
}

/*
for(j=0;j<grid->Ne;j++) {
  for(k=0;k<grid->Nkmax;k++) {
    if(grid->dzf[j][k]==0 && grid->mark[j]!=1)
printf("DZF = 0 at k=%d, marker=%d\n",k,grid->mark[j]);
  }
}
*/
}


/*
 * Function: UpdateFluxHeight
 * Usage: UpdateFluxHeight(grid,prop,phys,myproc);
 * ---------------------------------
 * Update the flux height use the upwind the scheme with u*, u^n-1, u^n-2
 * The reason is to ensure a positive layer height for the isopycnal coordinate
 * The flux height set by the SetFluxHeight is calculated by u^n-1 which need to be verified by
 * u^im which is used to calculate the flux into/out of a grid cell.
 * The new flux height will be used to further modified free surface height/ layer thickness
 * only work when vertcoord==2
 */
void VerifyFluxHeight(gridT *grid, propT *prop, physT *phys, int myproc)
{  
  int i, j, k, nc1, nc2;
  REAL dz_bottom, dzsmall=grid->dzsmall,h_uw,utmp,tmp;

  for(j=0;j<grid->Ne;j++) 
  {
    grid->hf[j]=0;
    //for(k=0; k<grid->Nkc[j];k++)
      //grid->dzf[j][k]=0;
  }

  for(j=0;j<grid->Ne;j++) {
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    if(nc1==-1) nc1=nc2;
    if(nc2==-1) nc2=nc1;

    for(k=0;k<grid->etop[j];k++)
      grid->dzf[j][k]=0;
 
    for(k=grid->etop[j];k<grid->Nke[j];k++){
      utmp=prop->imfac1*phys->u[j][k]+prop->imfac2*phys->u_old[j][k]+prop->imfac3*phys->u_old2[j][k];
      tmp=UpWind(utmp,grid->dzz[nc1][k],grid->dzz[nc2][k]);
      grid->dzf[j][k]=tmp;
    }

    /* This works with Wet/dry but not with cylinder case...*/
    if(grid->etop[j]==grid->Nke[j]-1) {
      // added part
     if(grid->mark[j]==2 && grid->dzf[j][k]<=0.01)
        grid->dzf[j][k]=0.01;
    }

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      if(grid->dzf[j][k]<=DRYCELLHEIGHT)
        grid->dzf[j][k]=0;
      
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      grid->hf[j]+=grid->dzf[j][k];
  }
}

/*
 * Function: TvdFluxHeight
 * Usage: TvdFluxHeight(gridT *grid, physT *phys, propT *prop, int TVD, MPI_Comm comm, int myproc)
 * ---------------------------------
 * calculate flux height use tvd scheme
 */
void TvdFluxHeight(gridT *grid, physT *phys, propT *prop, int TVD, MPI_Comm comm, int myproc)
{
  int jptr, j, iptr, i, ib, k;
  REAL z;

  // give values for the boundary flux height
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[4];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
        phys->boundary_tmp[jptr-grid->edgedist[2]][k]=grid->dzz[ib][k];
      }
  }

  // calculate tvd face values
  HorizontalFaceScalars(grid,phys,prop,grid->dzz,phys->u,phys->boundary_tmp,vert->dzfmeth,2,comm,myproc);
}

/*
 * Function: ComputeNodalData
 * Usage: ComputeNodalData(gridT *grid,int myproc)
 * ---------------------------------
 * prepare for the calculation of nodal relative vorticity
 * 1. get the type of each node vert->typep
 * 2. get the Ap for each node (type=1)
 * 3. get the Nccp for each edge
 * 4. get the Nkpmin for each node
 */
void ComputeNodalData(gridT *grid,int myproc)
{
  int j,jptr,p1,p2,nc1,nc2,tmp;
  
  // node type
  for(j=0;j<grid->Ne;j++)
  {
    p1=grid->edges[NUMEDGECOLUMNS*j];
    p2=grid->edges[NUMEDGECOLUMNS*j+1];
    // if not computational edge
    if(grid->mark[j]!=0 && grid->mark[j]!=5)
    {
      vert->typep[p1]=-1;
      vert->typep[p2]=-1;
    }
  }

  // computational edge CCNpe
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) 
  {
    j = grid->edgep[jptr];
    nc1=grid->grad[2*j];
    nc2=grid->grad[2*j+1];
    p1=grid->edges[NUMEDGECOLUMNS*j];
    p2=grid->edges[NUMEDGECOLUMNS*j+1];

    // node p1
    // the counter-clockwise direction of a node is (-(ye-yp),xe-xp)
    // the positive direction of a normal velocity is (xc1-xc2,yc1-yc2)
    // the inner production of p2 should be opposite to p1
    tmp=-(grid->xv[nc1]-grid->xv[nc2])*(grid->ye[j]-grid->yp[p1])+
    (grid->yv[nc1]-grid->yv[nc2])*(grid->xe[j]-grid->xp[p1]);
    if(tmp>0){
      vert->CCNpe[2*j]=1.0;
      vert->CCNpe[2*j+1]=-1.0;
    }
    else{
      vert->CCNpe[2*j]=-1.0;
      vert->CCNpe[2*j+1]=1.0;
    }
  }

  // node Ap and Nkpmin
  for(j=0;j<grid->Ne;j++)
  {
    p1=grid->edges[NUMEDGECOLUMNS*j];
    p2=grid->edges[NUMEDGECOLUMNS*j+1];
    
    if(vert->typep[p1]==1)
    {
      vert->Ap[p1]+=grid->dg[j]*grid->df[j]/4;
      if(grid->Nke[j]<vert->Nkpmin[p1])
        vert->Nkpmin[p1]=grid->Nke[j];
    }
    
    if(vert->typep[p2]==1)
    {
      vert->Ap[p2]+=grid->dg[j]*grid->df[j]/4;
      if(grid->Nke[j]<vert->Nkpmin[p2])
        vert->Nkpmin[p2]=grid->Nke[j];
    }
  }
}

/*
 * Function: ComputeRelativeVorticity
 * Usage: ComputeRelativeVorticity(gridT *grid, physT *phys, propT *prop, int myproc)
 * ---------------------------------
 * compute all relative vorticity
 * 1. cell-centered
 * 2. nodal 
 * 3. edge-center
 */
void ComputeRelativeVorticity(gridT *grid, physT *phys, propT *prop, int myproc)
{
  int i, j, k,p1,p2,Nkmin,jptr;

  // cell centered relative vorticity
  for(i=0;i<grid->Nc;i++)
    for(k=0;k<grid->Nk[i];k++)
      vert->f_r[i][k]=vert->dvdx[i][k]-vert->dudy[i][k];

  // nodal relative vorticity
  for(i=0;i<grid->Np;i++)
    for(k=0;k<grid->Nkp[i];k++)
      vert->f_rp[i][k]=0;

  for(j=0;j<grid->Ne;j++)
  {
    p1=grid->edges[NUMEDGECOLUMNS*j];
    p2=grid->edges[NUMEDGECOLUMNS*j+1];
    
    if(vert->typep[p1]==1)
      for(k=0;k<vert->Nkpmin[p1];k++)
        vert->f_rp[p1][k]+=grid->dg[j]*vert->CCNpe[2*j]*phys->u[j][k]/vert->Ap[p1];
    
    if(vert->typep[p2]==1)
      for(k=0;k<vert->Nkpmin[p2];k++)
        vert->f_rp[p2][k]+=grid->dg[j]*vert->CCNpe[2*j+1]*phys->u[j][k]/vert->Ap[p2];
  }

  // edge-centered velocity
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) 
  {
    j = grid->edgep[jptr];
    p1=grid->edges[NUMEDGECOLUMNS*j];
    p2=grid->edges[NUMEDGECOLUMNS*j+1];
    Nkmin=vert->Nkpmin[p1];
    if(Nkmin>vert->Nkpmin[p2])
      Nkmin=vert->Nkpmin[p2];

    // cell center to edge center
    for(k=grid->etop[j];k<grid->Nke[j];k++)
      vert->f_re[j][k]=InterpToFace(j,k,vert->f_r,phys->u,grid);

    // only use nodal for all nodes are type 1
    if(vert->typep[p1]==1 && vert->typep[p2]==1)
      for(k=grid->etop[j];k<Nkmin;k++)
        vert->f_re[j][k]=0.5*(vert->f_rp[p1][k]+vert->f_rp[p2][k]);
  }
}


void Generate_dzzolds_varationalvoid(gridT *grid, propT *prop, physT *phys)
{
  int i,k,j,nf,ne,dry,max_ind,neigh,num_bot_layers;
  REAL fac1,fac2,fac3,Ac,sum,sum_old,maxdzz,u_im,alpha_f,error,dzf;

  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;
 
  num_bot_layers = 0; //40
  REAL num_mid_layers = 0; //15
  REAL cutoffpoint;
    for(i=0;i<grid->Nc;i++) {
      for(k=grid->ctop[i];k<grid->Nk[i]-(num_bot_layers+num_mid_layers);k++) {
        for(nf=0;nf<grid->nfaces[i];nf++) {               
          ne = grid->face[i*grid->maxfaces+nf];

          u_im = phys->u[ne][k]; //n 

    //need old dzf      
    grid->dzzold[i][k]-=prop->dt*u_im*grid->dzf[ne][k]*
      grid->normal[i*grid->maxfaces+nf]*grid->df[ne]/grid->Ac[i];
        }

  // Fix the layer heights if smaller than vertdzmin
  if(grid->dzzold[i][k]<vert->vertdzmin)
    grid->dzzold[i][k]=vert->vertdzmin;
  }
    }

  // Hybrid layers 
  if(num_bot_layers){
      //update bottom layers
      cutoffpoint = 280; //270
      for(i=0;i<grid->Nc;i++){
      sum=0;
      for(k=0;k<grid->Nk[i]-(num_bot_layers+num_mid_layers);k++){
        sum+=grid->dzzold[i][k];
      }
      for(k=grid->Nk[i]-(num_bot_layers+num_mid_layers);k<(grid->Nk[i]-num_bot_layers);k++){
        grid->dzzold[i][k]=(1.0/num_mid_layers)*(cutoffpoint-sum+phys->h[i]);
      }
      for(k=grid->Nk[i]-(num_bot_layers);k<grid->Nk[i];k++){
        //grid->dzz[i][k]=(1.0/num_bot_layers)*(grid->dv[i]-cutoffpoint-phys->zB[i]);
        grid->dzzold[i][k]=grid->dzz[i][k]; //stay same? 
      }
      }

  }
}

