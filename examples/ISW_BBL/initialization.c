#include "math.h"
#include "fileio.h"
#include "suntans.h"
#include "initialization.h"
#include "sediments.h"
#include "wave.h"
#include "culvert.h"
#define sech 1/cosh

REAL myrand(void);

/*
 * Function: GetDZ
 * Usage: GetDZ(dz,depth,Nkmax,myproc);
 * ------------------------------------
 * Returns the vertical grid spacing in the array dz.
 *
 */
int GetDZ(REAL *dz, REAL depth, REAL localdepth, int Nkmax, int myproc) {
  int k, status,Nk_noiso=10;
  REAL z=0, dz0, r = GetValue(DATAFILE,"rstretch",&status),sum=0;
  REAL a0=120;
  REAL runsum=0, mindz, dzconst;
  //printf("getdz has been called \n");

  if(dz!=NULL) {
    if(r==1){
        for(k=0;k<Nkmax;k++)
          dz[k]=depth/Nkmax;
    }else if(r>1 && r<=1.1) {    
      dz[0] = depth*(r-1)/(pow(r,Nkmax)-1);
      if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[0]);
      for(k=1;k<Nkmax;k++) 
  dz[k]=r*dz[k-1];
    } else if(r>-1.1 && r<-1) {    
      
      //normal stretching from bottom
      //r=fabs(r);
      //dz[Nkmax-1] = depth*(r-1)/(pow(r,Nkmax)-1);
      //printf("dz for k = %d is %f \n", k, dz[Nkmax-1]);
      //if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[Nkmax-1]);
      //for(k=Nkmax-2;k>=0;k--){ 
      //dz[k]=r*dz[k+1];
	//printf("dz for k = %d is %f \n", k, dz[k]);
	//}
      
      
      //stretching up to a constant
      r=fabs(r);
      mindz = 0.2; //minimum res at bottom
      dzconst = 1.9985; //max res

      //mindz=0.1;

      mindz = 0.2; //minimum res at bottom                                                                                                                                                                       
      dzconst = 1.9985; //max res   

      dz[Nkmax-1] = mindz;
      sum = mindz;
      for(k=Nkmax-2;k>=0;k--) {
        dz[k]=r*dz[k+1];
        if(dz[k]>dzconst){
          dz[k]=dzconst;
        }
        sum+=dz[k];
      }

      //adjust every layer to make exactly sum to depth 

      REAL alpha = depth/sum;
      for(k=Nkmax-1;k>=0;k--) {
        dz[k]*=alpha;
        //printf("dz[k] is %f for k = %d \n", dz[k], k);
      }
      
  //   } else if(r>-1.1 && r<-1) {    
  //     r=fabs(r);
  //     dz[Nkmax/2-1] = (depth/4)*(r-1)/(pow(r,Nkmax/2)-1); //only stretch grid up to 1/4 of depth                                                            
  //     runsum=dz[Nkmax-1];
  //     if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[Nkmax-1]);
  //     for(k=Nkmax-2;k>=0;k--) {
  //       dz[k]=r*dz[k+1];
	// if(runsum+dz[k]>depth/4){
	//   dz[k]=dz[k+1];}//no more stretching
  //       runsum+=dz[k];
  //     }
  //     /*dz[Nkmax-1] = depth*(r-1)/(pow(r,Nkmax)-1);
  //     //if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[Nkmax-1]);
  //     //for(k=Nkmax-2;k>=0;k--) 
  //     dz[k]=r*dz[k+1];
  //     */
  //   } 
    } else if(r==2) {
      //use as marker for tanh grid, values found in matlab for D=300; 
      int numlayers_top = 173;
      REAL d_0 = 10.1;
      REAL dztop = 1.5;
      REAL dzbot = 0.1;
      REAL zcenter = -280;
      REAL tanhlayers = Nkmax-numlayers_top;
      REAL dzval;

      z=0;

      z = -300;
      for(k=Nkmax-1;k>numlayers_top-1;k--){
	dzval=dzbot+(dztop-dzbot)*0.5*(1 + tanh((z-zcenter)/d_0));
	dz[k]=roundf(dzval * 1000) / 1000; 
	z+=dz[k];
	//printf("dz = %f at z = %f and k = %d \n", dz[k], z, k);
      }
      
      
      //add a bit to the last three values in the tanh so we get to 300
      REAL remainder = numlayers_top*dztop+z;
      z=z-dz[k+3]-dz[k+1]-dz[k+2];
      dz[k+3]=dz[k+3]-remainder/3;
      dz[k+1]=dz[k+1]-remainder/3;
      dz[k+2]=dz[k+2]-remainder/3;
      z=z+dz[k+3]+dz[k+1]+dz[k+2];

      //printf("z = %f and k = %d \n", z, k);

      for(k=numlayers_top-1;k>=0;k--){
	dz[k]=dztop;
        z+=dz[k];
	//printf("dz = %f at z = %f and k = %d \n", dz[k], z, k);
      }
 
    } else {
      printf("Error in GetDZ when trying to create vertical grid:\n");
      printf("Absolute value of stretching parameter rstretch must  be in the range (1,1.1).\n");
      exit(1);
    }
  } else {
    r=fabs(r);
    if(r!=1)
      dz0 = depth*(r-1)/(pow(r,Nkmax)-1);
    else
      dz0 = depth/Nkmax;
    z = dz0;
    for(k=1;k<Nkmax;k++) {
      dz0*=r;
      z+=dz0;
      if(z>=localdepth) {
	return k;
      }
    }
  }
}
  
/*
 * Function: ReturnDepth
 * Usage: grid->dv[n]=ReturnDepth(grid->xv[n],grid->yv[n]);
 * --------------------------------------------------------
 * Helper function to create a bottom bathymetry.  Used in
 * grid.c in the GetDepth function when IntDepth is 0.
 *
 */
REAL ReturnDepth(REAL x, REAL y) {
  //return 100;
  return 300;
  //REAL L = 5e3, k = 2*PI/(L/4), a = 3;
  //return 300 - a*cos(k*x);
}
 
 
/*
  * Function: ReturnFreeSurface
  * Usage: grid->h[n]=ReturnFreeSurface(grid->xv[n],grid->yv[n]);
  * -------------------------------------------------------------
  * Helper function to create an initial free-surface. Used
  * in phys.c in the InitializePhysicalVariables function.
  *
  */
REAL ReturnFreeSurface(REAL x, REAL y, REAL d) {
  return 0;
}

/*
 * Function: ReturnSalinity
 * Usage: grid->s[n]=ReturnSalinity(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial salinity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnSalinity(REAL x, REAL y, REAL z) {
  return 0;
}

/*
 * Function: IsoReturnSalinity
 * Usage: grid->T[n]=IsoReturnSalinity(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial salinity field under iso 
 * pycnal coordinate. Used in phys.c in the 
 * InitializePhysicalVariables function.
 *
 */
REAL IsoReturnSalinity(REAL x, REAL y, REAL z, REAL zu, REAL zl, int i, int k) {
  int n, N=100;
  REAL s, sum_sdz=0, h=zu-zl, zint, dz=10, zc = 0.5*(zu+zl), dz_end, sum_dz;

  return ReturnSalinity(x,y,zc);
}

/*
 * Function: IsoReturnTemperature
 * Usage: grid->T[n]=IsoReturnTemperature(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial temperature field under iso 
 * pycnal coordinate. Used in phys.c in the 
 * InitializePhysicalVariables function.
 *
 */
REAL IsoReturnTemperature(REAL x, REAL y, REAL z, REAL depth, int i, int k) {
  return ReturnTemperature(x,y,z,depth);
}

/*
 * Function: ReturnTemperature
 * Usage: grid->T[n]=ReturnTemperaturegrid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial temperature field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnTemperature(REAL x, REAL y, REAL z, REAL depth) {
  return x/10000;
}

/*
 * Function: ReturnHorizontalVelocity
 * Usage: grid->u[n]=ReturnHorizontalVelocity(grid->xv[n],grid->yv[n],
 *                                            grid->n1[n],grid->n2[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial velocity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnHorizontalVelocity(REAL x, REAL y, REAL n1, REAL n2, REAL z) {
  // FILE *fid;
  
  // if((int)MPI_GetValue(DATAFILE,"init_u_from_file","ReadGrid",myproc)) {
  //   MPI_GetFile(filename,DATAFILE,"u_init_file","InitializePhysicalVariables",myproc);  
  //   //sprintf(str,"%s.%d",filename,myproc);  
  //   fid = MPI_FOpen(str,"r","InitializePhysicalVariables",myproc);

  //   for(j=0;j<grid->Ne;j++) {
  //     fread(u,sizeof(REAL),grid->Nkmax,fid);
  //   }
  //   fclose(fid);
  // }

  
  return 0.001*n1*sin(2*PI*x/10000);
}

/*
 * Function: ReturnSediment
 * Usage: SediC[Nsize][n][Nk]=ReturnSediment(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial sediment concentration field.  Used
 * in sediment.c IntitalizeSediment function
 *
 */
REAL ReturnSediment(REAL x, REAL y, REAL z, int sizeno) {
  return 0;
}

/*
 * Function: ReturnBedSedimentRatio
 * Usage: SediC[Nsize][n][Nk]=ReturnBedSedimentRatio(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial bed sediment concentration field.  Used
 * in sediment.c IntitalizeSediment function
 * the sum of ratio should be 1
 */
REAL ReturnBedSedimentRatio(REAL x, REAL y, int layer, int sizeno, int nsize) {
  return 1.00/nsize;
}

/*
 * Function: ReturnWindSpeed 
 * Usage: Uwind[n]=ReturnWindSpeed(grid->xv[n],grid->yv[n]);
 * ------------------------------------------------------------
 * Helper function to create an initial wind velocity field.  Used
 * in wave.c InitializeWave
 *
 */
REAL ReturnWindSpeed(REAL x, REAL y) {
  return 0;
}

/*
 * Function: ReturnWindDirection
 * Usage: Winddir[n]=ReturnWindDirection(grid->xv[n],grid->yv[n]);
 * ------------------------------------------------------------
 * Helper function to create an initial wind velocity field.  Used
 * in wave.c InitializeWave
 *
 */
REAL ReturnWindDirection(REAL x, REAL y) {
  return 0;
}

/*
 * Function: ReturnCulvertTop
 * usage: Culvertheight[i]=ReturnCulvertTop(grid->xv[i],grid->yv[i],myproc);
 * -------------------------------------------------------------------------------
 * provide the culvert size for culvert cell, for non culvert cell, assume Culvertheight[i]=EMPTY
 *
 */
REAL ReturnCulvertTop(REAL x, REAL y, REAL d)
{ 
  REAL sum;
  if(x>350 && x<650)
    return -1.2;
  else
    return INFTY;
}

/*
 * Function: ReturnMarshHeight
 * Usage: hmarsh[n]=ReturnMarshHeight(grid->xe[n],grid->ye[n]);
 * ------------------------------------------------------------
 * Helper function to create an initial hmarsh field.  Used
 * in marsh.c Interpmarsh
 *
 */
REAL ReturnMarshHeight(REAL x, REAL y)
{
  return 0;
}
/*
 * Function: ReturnMarshHeight
 * Usage: CdV[n]=ReturnMarshDragCoefficient(grid->xe[n],grid->ye[n]);
 * ------------------------------------------------------------
 * Helper function to create an initial CdV field.  Used
 * in marsh.c Interpmarsh
 *
 */
REAL ReturnMarshDragCoefficient(REAL x, REAL y)
{
  return 0.0;
}
/*
 * Function: ReturnSubgridPointDepth
 * Usage: Subgrid->dp[n]=ReturnSubgridPointDepth(subgrid->xp[n],subgrid->yp[n]);
 * ------------------------------------------------------------
 * Helper function to give the depth of each point for subgrid method
 *
 */
REAL ReturnSubgridPointDepth(REAL x, REAL y, REAL xv, REAL yv)
{
  return 1;
}

/*
 * Function: ReturnSubgridPointDepth
 * Usage: Subgrid->dp[n]=ReturnSubgridPointDepth(subgrid->xp[n],subgrid->yp[n]);
 * ------------------------------------------------------------
 * Helper function to give the depth of each point for subgrid method
 *
 */
REAL ReturnSubgridPointeDepth(REAL x, REAL y)
{
  return 1;
}



/*
 * Function: ReturnSubCellArea
 * Usage: calculate area for each subcell for each cell 
 * ------------------------------------------------------------
 * Helper function to give the area of each subcell for subgrid method
 *
 */
REAL ReturnSubCellArea(REAL x1, REAL y1, REAL x2, REAL y2, REAL x3, REAL y3, REAL h)
{
  return 0;
}

/*
 * Function: ReturnFluxHeight
 * Usage: calculate flux height for each subedge 
 * ------------------------------------------------------------
 * Helper function to give flux height of each subedge for subgrid method
 *
 */
REAL ReturnFluxHeight(REAL x1,REAL y1, REAL x2, REAL y2, REAL h)
{
  return 0;
}

/*
 * Function: ReturnSubgridErosionParameterizationEpslon
 * Usage: Subgrid->dp[n]=ReturnSubgridPointDepth(subgrid->xp[n],subgrid->yp[n]);
 * ------------------------------------------------------------
 * give the value of epslon when erosion parameterization is on
 *
 */
REAL ReturnSubgridErosionParameterizationEpslon(REAL x, REAL y)
{
  return 0;
}

REAL myrand(void) {
  return 2.0*(REAL)rand()/(REAL)RAND_MAX-1.0;
}
