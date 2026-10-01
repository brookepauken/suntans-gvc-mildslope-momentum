/*
 * File: check.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains functions that check for instability or
 * wetting and drying problems.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "check.h"
#include "suntans.h"
#include "stdio.h"
#include "timer.h"
#include "memory.h"
#include "util.h"

#define DASHES "----------------------------------------------------------------------\n"
#define CMAXSUGGEST 0.5

/*
 * Function: Check
 * Usage: Check(grid,phys,prop,myproc,numprocs,comm);
 * --------------------------------------------------
 * Check to make sure the run isn't blowing up.
 *
 */    
int Check(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, k, icu, kcu, icw, kcw, Nc=grid->Nc, Ne=grid->Ne, ih, is, ks, iu, ku, iw, kw, nc1, nc2,dry;
  int uflag=1, wflag=1, sflag=1, hflag=1, myalldone, alldone, progout;
  REAL C, CmaxU, CmaxW, allCmaxU, allCmaxW, dtsuggestU, dtsuggestW;

  icu=kcu=icw=kcw=ih=is=ks=iu=ku=iw=kw=0;

  for(i=0;i<Nc;i++) 
    if(phys->h[i]!=phys->h[i]) {
      hflag=0;
      ih=i;
      break;
    }

  for(i=0;i<Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      if(phys->s[i][k]!=phys->s[i][k]) {
        sflag=0;
        is=i;
        ks=k;
        break;
      }
    if(!sflag)
      break;
  }

  for(i=0;i<Ne;i++) {
    for(k=0;k<grid->Nke[i];k++)
      if(phys->u[i][k]!=phys->u[i][k]) {
        uflag=0;
        iu=i;
        ku=k;
        break;
      }
    if(!uflag)
      break;
  }

  for(i=0;i<Nc;i++) {
    for(k=0;k<grid->Nk[i];k++)
      if(phys->w[i][k]!=phys->w[i][k]) {
        wflag=0;
        iw=i;
        kw=k;
        break;
      }
    if(!wflag)
      break;
  }

  CmaxU=0;
  for(i=0;i<Ne;i++) 
    for(k=grid->etop[i];k<grid->Nke[i];k++) {
      C = fabs(phys->u[i][k])*prop->dt/grid->dg[i];
      if(C>CmaxU) {
        icu = i;
        kcu = k;
        CmaxU = C;
      }
    }

  /*CmaxU=0;
  for(i=0;i<Ne;i++) 
  {
    dry=0;
    nc1=grid->grad[2*i];
    nc2=grid->grad[2*i+1];
    if(nc1==-1)
      nc1=nc2;
    if(nc2==-1);
      nc2=nc1;
    if(grid->dzf[i][0]==0)
      dry=1;
    for(k=grid->etop[i];k<grid->Nke[i];k++) {
      C = fabs(phys->u[i][k])*prop->dt/grid->dg[i];
      if(C>CmaxU && !dry && grid->mark[i]!=3 && grid->mark[i]!=2) {
        icu = i;
        kcu = k;
        CmaxU = C;
      }
    }
  }*/

  CmaxW=0;
  for(i=0;i<Nc;i++) 
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      C = 0.5*fabs(phys->w[i][k]+phys->w[i][k+1])*prop->dt/grid->dzz[i][k];
      if(C>CmaxW && (grid->dv[i]+phys->h[i])>=prop->BUFFERHEIGHT && grid->dzz[i][k]>=prop->BUFFERHEIGHT) {
        icw = i;
        kcw = k;
        CmaxW = C;
      }
    }

  progout = (int)(prop->nsteps*(double)prop->ntprog/100);
  if(progout>0 && !(prop->n%progout)) {
    MPI_Reduce(&CmaxU,&allCmaxU,1,MPI_DOUBLE,MPI_MAX,0,comm);
    MPI_Reduce(&CmaxW,&allCmaxW,1,MPI_DOUBLE,MPI_MAX,0,comm);
  }

  if(myproc==0) {
    prop->CmaxU = allCmaxU;
    prop->CmaxW = allCmaxW;
  }

  myalldone=0;
  if(!uflag || !wflag || !sflag || !hflag || CmaxU>prop->Cmax) {
    printf(DASHES);
    printf("Time step %d: Processor %d, Run is blowing up!\n",prop->n,myproc);

    if(CmaxU>prop->Cmax) {
      nc1 = grid->grad[2*icu];
      nc2 = grid->grad[2*icu+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;

      dtsuggestU = CMAXSUGGEST*grid->dg[icu]/fabs(phys->u[icu][kcu]);

      printf("Horizontal Courant number problems:\n");
      printf("  Grid indices: j=%d k=%d (Nke=%d)\n", icu, kcu, grid->Nke[icu]);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xe[icu],grid->ye[icu],
          0.5*(DepthFromDZ(grid,phys,nc1,kcu)+DepthFromDZ(grid,phys,nc2,kcu)));
      printf("  Free-surface heights (on either side): %.3e, %.3e\n",phys->h[nc1],phys->h[nc2]);
      printf("  Layer heights (on either side): \n");
      for(k=0;k<grid->Nkmax;k++) {
	printf("     k=%d, %.3e, %.3e\n",k,grid->dzz[nc1][k],grid->dzz[nc2][k]);
      }
      printf("  Depths (on either side): %.3e, %.3e\n",grid->dv[nc1],grid->dv[nc2]);
      printf("  Bottom slope: %.3e\n",fabs(grid->dv[nc1]-grid->dv[nc2])/grid->dg[icu]);
      printf("  Salinity (on either side): %.3e, %.3e\n",phys->s[nc1][kcu],phys->s[nc2][kcu]);
      printf("  Flux-face height = %.3e, CdB = %.3e\n",grid->dzf[icu][kcu],phys->CdB[icu]);
      printf("  Umax = %.3e\n",phys->u[icu][kcu]);
      printf("  Horizontal grid spacing grid->dg[%d] = %.3e\n",icu,grid->dg[icu]);
      printf("  Horizontal Courant number is CmaxU = %.2f.\n",CmaxU);
      printf("  You specified a maximum of %.2f in suntans.dat\n",prop->Cmax);
      printf("  Your time step size is %.2f.\n",prop->dt);
      printf("  Reducing it to at most %.2f (For C=0.5) might solve this problem.\n",dtsuggestU);
    }

    if(!uflag) {
      printf("Problem with U (U=NaN):\n");
      printf("  Grid indices: j=%d k=%d (Nke=%d)\n", iu, ku, grid->Nke[iu]);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xe[iu],grid->ye[iu],
          0.5*(DepthFromDZ(grid,phys,grid->grad[2*iu],ku)+
            DepthFromDZ(grid,phys,grid->grad[2*iu+1],ku)));	
    }

    if(CmaxW>prop->Cmax && prop->thetaM<0.5 && prop->vertcoord!=2 && prop->vertcoord!=4) {
      dtsuggestW = CMAXSUGGEST*grid->dzz[icw][kcw]/fabs(0.5*(phys->w[icw][kcw]+phys->w[icw][kcw]));

      printf("Vertical Courant number problems:\n");
      printf("  Grid indices: i=%d k=%d (Nkc=%d)\n", icw, kcw, grid->Nkc[icw]);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[icw],grid->yv[icw],DepthFromDZ(grid,phys,icw,kcw));
      printf("  Free-surface height: %.3e\n",phys->h[icw]);
      printf("  Depth: %.3e\n",grid->dv[icw]);
      printf("  Wmax = %.3e (located half-way between faces)\n",0.5*(phys->w[icw][kcw]+phys->w[icw][kcw+1]));
      printf("  Vertical grid spacing dz = %.3e\n",grid->dzz[icw][kcw]);
      printf("  Vertical Courant number is CmaxW = %.2f.\n",CmaxW);
      printf("  You specified a maximum of %.2f in suntans.dat\n",prop->Cmax);
      printf("  Your time step size is %.2f.\n",prop->dt);
      printf("  Reducing it to at most %.2f (For C=0.5) might solve this problem.\n",dtsuggestW);
    }

    if(!wflag) {
      printf("Problem with W (W=NaN):\n");
      printf("  Grid indices: i=%d k=%d (Nkc=%d)\n", iw, kw, grid->Nkc[iw]);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[iw],grid->yv[iw],DepthFromDZ(grid,phys,iw,kw));
    }

    if(!sflag) {
      printf("Problem with the scalar s (s=NaN):\n");
      printf("  Grid indices: i=%d k=%d (Nkc=%d)\n", is, ks, grid->Nkc[is]);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e dv =%.3e h = %.3e\n",grid->xv[is],grid->yv[is],DepthFromDZ(grid,phys,is,ks),grid->dv[is],phys->h[is]);
    }

    if(!hflag) {
      printf("Problem with the free surface (h=NaN):\n");
      printf("  Grid index: i=%d\n", ih);
      printf("  Location: x=%.3e, y=%.3e\n",grid->xv[ih],grid->yv[ih]);
    }
    printf(DASHES);

    myalldone=1;
  }

  MPI_Reduce(&myalldone,&alldone,1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&alldone,1,MPI_INT,0,comm);

  return alldone;
}

/*
 * Function: CheckDZ
 * Usage: CheckDZ(grid,phys,prop,myproc,numprocs,comm);
 * ----------------------------------------------------
 * Check to make sure the vertical grid spacing is >= 0 and that the free surface is
 * not crossing through cells when wetdry != 0.
 *
 */    
int CheckDZ(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm)
{
  int i, k, iz, kz, iw, kw, Nc=grid->Nc, Ne=grid->Ne,nf,ne;
  int zflag=1, wflag=1, myalldone, alldone;
  REAL C, u_im;

  iw=kw=iz=kz=0;

  for(i=0;i<Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      if(grid->dzz[i][k]<=0) {
        zflag=0;
        iz=i;
        kz=k;
        break;
      }
    if(!zflag)
      break;
  }

  if(!prop->wetdry)
    for(i=0;i<Nc;i++) {
      if(grid->ctop[i]!=grid->ctopold[i]) {
        wflag=0;
        iw=i;
        kw=k;
        break;
      }
    }

  myalldone=0;
  if(!zflag || !wflag) {
    printf(DASHES);
    printf("Time step %d: Processor %d, Wetting and drying problems!\n",prop->n,myproc);

    if(!zflag) {
      printf("Problems with the vertical grid spacing:\n");
      printf("  Grid indices: i=%d k=%d\n", iz, kz);
      printf("  Location: x=%.3e, y=%.3e, z=%.3e\n",grid->xv[iz],grid->yv[iz],
          0.5*(DepthFromDZ(grid,phys,grid->grad[2*iz],kz)+
            DepthFromDZ(grid,phys,grid->grad[2*iz+1],kz)));
      printf("  Vertical grid spacing = %.3e <= 0.\n",grid->dzz[iz][kz]);
      printf("Courant number at faces:\n");
      for(nf=0;nf<grid->nfaces[iz];nf++)
      {
        ne=grid->face[iz*grid->maxfaces+nf];
        if(ne!=-1 && kz<grid->Nke[ne] && kz>=grid->etop[ne]){
	  u_im=prop->imfac1*phys->u[ne][kz]
	    +prop->imfac2*phys->u_old[ne][kz]
	    +prop->imfac3*phys->u_old2[ne][kz];
	  C=u_im*prop->dt/grid->dg[ne];
          printf("\tC=%.3e at k=%d edge index=%d cell index=%d\n",C,kz,ne,iz);
        }
      }
    }

    if(!wflag) {
      printf("Cells are wetting and drying although this is not allowed.\n");
      printf("Wetting and drying is only allowed if wetdry=1 in suntans.dat.\n");
      printf("  Grid indices: j=%d k=%d\n", iw, kw);
      printf("  Location: x=%.3e, y=%.3e\n",grid->xv[iw],grid->yv[iw]);
    }
    printf(DASHES);

    myalldone=1;
  }

  MPI_Reduce(&myalldone,&alldone,1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&alldone,1,MPI_INT,0,comm);

  return alldone;
}

/*
 * Function: Progress
 * Usage: Progress(prop,myproc);
 * -----------------------------
 * Output the progress of the calculation to the terminal.
 *
 */
void Progress(propT *prop, int myproc, int numprocs) 
{
  int progout, prog;
  char filename[BUFFERLENGTH];
  FILE *fid;
  REAL timeperstep = (Timer()-t_start)/(prop->n-prop->nstart);
  REAL t_sim, t_rem;

  MPI_GetFile(filename,DATAFILE,"ProgressFile","Progress",myproc);

  if(myproc==0) {
    fid = fopen(filename,"w");
    fprintf(fid,"From %d On %d of %d, t=%.2f (%d%% Complete, %d output)",
        prop->nstart,prop->n,prop->nstart+prop->nsteps,prop->rtime,100*(prop->n-prop->nstart)/prop->nsteps,
        1+(prop->n-prop->nstart)/prop->ntout);      
    fclose(fid);
  }

  if(myproc==0 && prop->ntprog>0 && VERBOSE>0) {
    progout = (int)(prop->nsteps*(double)prop->ntprog/100);
    prog=(int)(100.0*(double)(prop->n-prop->nstart)/(double)prop->nsteps);
    if(progout>0)
      if(!(prop->n%progout)) {
        if(prop->nonlinear) {
          printf("%d%% at %.1es. CmaxU=%.2e, CmaxW=%.2e, %.2e s/step; %.2f s remaining.\n",
              prog,prop->rtime, prop->CmaxU,prop->CmaxW,timeperstep,timeperstep*(prop->nsteps+prop->nstart-prop->n));
        } else {
          printf("%d%% at %.1es. CmaxU=%.2e, %.2e s/step; %.2f s remaining.\n",
              prog,prop->rtime, prop->CmaxU,timeperstep,timeperstep*(prop->nsteps+prop->nstart-prop->n));	  
        }
      }
    if(prop->n==prop->nsteps+prop->nstart) {
      t_sim = Timer()-t_start;
      t_rem = t_sim
        -t_nonhydro-t_predictor-t_source-t_transport-t_turb-t_io-t_check;

      printf("Total simulation time: %.2f s\n",t_sim);
      printf("Average per time step: %.2e s\n",t_sim/prop->nsteps);
      printf("Timing Summary:\n");
      if(prop->nonhydrostatic)
        printf("  Nonhydrostatic pressure: %.2f s (%.2f%%)\n",t_nonhydro,
            100*t_nonhydro/t_sim);
      printf("  Free surface and vertical friction: %.2f s (%.2f%%)\n",t_predictor,
          100*t_predictor/t_sim);
      printf("  Explicit terms: %.2f s (%.2f%%)\n",t_source,
          100*t_source/t_sim);
      if(prop->beta || prop->gamma)
        printf("  Scalar transport: %.2f s (%.2f%%)\n",t_transport,
            100*t_transport/t_sim);
      if(prop->turbmodel)
        printf("  Turbulence: %.2f s (%.2f%%)\n", t_turb,
            100*t_turb/t_sim);
      printf("  Bounds checking: %.2f s (%.2f%%)\n", t_check,
          100*t_check/t_sim);
      printf("  I/O: %.2f s (%.2f%%)\n", t_io,
          100*t_io/t_sim);
      printf("  Remainder: %.2f s (%.2f%%)\n", t_rem,
          100*t_rem/t_sim);
      if(numprocs>1) {
        printf("  Communication time: %.2e s\n",t_comm);
        printf("  Computation/Communication: %.2e\n",(t_sim-t_comm)/t_comm);
      }
    }
  }
}

/*
 * Function: MemoryStats
 * Usage: MemoryStats(myproc,numprocs,comm);
 * -----------------------------------------
 * Print out statistics on total memory and grid points.
 *
 */
void MemoryStats(gridT *grid, int myproc, int numprocs, MPI_Comm comm) {
  int i, ncells, allncells;
  unsigned AllSpace, TotSpacekb = TotSpace>>10;

  ncells=0;
  for(i=0;i<grid->Nc;i++)
    ncells+=grid->Nk[i];

  MPI_Reduce(&TotSpacekb,&(AllSpace),1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&AllSpace,1,MPI_INT,0,comm);
  MPI_Reduce(&ncells,&(allncells),1,MPI_INT,MPI_SUM,0,comm);
  MPI_Bcast(&allncells,1,MPI_INT,0,comm);

  if(numprocs>0)
    printf("Processor %d,  Total memory: %u Mb, %d cells\n",
        myproc,TotSpacekb>>10,ncells);
  if(myproc==0) 
    printf("All processors: %u Mb, %d cells (%d bytes/cell)\n",
        AllSpace>>10,allncells,
        (int)(1024.0*(REAL)AllSpace/(REAL)allncells));
}

/*
 * Function: CheckDivergence
 * Usage: CheckDivergence();
 * -------------------------
 * Checks the divergence and returns an error if divergence is not satisfied.
 *
 */
int CheckDivergence(REAL **dzz, REAL **dzz_old, REAL **u_n_plus_1, REAL **u_n, REAL **u_n_minus_1,
		    REAL **w_im, gridT *grid, propT *prop, REAL tolerance, int myproc) {
  int i, k, iptr, flag, nf, ne, div_local_count, div_da_count;
  REAL div_da, div_da_max, div_local, div_local_max, u_im;

  div_da_max = 0;
  div_local_max = 0;

  div_local_count=0;
  div_da_count=0;  

  for(iptr=grid->celldist[0];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];

    // Only check for divergence in non-boundary cells (flag=0)
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

	  u_im = prop->imfac1*u_n_plus_1[ne][k] + prop->imfac2*u_n[ne][k] + prop->imfac3*u_n_minus_1[ne][k];
	  div_local+=u_im*grid->dzf[ne][k]*grid->normal[i*grid->maxfaces+nf]*grid->df[ne];
	}
	div_local+=grid->Ac[i]*(dzz[i][k]-dzz_old[i][k])/prop->dt;	  	  
	div_local+=grid->Ac[i]*(w_im[i][k]-w_im[i][k+1]);
	div_da+=div_local;
	
	if(k>=grid->ctop[i]) {
	  if(fabs(div_local)>tolerance) {
	    div_local_max = Max(div_local_max,fabs(div_local));
	    div_local_count++;
	  }
	}
      }
      if(fabs(div_da)>tolerance) {
	div_da_max = Max(div_da_max,fabs(div_da));
	div_da_count++;
      }
    }
  }
  
  if(div_local_count>0)
    printf("Proc %d, Max local divergence = %.3e, count exceeding %.3e = %d\n",
	   myproc,div_local_max,tolerance,div_local_count);
  if(div_da_count>0)
    printf("Proc: %d, Max depth-averaged divergence = %.3e, count exceeding %.3e = %d\n",
	   myproc,div_da_max,tolerance,div_da_count);
}


/*
 * Function: CheckGrid
 * Usage: CheckGrid
 * ----------------
 * Check grid quality.
 *
 */
void CheckGrid(gridT *grid, int myproc, MPI_Comm comm) {
  int i, j, k, ne, nf, m, ind, ind1, ind2, count1, count2, count3, normal, nc1, nc2;
  REAL sum, sum_df1, sum_df2, sum_df1_max = 0, sum_df2_max = 0;
  REAL mag, Ax, Ay, Axy;
  REAL tol = 1e-6;

  if(VERBOSE>2 && myproc==0) {
    printf("Checking grid geometry...\n");
  }
  for(i=0;i<grid->Nc;i++) {
    if(grid->Ac[i]==0) printf("Error! Ac[%d]=0!\n",i);

    for(nf=0;nf<grid->nfaces[i];nf++) {
      if(grid->face[i*grid->maxfaces+nf]<0 ||
	 grid->face[i*grid->maxfaces+nf]>=grid->Ne) {
	printf("Error! Problem with face pointer = %d (should satisfy 0<=face<Ne(=%d).\n",
	       grid->face[i*grid->maxfaces+nf],grid->Ne);
      }
      if(grid->neigh[i*grid->maxfaces+nf]<-1 ||
	 grid->neigh[i*grid->maxfaces+nf]>=grid->Nc) {
	printf("Error! Problem with neigh pointer = %d (should satisfy -1<=neigh<Nc(=%d).\n",
	       grid->neigh[i*grid->maxfaces+nf],grid->Nc);
      }
      if(grid->normal[i*grid->maxfaces+nf]!=-1 &&
	 grid->normal[i*grid->maxfaces+nf]!=1) {
	printf("Error! Problem with normal pointer = %d (should satisfy normal=+/-1.\n",
	       grid->normal[i*grid->maxfaces+nf]);
      }
      if(grid->def[i*grid->maxfaces+nf]<0 || IsNan(grid->def[i*grid->maxfaces+nf])) {
	printf("Error! Problem with def = %f.\n",grid->def[i*grid->maxfaces+nf]);
      }
    }
  }

  for(j=0;j<grid->Ne;j++) {
    if(grid->grad[2*j]<-1 || grid->grad[2*j]>=grid->Nc ||
       grid->grad[2*j+1]<-1 || grid->grad[2*j+1]>=grid->Nc ||
       (grid->grad[2*j]==-1 && grid->grad[2*j+1]==-1)) {
      printf("Error! grad1=%d grad2=%d (should satisfy -1<=grad<Nc(=%d).\n",
	     grid->grad[2*j],grid->grad[2*j+1],grid->Nc);
    }
    if(grid->gradf[2*j]<-1 || grid->gradf[2*j+1]>=grid->maxfaces ||
       grid->gradf[2*j+1]<-1 || grid->gradf[2*j+1]>=grid->maxfaces) {
      printf("Error! gradf1=%d, gradf2=%d, should satisfy (-1<=gradf<%d).\n",
	     grid->gradf[2*j],grid->gradf[2*j+1],grid->maxfaces);
    }
    mag=grid->n1[j]*grid->n1[j]+grid->n2[j]*grid->n2[j];
    if(grid->n1[j]>1.0 || grid->n1[j]<-1.0 ||
       grid->n2[j]>1.0 || grid->n2[j]<-1.0 ||
       fabs(mag-1.0)>1e-6) {
      printf("Error! n1=%f, n2=%f, |n1^2+n2^2-1|=%.2e. Should satisfy -1<=n<=+1.\n",
	     grid->n1[j],grid->n2[j],fabs(mag-1.0));
    }
  }

  // Check grid->grad
  /*
  for(j=0;j<grid->Ne;j++) {
    ind=grid->grad[2*j];
    if(ind==-1 || ind>=grid->Nc) {
      printf("Error! grad=%d must satisfy 0<=gradf<%d.\n",
           ind,grid->Nc);
      exit(1);
    }
  }
  */

  /* Check whether sum(n1 df normal) == 0 and sum(n2 df normal) == 0 */
  count1 = count2 = 0;
  for(i=0;i<grid->Nc;i++) {

    sum_df1=0;
    sum_df2=0;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      normal = grid->normal[i*grid->maxfaces+nf];

      sum_df1+=grid->n1[ne]*grid->df[ne]*normal;
      sum_df2+=grid->n2[ne]*grid->df[ne]*normal;
    }
    if(fabs(sum_df1) > sum_df1_max) {
      sum_df1_max = fabs(sum_df1);
    }
    if(fabs(sum_df1) > tol) {
      printf("%f %f 0 0 %e 0\n",grid->xv[i],grid->yv[i],fabs(sum_df1));
      /* Useful for debugging in Matlab
      printf("%f %f 0 0 %e 0\n",grid->xv[i],grid->yv[i],fabs(sum_df1));

      printf("n1 = [");
      for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      printf("%f ",grid->n1[ne]);
      }
      printf("];\n");
      
      printf("n2 = [");
      for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      printf("%f ",grid->n2[ne]);
      }
      printf("];\n");

      printf("xe = [");
      for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      printf("%f ",grid->xe[ne]);
      }
      printf("];\n");

      printf("ye = [");
      for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      printf("%f ",grid->ye[ne]);
      }
      printf("];\n");

      printf("df = [");
      for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];
      printf("%f ",grid->df[ne]);
      }
      printf("];\n");

      printf("normal = [");
      for(nf=0;nf<grid->nfaces[i];nf++) {
      printf("%d ",grid->normal[i*grid->maxfaces+nf]);
      }
      printf("];\n");                  
      printf("\n\n\n");
      */
      
      ind1 = i;
      count1++;
    }
    
    if(fabs(sum_df2) > sum_df2_max) {
      sum_df2_max = fabs(sum_df2);
    }
    if(fabs(sum_df2) > tol) {
      /* Useful for debugging in Matlab      
      printf("0 0 %f %f 0 %e\n",grid->xv[i],grid->yv[i],fabs(sum_df2));
      */
      ind2 = i;
      count2++;
    }    
  }

  if(VERBOSE > 1) {
    if(count1>0) {
      printf("Warning: Grid error on proc %d: \n",myproc);
      printf("\tmax(Sum n1 df) = %.3e at (x,y)=(%.2e,%.2e) count>%.1e: %d\n",
	     sum_df1_max,grid->xv[ind1],grid->yv[ind1],tol,count1);
    }
    if(count2>0) {
      printf("Warning: Grid error on proc %d: \n",myproc);
      printf("\tmax(Sum n2 df) = %.3e at (x,y)=(%.2e,%.2e) count>%.1e: %d\n",
	     sum_df2_max,grid->xv[ind2],grid->yv[ind2],tol,count2);
    }    
  }

  /* Check whether Sum(n1^2 def df) = Ac and Sum(n1^2 N df) = Ac, which
     is needed for Perot's method. */
  count1=count2=count3=0;
  for(i=0;i<grid->Nc;i++) {

    Ax=Ay=Axy=0;
    for(nf=0;nf<grid->nfaces[i];nf++) {
      ne = grid->face[i*grid->maxfaces+nf];

      Ax+=pow(grid->n1[ne],2.0)*grid->def[i*grid->maxfaces+nf]*grid->df[ne];
      Ay+=pow(grid->n2[ne],2.0)*grid->def[i*grid->maxfaces+nf]*grid->df[ne];
      Axy+=grid->n1[ne]*grid->n2[ne]*grid->def[i*grid->maxfaces+nf]*grid->df[ne];            
    }
    if(fabs(Ax-grid->Ac[i])/grid->Ac[i]>tol) {
      printf("%f %f 0 0 0 \n",grid->xv[i],grid->yv[i]);
      printf("Ax = %f, Ac = %f\n",Ax,grid->Ac[i]);
      count1++;
    }
    if(fabs(Ay-grid->Ac[i])/grid->Ac[i]>tol) {
      //      printf("Ay = %f, Ac = %f\n",Ay,grid->Ac[i]);
      count2++;     
    }
    if(fabs(Axy)/grid->Ac[i]>tol) {
      //      printf("Axy = %.2e\n",Axy);
      count3++;
    }    
  }

  if(VERBOSE > 1) {
    if(count1>0 || count2>0 || count3>0) {
      printf("Warning: Grid error on proc %d: \n",myproc);
    }
    if(count1>0) {
      printf("\t Sum(n1^2 def df) is not equal to Ac - count>%.1e: %d\n",tol,count1);
    }
    if(count2>0) {
      printf("\t Sum(n2^2 def df) is not equal to Ac - count>%.1e: %d\n",tol,count2);      
    }
    if(count3>0) {
      printf("\t Sum(n1 n2 def df) is not equal to 0 - count>%.1e: %d\n",tol,count3);      
    }
  }
}
