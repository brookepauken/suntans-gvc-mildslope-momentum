/*
 * Boundaries test file.
 *
 */
#include "boundaries.h"
#include "sediments.h"
#include "wave.h"
#include "subgrid.h"
#include "vertcoordinate.h"
static void SetUVWH(gridT *grid, physT *phys, propT *prop, int ib, int j, int boundary_index, REAL boundary_flag);

/*
 * Function: OpenBoundaryFluxes
 * Usage: OpenBoundaryFluxes(q,ubnew,ubn,grid,phys,prop);
 * ----------------------------------------------------
 * This will update the boundary flux at the edgedist[2] to edgedist[3] edges.
 * 
 * Note that phys->uold,vold contain the velocity at time step n-1 and 
 * phys->uc,vc contain it at time step n.
 *
 * The radiative open boundary condition does not work yet!!!  For this reason c[k] is
 * set to 0
 *
 */
void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop) {
  int j, jptr, ib, k, forced;
  REAL *uboundary = phys->a, **u = phys->uc, **v = phys->vc, **uold = phys->uold, **vold = phys->vold;
  REAL z, c0, c1, C0, C1, dt=prop->dt, u0, u0new, uc0, vc0, uc0old, vc0old, ub0;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];

    for(k=grid->etop[j];k<grid->Nke[j];k++) 
      ub[j][k]=phys->boundary_u[jptr-grid->edgedist[2]][k]*grid->n1[j]+
	phys->boundary_v[jptr-grid->edgedist[2]][k]*grid->n2[j];
  }
}

/*
 * Function: BoundaryScalars
 * Usage: BoundaryScalars(boundary_s,boundary_T,grid,phys,prop);
 * -------------------------------------------------------------
 * This will set the values of the scalars at the open boundaries.
 * 
 */
void BoundaryScalars(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int jptr, j, iptr, i, ib, k;
  REAL z;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      z=phys->h[ib];
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
	z-=grid->dzz[ib][k]/2.0;
	if(z<-7.5 && z>-12.5)
	  phys->boundary_T[jptr-grid->edgedist[2]][k]=1;
	else
	  phys->boundary_T[jptr-grid->edgedist[2]][k]=0;	  
	phys->boundary_s[jptr-grid->edgedist[2]][k]=phys->s[ib][k];
	z-=grid->dzz[ib][k]/2.0;	
      }
  }

  // At the ocean boundary
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    
    for(k=0;k<grid->ctop[i];k++) {
      phys->s[i][k]=0;
      phys->T[i][k]=0;
    } 
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      phys->s[i][k]=0;
      phys->T[i][k]=0;
    }
  }
}

/*
 * Function: BoundaryVelocities
 * Usage: BoundaryVelocities(grid,phys,prop,myproc);
 * -------------------------------------------------
 * This will set the values of u,v,w, and h at the boundaries.
 * 
 */
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int jptr, j, ib, k,iptr,i,boundary_index;
  REAL z;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      phys->boundary_u[jptr-grid->edgedist[2]][k]=0;
      phys->boundary_v[jptr-grid->edgedist[2]][k]=0;
      phys->boundary_w[jptr-grid->edgedist[2]][k]=0;
    }
  }

  // At the ocean boundary
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    phys->h[i]=0;
  }


}

static void SetUVWH(gridT *grid, physT *phys, propT *prop, int ib, int j, int boundary_index, REAL boundary_flag) {
  int k;

  if(boundary_flag==open) {
    phys->boundary_h[boundary_index]=phys->h[ib];
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      phys->boundary_u[boundary_index][k]=phys->uc[ib][k];
      phys->boundary_v[boundary_index][k]=phys->vc[ib][k];
      phys->boundary_w[boundary_index][k]=0.5*(phys->w[ib][k]+phys->w[ib][k+1]);
    }
  } else {
    phys->boundary_h[boundary_index]=prop->amp*fabs(cos(prop->omega*prop->rtime));
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      phys->boundary_u[boundary_index][k]=phys->u[j][k]*grid->n1[j];
      phys->boundary_v[boundary_index][k]=phys->u[j][k]*grid->n2[j];
      phys->boundary_w[boundary_index][k]=0.5*(phys->w[ib][k]+phys->w[ib][k+1]);
    }
  }
}
	
/*
 * Function: WindStress
 * Usage: WindStress(grid,phys,prop,myproc);
 * -----------------------------------------
 * Set the wind stress.
 *
 */
void WindStress(gridT *grid, physT *phys, propT *prop, metT *met, int myproc) {
  int j, jptr;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];
    
    //phys->tau_T[j]=grid->n2[j]*prop->tau_T;
    phys->tau_T[j]=0;
    //phys->tau_B[j]=0;
    //printf("c is this: %f\n", prop->cwave);
    phys->tau_B[j]=grid->n1[j]*prop->tau_T;
    //phys->tau_B[j]=grid->n1[j]*cos(PI)*prop->cwave*prop->cwave;
  }
}

//added part
/*
 * Function: BoundarySediment
 * Usage: BoundarySediment(boundary_s,boundary_T,grid,phys,prop);
 * -------------------------------------------------------------
 * This will set the values of the suspended sediment concentration
 * at the open boundaries.
 * 
 */
void BoundarySediment(gridT *grid, physT *phys, propT *prop) {
  int jptr, j, ib, k,nosize,i,iptr;
  REAL z;

  // At the upstream boundary
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    j=grid->edgep[jptr];
    ib=grid->grad[2*j];
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
        sediments->boundary_sediC[nosize][jptr-grid->edgedist[2]][k]=0;
      }
    }
  }

  // At the ocean boundary
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      for(k=0;k<grid->ctop[i];k++) {
        sediments->SediC[nosize][i][k]=0;
      } 
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        sediments->SediC[nosize][i][k]=0;
      }
    }
  }
}
/*
 * Function: WindSpeedandDirection
 * usage: calculate Wind field when Wind is not constant
 * -------------------------------------------------
 * calculate Uwind and Winddir
 *
 */
void FetchWindSpeedandDirection(gridT *grid, propT *prop, int myproc){
   int i;
   for(i=0;i<grid->Nc;i++){
     wave->Uwind[i]=0;
     wave->Winddir[i]=0;
   }
}
void InitBoundaryData(propT *prop, gridT *grid, int myproc, MPI_Comm comm){}
void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc, MPI_Comm comm){}
void UserDefinedFunction(gridT *grid, physT *phys, propT *prop,int myproc)
{
  /*
  static FILE *cfid = NULL;
  static REAL *ds = NULL;
  if(prop->vertcoord==4 && prop->ntout>0 &&
     (!(prop->n % prop->ntout) || prop->n==1+prop->nstart)) {
    REAL uim, sum, sums, alpha, R, Rmax=0.0;
    REAL l[5]={0,0,0,0,0};
    int i, k, iptr, nf2, ne2;

    if(cfid==NULL) { char fn[BUFFERLENGTH];
      sprintf(fn,"%s/consist.dat.%d",DATADIR,myproc); cfid=fopen(fn,"w"); }
    if(ds==NULL) ds = (REAL *)malloc(grid->Nkmax*sizeof(REAL));

    for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++) {
      i = grid->cellp[iptr];

      // --- R: free-surface residual with the corrected velocity ---
      sum = 0.0;
      for(nf2=0; nf2<grid->nfaces[i]; nf2++) {
        ne2 = grid->face[i*grid->maxfaces+nf2];
        for(k=grid->etop[ne2]; k<grid->Nke[ne2]; k++) {
          uim = prop->imfac1*phys->u[ne2][k]
              + prop->imfac2*phys->u_old[ne2][k]
	    + prop->imfac3*phys->u_old2[ne2][k];
          sum += uim*grid->dzf[ne2][k]*grid->df[ne2]
	    *grid->normal[i*grid->maxfaces+nf2];
        }
      }
      R = (phys->h[i]-phys->h_old[i])/prop->dt + sum/grid->Ac[i];
      l[0] += grid->Ac[i]*R*R;
      l[1] += grid->Ac[i];
      if(fabs(R)>Rmax) Rmax = fabs(R);

      // --- dzz*: transport with corrected u, then floor, then Step 3 --- 
      for(k=grid->ctop[i]; k<grid->Nk[i]; k++) {
        ds[k] = grid->dzzold[i][k];
        for(nf2=0; nf2<grid->nfaces[i]; nf2++) {
          ne2 = grid->face[i*grid->maxfaces+nf2];
          uim = prop->imfac1*phys->u[ne2][k]
              + prop->imfac2*phys->u_old[ne2][k]
	    + prop->imfac3*phys->u_old2[ne2][k];
          ds[k] -= prop->dt*uim*grid->dzf[ne2][k]
	    * grid->normal[i*grid->maxfaces+nf2]*grid->df[ne2]/grid->Ac[i];
        }
        if(ds[k] < vert->vertdzmin) ds[k] = vert->vertdzmin;
      }
      sums = 0.0;
      for(k=grid->ctop[i]; k<grid->Nk[i]; k++) sums += ds[k];
      alpha = (phys->h[i]+grid->dv[i]-phys->zB[i])/sums;
      for(k=grid->ctop[i]; k<grid->Nk[i]; k++) {
        ds[k] *= alpha;
        l[2] += grid->Ac[i]*(ds[k]-grid->dzz[i][k])*(ds[k]-grid->dzz[i][k]);
        l[3] += grid->Ac[i]*grid->dzz[i][k]*grid->dzz[i][k];
      }
    }
    l[4] = Rmax;
    fprintf(cfid,"%d %.16e %.16e %.16e %.16e %.16e\n",
            prop->n, l[0], l[1], l[2], l[3], l[4]);
    fflush(cfid);
  }
  */
}
