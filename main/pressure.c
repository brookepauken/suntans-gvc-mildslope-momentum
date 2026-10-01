/*
 * File: pressure.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains functions needed to solve for the nonhydrostatic pressure.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "pressure.h"
#include "sendrecv.h"
#include "memory.h"
#include "vertcoordinate.h"
#include "subgrid.h"
#include "rand.h"

/*
 * Private Function declarations.
 *
 */
static void ConditionQ(REAL **x, gridT *grid, physT *phys, propT *prop, int myproc, 
		       MPI_Comm comm);
static void Preconditioner(REAL **x, REAL **xc, REAL **coef, gridT *grid, physT *phys, 
			   propT *prop);
static void GuessQ(REAL **q, REAL **wold, REAL **w, gridT *grid, physT *phys, 
		   propT *prop, int myproc, int numprocs, MPI_Comm comm);
static void OperatorQC(REAL **coef, REAL **fcoef, REAL **x, REAL **y, REAL **c, 
		       gridT *grid, physT *phys, propT *prop);
static void QCoefficients(REAL **coef, REAL **fcoef, REAL **c, gridT *grid, 
			  physT *phys, propT *prop);
static void OperatorQ(REAL **coef, REAL **x, REAL **y, REAL **c, gridT *grid, 
		      physT *phys, propT *prop);
static void OperatorQSlope(REAL **x, REAL **y, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void Compute_dqdxi_Face(REAL **dqdxi_face, REAL **x, gridT *grid, physT *phys, int myproc, MPI_Comm comm);
static REAL InnerProduct3(REAL **x, REAL **y, gridT *grid, int myproc, int numprocs, 
			  MPI_Comm comm);

/*************************************************************************/
/*                                                                       */
/* Public functions                                                      */
/*                                                                       */
/*************************************************************************/

void BiCGSolveQ(REAL **q, REAL **src, REAL **c, gridT *grid, physT *phys,
		propT *prop, int myproc, int numprocs, MPI_Comm comm) {
  int i, k, iptr, n, niters;
  REAL **x, **u, **y, **r, **p, **z, **r0, **v, **t, **s, eps, eps0;
  REAL rho0, rho, alpha, omega, beta, pmax, pmin, mypmax,mypmin;

  int use_operator_slope = 1; // Use the OperatorQSlope function. For debugging

  u = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "BiCGSolveQ");
  y = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "BiCGSolveQ");    
  z = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "BiCGSolveQ");
  v = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "BiCGSolveQ");  
  s = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "BiCGSolveQ");  
  t = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "BiCGSolveQ");  
  r0 = (REAL **)SunMalloc(grid->Nc*sizeof(REAL *), "BiCGSolveQ");  
  for(i=0;i<grid->Nc;i++) {
    u[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"BiCGSolveQ");    
    y[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"BiCGSolveQ");    
    z[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"BiCGSolveQ");
    v[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"BiCGSolveQ");    
    s[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"BiCGSolveQ");
    t[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"BiCGSolveQ");
    r0[i] = (REAL *)SunMalloc(grid->Nk[i]*sizeof(REAL),"BiCGSolveQ");
  }
  x = q;
  r = phys->stmp3;
  p = src;

  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      r[i][k]=x[i][k]=0;
  }

  niters = prop->qmaxiters;

  // Compute min/max of source to create the random array r0
  /*
  mypmax=-INFTY;
  mypmin=INFTY;  
  for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      if(p[i][k]<=mypmin)
	mypmin=p[i][k];
      if(p[i][k]>=mypmax)
	mypmax=p[i][k];
    }
  }
  MPI_Reduce(&mypmax,&(pmax),1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Bcast(&pmax,1,MPI_DOUBLE,0,comm);
  MPI_Reduce(&mypmin,&(pmin),1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Bcast(&pmin,1,MPI_DOUBLE,0,comm);    
  */
  
  for(i=0;i<grid->Nc;i++) {
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      z[i][k] = 0;
      v[i][k] = 0;
      s[i][k] = 0;
      t[i][k] = 0;
      r0[i][k] = 0;
    }
  }

  // Create the coefficients for the preconditioner and
  // OperatorQ()
  QCoefficients(phys->wtmp,phys->qtmp,c,grid,phys,prop);  
      
  ISendRecvCellData3D(x, grid, myproc, comm);
  if(use_operator_slope)
    OperatorQSlope(x,z,grid,phys,prop,myproc,comm);
  else 
    OperatorQ(phys->wtmp,x,z,c,grid,phys,prop);

  for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      p[i][k] = p[i][k]-z[i][k]; 
      r[i][k] = p[i][k];  
      r0[i][k] = r[i][k];
      //      r0[i][k]=frand(pmin,pmax);
    }
  }
   
  for(iptr=grid->celldist[1]; iptr<grid->celldist[2]; iptr++){
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {    
      p[i][k] = 0;
    }
  }
  eps0=eps=sqrt(InnerProduct3(r,r,grid,myproc,numprocs,comm));

  if (eps > prop->qepsilon){
    if(!prop->resnorm) eps0 = 1;

    // Main BiCGStab Iteration
    for(n = 0; n < niters && eps!=0; n++){ 
      rho0 = InnerProduct3(r0,r,grid,myproc,numprocs,comm);

      // Preconditioner: y = inv(K)*p
      // v = Ay instead of v = Ap
      if(use_operator_slope) {
	if(prop->qprecond==2) {
	  Preconditioner(p,y,phys->wtmp,grid,phys,prop);

	  ISendRecvCellData3D(y, grid, myproc, comm);	  
	  OperatorQSlope(y,v,grid,phys,prop,myproc,comm);      
	} else {
	  ISendRecvCellData3D(p, grid, myproc, comm);
	  OperatorQSlope(p,v,grid,phys,prop,myproc,comm);
	}
      } else {
	ISendRecvCellData3D(p, grid, myproc, comm);		
	OperatorQ(phys->wtmp,p,v,c,grid,phys,prop);
      }
      alpha =  rho0/InnerProduct3(r0,v,grid,myproc,numprocs,comm);

      for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
      	i = grid->cellp[iptr];

	for(k=grid->ctop[i];k<grid->Nk[i];k++) {    		
	  s[i][k] = r[i][k]-alpha*v[i][k];
	}
      }

      eps = sqrt(InnerProduct3(s, s, grid, myproc, numprocs, comm));

      if(VERBOSE>3 && myproc==0) printf("BiCGSolve Iteration (1): %d, resid=%e.\n",n,eps);
      if(eps<prop->qepsilon){
	break;
      }

      // Preconditioner: z = inv(K)*s
      // t = Az instead of t = As
      // Preconditioner: u = inv(K)*t
      if(use_operator_slope) {
	if(prop->qprecond==2) {
	  Preconditioner(s,z,phys->wtmp,grid,phys,prop);

	  ISendRecvCellData3D(z, grid, myproc, comm);	  	  
	  OperatorQSlope(z,t,grid,phys,prop,myproc,comm);
	  Preconditioner(t,u,phys->wtmp,grid,phys,prop);      
	} else {
	  ISendRecvCellData3D(s, grid, myproc, comm);
	  OperatorQSlope(s,t,grid,phys,prop,myproc,comm);
	}
      } else {
	ISendRecvCellData3D(s, grid, myproc, comm);	
	OperatorQ(phys->wtmp,s,t,c,grid,phys,prop);
      }

      // Preconditioner: omega = ip(inv(K)*t, inv(K)*s)/ip(inv(K)*t,inv(K)*t);
      if(prop->qprecond==2) {
	omega = InnerProduct3(u, z, grid, myproc, numprocs, comm)/
	  InnerProduct3(u, u, grid, myproc, numprocs, comm);
      } else {
	omega = InnerProduct3(t, s, grid, myproc, numprocs, comm)/
	  InnerProduct3(t, t, grid, myproc, numprocs, comm);
      }

      for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
      	i = grid->cellp[iptr];

	if(prop->qprecond==2) {
	  for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	    // Preconditioner: x += alpha*y + omega*z
	    x[i][k] += alpha*y[i][k] + omega*z[i][k];
	  }
	} else {
	  for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	    x[i][k] += alpha*p[i][k] + omega*s[i][k];
	  }
	}
	
	for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	  r[i][k] = s[i][k] - omega*t[i][k];
	}
      }

      eps = sqrt(InnerProduct3(r, r, grid, myproc, numprocs, comm));      
      if(VERBOSE>3 && myproc==0) printf("BiCGSolve Iteration (2): %d, resid=%e\n",n,eps);
      if(eps < prop->qepsilon) {
	break;
      }
      
      rho = InnerProduct3(r, r0, grid, myproc, numprocs, comm);      
      beta = alpha/omega*rho/rho0;

      for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
	i = grid->cellp[iptr];

	for(k=grid->ctop[i];k<grid->Nk[i];k++) {    	
	  p[i][k] = r[i][k]+beta*(p[i][k]-omega*v[i][k]); 
	}
      }

      if(sqrt(fabs(rho))<1e-6) {
	for(iptr=grid->celldist[0]; iptr<grid->celldist[1]; iptr++){
	  i = grid->cellp[iptr];
	  
	  for(k=grid->ctop[i];k<grid->Nk[i];k++) {
	    r0[i][k] = r[i][k];
	    p[i][k] = r[i][k];
	  }
	}
      }
      //      if(beta!=beta) exit(1);
    }
  }

  if(myproc==0 && VERBOSE>2) {
    if(n==niters) {
      printf("Warning... Step %d, pressure iteration not converging after %d iterations! RES=%e > %.2e\n",
	     prop->n, n, eps, prop->qepsilon);
    } else {
      printf("Step %d, BiCGSolve pressure converged after %d iterations, rsdl = %e < %e\n",
	     prop->n, n, eps, prop->qepsilon);
    }
  }
  ISendRecvCellData3D(x, grid, myproc, comm);

  for(i=0;i<grid->Nc;i++) {
    SunFree(u[i], grid->Nk[i]*sizeof(REAL), "BiCGSolveQ");
    SunFree(y[i], grid->Nk[i]*sizeof(REAL), "BiCGSolveQ");    
    SunFree(z[i], grid->Nk[i]*sizeof(REAL), "BiCGSolveQ");
    SunFree(v[i], grid->Nk[i]*sizeof(REAL), "BiCGSolveQ");
    SunFree(s[i], grid->Nk[i]*sizeof(REAL), "BiCGSolveQ");
    SunFree(t[i], grid->Nk[i]*sizeof(REAL), "BiCGSolveQ");
    SunFree(r0[i], grid->Nk[i]*sizeof(REAL), "BiCGSolveQ");
  }
  SunFree(u, grid->Nc*sizeof(REAL *), "BiCGSolveQ");
  SunFree(y, grid->Nc*sizeof(REAL *), "BiCGSolveQ");  
  SunFree(z, grid->Nc*sizeof(REAL *), "BiCGSolveQ");
  SunFree(v, grid->Nc*sizeof(REAL *), "BiCGSolveQ");
  SunFree(s, grid->Nc*sizeof(REAL *), "BiCGSolveQ");
  SunFree(t, grid->Nc*sizeof(REAL *), "BiCGSolveQ");
  SunFree(r0, grid->Nc*sizeof(REAL *), "BiCGSolveQ");  
}

/*
 * Function: CGSolveQ
 * Usage: CGSolveQ(q,src,c,grid,phys,prop,myproc,numprocs,comm);
 * -------------------------------------------------------------
 * Solve for the nonhydrostatic pressure with the preconditioned
 * conjugate gradient algorithm.
 *
 * The preconditioner stores the diagonal preconditioning elements in 
 * the temporary c array.
 *
 * This function replaces q with x and src with p.  phys->uc and phys->vc
 * are used as temporary arrays as well to store z and r.
 *
 */
void CGSolveQ(REAL **q, REAL **src, REAL **c, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, MPI_Comm comm) {

  int i, iptr, k, n, niters;

  REAL **x, **r, **rtmp, **p, **z, mu, nu, alpha, alpha0, eps, eps0;

  z = phys->stmp2;
  x = q;
  r = phys->stmp3;
  rtmp = phys->stmp4;
  p = src;

  // Compute the preconditioner and the preconditioned solution
  // and send it to neighboring processors
  if(prop->qprecond==1) {
    ConditionQ(c,grid,phys,prop,myproc,comm);
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        p[i][k]/=c[i][k];
        x[i][k]*=c[i][k];
      }
    }
  }
  ISendRecvCellData3D(x,grid,myproc,comm);

  niters = prop->qmaxiters;

  // Create the coefficients for the operator
  QCoefficients(phys->wtmp,phys->qtmp,c,grid,phys,prop);

  // Initialization for CG
  if(prop->qprecond==1) OperatorQC(phys->wtmp,phys->qtmp,x,z,c,grid,phys,prop);
  else OperatorQ(phys->wtmp,x,z,c,grid,phys,prop);
  //  else OperatorQSlope(x,z,grid,phys,prop,myproc,comm);  

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      r[i][k] = p[i][k]-z[i][k];
  }

  if(prop->qprecond==2) {
    Preconditioner(r,rtmp,phys->wtmp,grid,phys,prop);
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
        p[i][k] = rtmp[i][k];
    }
    alpha = alpha0 = InnerProduct3(r,rtmp,grid,myproc,numprocs,comm);
  } else {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
        p[i][k] = r[i][k];
    }
    alpha = alpha0 = InnerProduct3(r,r,grid,myproc,numprocs,comm);
  }
  if(!prop->resnorm) alpha0 = 1;

  if(prop->qprecond==2)
    eps=eps0=InnerProduct3(r,r,grid,myproc,numprocs,comm);
  else
    eps=eps0=alpha0;

  // Iterate until residual is less than prop->qepsilon
  for(n=0;n<niters && eps!=0;n++) {

    ISendRecvCellData3D(p,grid,myproc,comm);
    if(prop->qprecond==1) OperatorQC(phys->wtmp,phys->qtmp,p,z,c,grid,phys,prop);
    else OperatorQ(phys->wtmp,p,z,c,grid,phys,prop);
    //    else OperatorQSlope(p,z,grid,phys,prop,myproc,comm);
    
    mu = 1/alpha;
    nu = alpha/InnerProduct3(p,z,grid,myproc,numprocs,comm);

    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        x[i][k] += nu*p[i][k];
        r[i][k] -= nu*z[i][k];
      }
    }
    if(prop->qprecond==2) {
      Preconditioner(r,rtmp,phys->wtmp,grid,phys,prop);
      alpha = InnerProduct3(r,rtmp,grid,myproc,numprocs,comm);
      mu*=alpha;
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];

        for(k=grid->ctop[i];k<grid->Nk[i];k++) 
          p[i][k] = rtmp[i][k] + mu*p[i][k];
      }
    } else {
      alpha = InnerProduct3(r,r,grid,myproc,numprocs,comm);
      mu*=alpha;
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];

        for(k=grid->ctop[i];k<grid->Nk[i];k++) 
          p[i][k] = r[i][k] + mu*p[i][k];
      }
    }

    if(prop->qprecond==2)
      eps=InnerProduct3(r,r,grid,myproc,numprocs,comm);
    else
      eps=alpha;

    //    if(myproc==0) printf("%d %e %e %e %e\n",myproc,mu,nu,eps0,sqrt(eps/eps0));
    if(VERBOSE>2 && myproc==0) printf("CGSolve Pressure Iteration: %d, resid=%e\n",n,sqrt(eps/eps0));
    if(sqrt(eps/eps0)<prop->qepsilon || eps0==0 || eps0!=eps0) 
      break;
  }

  if(eps0==0 || eps0!=eps0) {
    if(myproc==0) printf("Error..Time step %d, norm of pressure source is %f\n",prop->n,eps0);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  
  if(myproc==0 && VERBOSE>2) {
    if(n==niters)  printf("Warning... Time step %d, Pressure iteration not converging after %d steps! RES=%e > %.2e\n",
			  prop->n,n,sqrt(eps/eps0),prop->qepsilon);
    else printf("Time step %d, CGSolve pressure converged after %d iterations, res=%e < %.2e\n",
		prop->n,n,sqrt(eps/eps0),prop->qepsilon);
  }

  // Rescale the preconditioned solution 
  if(prop->qprecond==1) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];

      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
        x[i][k]/=c[i][k];
    }
  }

  // Send the solution to the neighboring processors
  ISendRecvCellData3D(x,grid,myproc,comm);
}

/*
 * Function: Corrector
 * Usage: Corrector(qc,grid,phys,prop,myproc,numprocs,comm);
 * ---------------------------------------------------------
 * Correct the horizontal velocity field with the pressure correction.
 * Do not correct velocities for which D[j]==0 since these are boundary
 * cells.
 *
 * After correcting the horizontal velocity, update the total nonhydrostatic
 * pressure with the pressure correction.
 *
 */
void Corrector(REAL **qc, gridT *grid, physT *phys, propT *prop, int myproc,
	       int numprocs, MPI_Comm comm) {

  int i, iptr, j, jptr, k, nc1, nc2;
  REAL **dqdxi_face = vert->uf;

  if(prop->include_slope_terms) {
    // Compute the cell-centered gradients dqdx and dqdy    
    // the first integer is the direction (0=x,1=y), the second zero means this is cell-centered (1=w-centered)
    ComputeCellAveragedHorizontalGradientCell(vert->dqdx,0,qc,0,grid,prop,phys,myproc);
    ComputeCellAveragedHorizontalGradientCell(vert->dqdy,1,qc,0,grid,prop,phys,myproc);

    Compute_dqdxi_Face(dqdxi_face,qc,grid,phys,myproc,comm);
  }
  
  // Correct the horizontal velocity only if this is not a boundary point.
  // no boundary points are corrected for non-hydrostatic pressure!!!
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];
    
    if(phys->D[j]!=0 && grid->etop[j]<grid->Nke[j]-1) {
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
        phys->u[j][k]-=prop->dt*(qc[nc1][k]-qc[nc2][k])/grid->dg[j];

	// This is the dt/J dz/dn dq/dzeta term
	if(prop->include_slope_terms)
	  phys->u[j][k]+=prop->dt*dqdxi_face[j][k]*(vert->zc[nc1][k]-vert->zc[nc2][k])/grid->dg[j];
      }
    }
  }

  // Correct the vertical velocity
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr]; 
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
      phys->w[i][k]-=2.0*prop->dt/(grid->dzz[i][k-1]+grid->dzz[i][k])*
        (qc[i][k-1]-qc[i][k]);
    phys->w[i][grid->ctop[i]]+=2.0*prop->dt/grid->dzz[i][grid->ctop[i]]*
      qc[i][grid->ctop[i]];
  }

  // Update the pressure depending on value of pressureMethod:
  // pressureMethod == 0 : Projection
  // pressureMethod == 1 : Correction (default)
  if(prop->pressureMethod==1) {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      if(grid->ctop[i]<grid->Nk[i]-1)
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  phys->q[i][k]+=qc[i][k];
    }
  } else {
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      if(grid->ctop[i]<grid->Nk[i]-1)
	for(k=grid->ctop[i];k<grid->Nk[i];k++)
	  phys->q[i][k]=qc[i][k];
    }
  }
}

/*
 * Function: ComputeQSource
 * Usage: ComputeQSource(src,grid,phys,prop,myproc,numprocs);
 * ----------------------------------------------------------
 * Compute the source term for the nonhydrostatic pressure by computing
 * the divergence of the predicted velocity field, which is in phys->u and
 * phys->w.  The upwind flux face heights are used in order to ensure
 * consistency with continuity.
 *
 */
void ComputeQSource(REAL **src, gridT *grid, physT *phys, propT *prop, int myproc, int numprocs) {

  int i, iptr, j, jptr, k, nf, ne, nc1, nc2;
  REAL *ap=phys->a, *am=phys->b, thetafactor=(1-prop->theta)/prop->theta,fac1,fac2,fac3,flux;

  // new implicit method
  fac1=prop->imfac1;
  fac2=prop->imfac2;
  fac3=prop->imfac3;

  // Should enforce continuity at step n+1
  fac1=1;
  fac2=0;
  fac3=0;
  
  // for each cell
  for(i=0;i<grid->Nc;i++){
    // initialize to zero
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      src[i][k] = 0;
  }

  // for each computational cell
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    // get the cell pointer
    i = grid->cellp[iptr];
    /* VERTICAL CONTRIBUTION */
    // over all cells that are defined to a depth
    if(prop->vertcoord==1)
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        // compute the vertical contributions to the source term
        if(!prop->subgrid)
          src[i][k] = grid->Ac[i]*(phys->w[i][k]-phys->w[i][k+1])
            +fac2/fac1*grid->Ac[i]*(phys->w_old[i][k]-phys->w_old[i][k+1])+fac3/fac1*grid->Ac[i]*(phys->w_old2[i][k]-phys->w_old2[i][k+1]);
        else
          src[i][k] = subgrid->Acveff[i][k]*phys->w[i][k]-subgrid->Acveff[i][k+1]*phys->w[i][k+1]
            +fac2/fac1*(subgrid->Acveffold[i][k]*phys->w_old[i][k]-subgrid->Acveffold[i][k+1]*phys->w_old[i][k+1])
            +fac3/fac1*(subgrid->Acveffold2[i][k]*phys->w_old2[i][k]-subgrid->Acveffold2[i][k+1]*phys->w_old2[i][k+1]);
      }
    else
      // modification for the new generalized vertical coordinate
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        // compute the vertical contributions to the source term
        if(!prop->subgrid)
          src[i][k] = grid->Ac[i]*(vert->U3[i][k]-vert->U3[i][k+1])+
                fac2/fac1*grid->Ac[i]*(vert->U3_old[i][k]-vert->U3_old[i][k+1])+
                fac3/fac1*grid->Ac[i]*(vert->U3_old2[i][k]-vert->U3_old2[i][k+1]);
        else
          src[i][k] =subgrid->Acveff[i][k]*vert->U3[i][k]-subgrid->Acveff[i][k+1]*vert->U3[i][k+1]+
             fac2/fac1*(subgrid->Acveffold[i][k]*vert->U3_old[i][k]-subgrid->Acveffold[i][k+1]*vert->U3_old[i][k+1])+
             fac3/fac1*(subgrid->Acveffold2[i][k]*vert->U3_old2[i][k]-subgrid->Acveffold2[i][k+1]*vert->U3_old2[i][k+1]);
      }        

    /* HORIZONTAL CONTRIBUTION */
    // over each face to get the horizontal contributions to the source term
    // no change is needed for the new generalized vertical coordinate
    for(nf=0;nf<grid->nfaces[i];nf++) {

      // get the edge pointer
      ne = grid->face[i*grid->maxfaces+nf];
      // for each of the defined edges over depth
      for(k=grid->ctop[i];k<grid->Nke[ne];k++) 
        // compute the horizontal source term via the (D_H)(u^*)
        src[i][k]+=(phys->u[ne][k]+fac2/fac1*phys->u_old[ne][k]+fac3/fac1*
          phys->u_old2[ne][k])*grid->dzf[ne][k]*grid->normal[i*grid->maxfaces+nf]*grid->df[ne];
    }
   
 
    // divide final result by dt
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      src[i][k]/=prop->dt;
  }

  // D[j] is used in OperatorQ, and it must be zero to ensure no gradient
  // at the hydrostatic faces.
  // this can artificially be used to control the boundary condition
  for(j=0;j<grid->Ne;j++) {
    phys->D[j]=grid->df[j]/grid->dg[j];
  }
} 

/*************************************************************************/
/*                                                                       */
/* Private functions                                                     */
/*                                                                       */
/*************************************************************************/

/*
 * Function: ConditionQ
 * Usage: ConditionQ(x,grid,phys,prop,myproc,comm);
 * ------------------------------------------------
 * Compute the magnitude of the diagonal elements of the coefficient matrix
 * for the pressure-Poisson equation and place it into x after taking its square root.
 *
 */
static void ConditionQ(REAL **x, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {

  int i, iptr, k, ne, nf, nc, kmin, warn=0;
  REAL *a = phys->a;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      x[i][k]=0;

    for(nf=0;nf<grid->nfaces[i];nf++) 
      if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1) {

        ne = grid->face[i*grid->maxfaces+nf];

        if(grid->ctop[nc]>grid->ctop[i])
          kmin = grid->ctop[nc];
        else
          kmin = grid->ctop[i];

        for(k=kmin;k<grid->Nke[ne];k++) 
          x[i][k]+=grid->dzz[i][k]*phys->D[ne];
      }

    a[grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]];
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
      a[k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1]);

    for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++) 
      x[i][k]+=(a[k]+a[k+1]);

    if(grid->ctop[i]<grid->Nk[i]-1) {
      // Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
      k=grid->ctop[i];
      x[i][k]+=2*a[k]+a[k+1];

      // Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
      k=grid->Nk[i]-1;
      x[i][k]+=a[k];
    }
  }

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      if(x[i][k]<=0) {
        x[i][k]=1;
        warn=1;
      }
      x[i][k]=sqrt(x[i][k]);
    }
  }
  if(WARNING && warn) printf("Warning...invalid preconditioner!\n");

  // Send the preconditioner to the neighboring processors.
  ISendRecvCellData3D(x,grid,myproc,comm);
}

/*
 * Function: Preconditioner
 * Usage: Preconditioner(x,xc,coef,grid,phys,prop);
 * ------------------------------------------------
 * Multiply the vector x by the inverse of the preconditioner M with
 * xc = M^{-1} x
 *
 */
static void Preconditioner(REAL **x, REAL **xc, REAL **coef, gridT *grid, physT *phys, propT *prop) {
  int i, iptr, k, nf, ne, nc, kmin;
  REAL *a = phys->a, *b = phys->b, *c = phys->c, *d = phys->d;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i=grid->cellp[iptr];

    if(grid->ctop[i]<grid->Nk[i]-1) {
      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++) {
        a[k]=coef[i][k];
        b[k]=-coef[i][k]-coef[i][k+1];
        c[k]=coef[i][k+1];
        d[k]=x[i][k];
      }

      // Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
      k=grid->ctop[i];
      b[k]=-2*coef[i][k]-coef[i][k+1];
      c[k]=coef[i][k+1];
      d[k]=x[i][k];

      // Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
      k=grid->Nk[i]-1;
      a[k]=coef[i][k];
      b[k]=-coef[i][k];
      d[k]=x[i][k];

      TriSolve(&(a[grid->ctop[i]]),&(b[grid->ctop[i]]),&(c[grid->ctop[i]]),
          &(d[grid->ctop[i]]),&(xc[i][grid->ctop[i]]),grid->Nk[i]-grid->ctop[i]);
      //      for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      //	xc[i][k]=x[i][k];
    } else 
      xc[i][grid->ctop[i]]=-0.5*x[i][grid->ctop[i]]/coef[i][grid->ctop[i]];
  }
}

/*
 * Function: OperatorQC
 * Usage: OperatorQC(coef,fcoef,x,y,c,grid,phys,prop);
 * ---------------------------------------------------
 * Given a vector x, computes the left hand side of the nonhydrostatic pressure
 * Poisson equation and places it into y with y = L(x) for the preconditioned
 * solver.
 *
 * The coef array contains coefficients for the vertical derivative terms in the operator
 * while the fcoef array contains coefficients for the horizontal derivative terms.  These
 * are computed before the iteration in QCoefficients. The array c stores the preconditioner.
 *
 */
static void OperatorQC(REAL **coef, REAL **fcoef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop) {
  int i, iptr, k, ne, nf, nc, kmin, kmax;
  REAL *a = phys->a;

  // sum over all computational cells
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    // over the depth of defined cells
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      y[i][k]=-x[i][k];

    // over each face
    for(nf=0;nf<grid->nfaces[i];nf++) {
      if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1) {

        ne = grid->face[i*grid->maxfaces+nf];

        if(grid->ctop[nc]>grid->ctop[i])
          kmin = grid->ctop[nc];
        else
          kmin = grid->ctop[i];

        // now sum over depth too
        for(k=kmin;k<grid->Nke[ne];k++) 
          y[i][k]+=x[nc][k]*fcoef[i*grid->maxfaces+nf][k];
      }
    }

    // over depth now
    for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
      y[i][k]+=coef[i][k]*x[i][k-1]+coef[i][k+1]*x[i][k+1];

    // apply top/bottom bcs in vertical since we have no flux on sides
    if(grid->ctop[i]<grid->Nk[i]-1) {
      // Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
      k=grid->ctop[i];
      y[i][k]+=coef[i][k+1]*x[i][k+1];

      // Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
      k=grid->Nk[i]-1;
      y[i][k]+=coef[i][k]*x[i][k-1];
    }
  }
}

/*
 * Function: QCoefficients
 * Usage: QCoefficients(coef,fcoef,c,grid,phys,prop);
 * --------------------------------------------------
 * Compute coefficients for the pressure-Poisson equation.  fcoef stores
 * coefficients at the vertical flux faces while coef stores coefficients
 * at the horizontal faces for vertical derivatives of q.
 *
 */
static void QCoefficients(REAL **coef, REAL **fcoef, REAL **c, gridT *grid, 
    physT *phys, propT *prop) {

  int i, iptr, k, kmin, nf, nc, ne;

  // if we want to use the preconditioner
  if(prop->qprecond==1) {
    // over all computational cells
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      coef[i][grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]]/c[i][grid->ctop[i]];
      if(prop->subgrid)
        coef[i][grid->ctop[i]]=subgrid->Acveff[i][grid->ctop[i]]/grid->dzz[i][grid->ctop[i]]/c[i][grid->ctop[i]];
      // over all the compuational cell depths
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
        if(!prop->subgrid)
          // compute the cell coefficients
          coef[i][k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1])/(c[i][k]*c[i][k-1]);
        else
          coef[i][k] = 2*subgrid->Acveff[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1])/(c[i][k]*c[i][k-1]);          
      // over all the faces 
      for(nf=0;nf<grid->nfaces[i];nf++) 
        if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1) {

          ne = grid->face[i*grid->maxfaces+nf];

          if(grid->ctop[nc]>grid->ctop[i])
            kmin = grid->ctop[nc];
          else
            kmin = grid->ctop[i];

          // over the edge depths
          for(k=kmin;k<grid->Nke[ne];k++) 
            // compute the face coefficients
            fcoef[i*grid->maxfaces+nf][k]=grid->dzz[i][k]*phys->D[ne]/(c[i][k]*c[nc][k]);
        }
    }
  }
  // no preconditioner
  else {
    /*for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++)  */
    // over each and every cell
    for(i=0;i<grid->Nc;i++) {
      //      i = grid->cellp[iptr];

      // compute the cell coeficient for the top cell
      coef[i][grid->ctop[i]]=grid->Ac[i]/grid->dzz[i][grid->ctop[i]];
      if(prop->subgrid)
        coef[i][grid->ctop[i]]=subgrid->Acveff[i][grid->ctop[i]]/grid->dzz[i][grid->ctop[i]];
      // compute the coefficients for the rest of the cells (towards bottom)
      // where dzz is averaged to the face 
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) 
        if(!prop->subgrid)
          // compute the cell coefficients
          coef[i][k] = 2*grid->Ac[i]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
        else
          coef[i][k] = 2*subgrid->Acveff[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1]);      
    }
  }
  // where is fcoef solved for???
}

/*
 * Function: OperatorQ
 * Usage: OperatorQ(coef,x,y,c,grid,phys,prop);
 * --------------------------------------------
 * Given a vector x, computes the left hand side of the nonhydrostatic pressure
 * Poisson equation and places it into y with y = L(x) for the non-preconditioned
 * solver.
 *
 * The coef array contains coefficients for the vertical derivative terms in the operator.
 * This is computed before the iteration in QCoefficients. The array c stores the preconditioner.
 * The preconditioner stored in c is not used.
 *
 */
// This one is for debugging to check whether OperatorQSlope() works (with no send/recv)
static void OperatorQ0(REAL **coef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop) {
  MPI_Comm comm;
  int myproc;
  OperatorQSlope(x,y,grid,phys,prop,myproc,comm);
}
static void OperatorQ(REAL **coef, REAL **x, REAL **y, REAL **c, gridT *grid, physT *phys, propT *prop) {

  int i, iptr, k, ne, nf, nc, kmin, kmax;

  // over each computational cell
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    // over cells that exist and aren't cut off by bathymetry
    for(k=grid->ctop[i];k<grid->Nk[i];k++) 
      y[i][k]=0;

    // over each face
    for(nf=0;nf<grid->nfaces[i];nf++) 
      // we only apply this to non-boundary cells 
      // (meaning phys->D implicitly 0 at this point)
      if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1) {

        ne = grid->face[i*grid->maxfaces+nf];

        // determine the minimal k based on each cell over edge
        if(grid->ctop[nc]>grid->ctop[i])
          kmin = grid->ctop[nc];
        else
          kmin = grid->ctop[i];

        // create summed contributions for cell along the vertical
        for(k=kmin;k<grid->Nke[ne];k++) 
          y[i][k]+=(x[nc][k]-x[i][k])*grid->dzf[ne][k]*phys->D[ne];
      }
    // otherwise y[i][...] = 0 (boundary cell with no gradient in horiz.)

    // over all the depth (z dependence on Laplacian)
    for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
      // we are adding on the vertical diffusional part 
      y[i][k]+=coef[i][k]*x[i][k-1]-(coef[i][k]+coef[i][k+1])*x[i][k]+coef[i][k+1]*x[i][k+1];

    // apply top and bottom BCs in vertical, no flux on sides
    if(grid->ctop[i]<grid->Nk[i]-1) {
      // Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
      k=grid->ctop[i];
      y[i][k]+=(-2*coef[i][k]-coef[i][k+1])*x[i][k]+coef[i][k+1]*x[i][k+1];

      // Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
      k=grid->Nk[i]-1;
      y[i][k]+=coef[i][k]*x[i][k-1]-coef[i][k]*x[i][k];
    } 
    else
      y[i][grid->ctop[i]]-=2.0*coef[i][grid->ctop[i]]*x[i][grid->ctop[i]];
  }
}

/*
 * Function: Compute_dqdxi_Face
 * Usage: Compute_dqdxi_Face(dqdxi_face,x,grid,phys,myproc,comm);
 * --------------------------------------------------------------
 * Compute 1/J dq/dxi on the faces. This requires first computing dq/dxi on the cell centers with
 * vertical interpolation and then interpolating to the faces with horizontal interpolation.
 *
 */
void Compute_dqdxi_Face(REAL **dqdxi_face, REAL **x, gridT *grid, physT *phys, int myproc, MPI_Comm comm) {
  int i, iptr, j, jptr, k, nc1, nc2;
  REAL def1, def2, dzz_kmh, dzz_kph, dqdxi_kmh, dqdxi_kph, q_kmh, q_kph;
  REAL slope_mag_squared, dqdx_kph, dqdy_kph;
  REAL **dqdxi_Center = phys->stmp4; // Temporary array
  REAL *dzh = SunMalloc((grid->Nkmax+1)*sizeof(REAL),"Compute_dqdxi_Face"), *dz;
  
  // First comptue 1/J dq/dxi at cell centers
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    // Delta z at the w faces
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
      dzh[k]=0.5*(grid->dzz[i][k-1]+grid->dzz[i][k]);
    dzh[grid->ctop[i]]=grid->dzz[i][grid->ctop[i]];
    dzh[grid->Nk[i]]=grid->dzz[i][grid->Nk[i]-1];    

    // Delta z at cell centers (easier/cleaner to access dz[k] than grid->dzz[i][k])
    dz = grid->dzz[i];
    
    for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++) {
      q_kmh = (dz[k]*x[i][k-1]+dz[k-1]*x[i][k])/(dz[k]+dz[k-1]);   // q at k-1/2
      q_kph = (dz[k]*x[i][k+1]+dz[k+1]*x[i][k])/(dz[k]+dz[k+1]);   // q at k+1/2

      // Derivative at center is based on top and bottom faces
      // Note that for equispaced grids this is central differencing (q_{k-1} - q_{k+1})/(2 dz)
      dqdxi_Center[i][k]=(q_kmh - q_kph)/dz[k];
    }

    // Top boundary q=0. Same as above but with q_kmh = 0.
    k=grid->ctop[i];
    q_kmh = 0; // q at surface
    q_kph = (dz[k]*x[i][k+1]+dz[k+1]*x[i][k])/(dz[k]+dz[k+1]);   // q at k+1/2    
    dqdxi_Center[i][k]=(q_kmh - q_kph)/dz[k];

    // Bottom boundary is averaged dq/dxi at top and bottom face, but at bottom face,
    // dq/dxi = (dzdx q_x + dzdy q_y)/(1+slope^2)
    k=grid->Nk[i]-1;
    dqdxi_kmh = (x[i][k-1]-x[i][k])/dzh[k];

    // (dz/dx)^2 + (dz/dy)^2 defined at bottom face
    slope_mag_squared =  pow(vert->dzdx[i][k+1],2.0)+pow(vert->dzdy[i][k+1],2.0);
    
    // Need to extrapolate to obtain dqdx and dqdy at the bottom face since
    // dq/dx and dq/dy are defined at cell centers. This is AB2 for equispaced grids,
    // i.e. dqdx_kph = 1.5*dqdx_k-0.5*dqdx_{k-1}
    dqdx_kph = vert->dqdx[i][k]*(2.0*dz[k]+dz[k-1])/(dz[k]+dz[k-1])-
      vert->dqdx[i][k-1]*dz[k]/(dz[k]+dz[k-1]);
    dqdy_kph = vert->dqdy[i][k]*(2.0*dz[k]+dz[k-1])/(dz[k]+dz[k-1])-
      vert->dqdy[i][k-1]*dz[k]/(dz[k]+dz[k-1]);
    dqdxi_kph = (vert->dzdx[i][k+1]*dqdx_kph+vert->dzdy[i][k+1]*dqdy_kph)/(1+slope_mag_squared);
   
    dqdxi_Center[i][k]=0.5*(dqdxi_kmh + dqdxi_kph);
  }

  // Comment this to set dqdxi_Center = 0
  /*
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      dqdxi_Center[i][k]=0;
    }
  }
  */
  ISendRecvCellData3D(dqdxi_Center, grid, myproc, comm);  
  
  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
    j = grid->edgep[jptr];
    nc1 = grid->grad[2*j];
    nc2 = grid->grad[2*j+1];

    if(nc2==-1)
      nc2=nc1;
    if(nc1==-1) { printf("ERROR in compute_dqdxi_face!\n"); exit(1); }

    def1=0.5*grid->dg[j];
    def2=0.5*grid->dg[j];

    // Horizontal interpolation to get face value from cell centers
    for(k=0;k<grid->Nkmax;k++)
      dqdxi_face[j][k]=(def2*dqdxi_Center[nc1][k]+def1*dqdxi_Center[nc2][k])/grid->dg[j];
  }

  SunFree(dzh, (grid->Nkmax+1)*sizeof(REAL),"Compute_dqdxi_Face");
}

/*
 * Function: OperatorQSlope
 * Usage: OperatorQSlope(x,z,grid,phys,prop,myproc,comm);
 * ------------------------------------------------------
 * Compute the Laplacian of x and output it into y, i.e. y = L(x), where
 *
 * L(x) = Sum_{faces} (dq/dn - 1/J dz/dn dq/dxi)_{face} h_f df N
 *      + A_c (F_{k-1/2} - F_{k+1/2})
 *      + A_c/J_{k-1/2} (1 + s_x^2 + s_y^2)_{k-1/2} (x_{k-1} - x_k)
 *      - A_c/J_{k+1/2} (1 + s_x^2 + s_y^2)_{k+1/2} (x_{k} - x_{k+1})
 *
 * where F_{k-1/2} = - dz/dx dq/dx - dz/dy dq/dy
 *
 */
static void OperatorQSlope(REAL **x, REAL **y, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int i, j, iptr, jptr, k, ne, nf, nc, nc1, nc2;
  REAL dqdx_kmh, dqdy_kmh;
  REAL *dzh = SunMalloc((grid->Nkmax+1)*sizeof(REAL),"OperatorQSlope");
  REAL *Fq =  SunMalloc((grid->Nkmax+1)*sizeof(REAL),"OperatorQSlope");
  REAL *slope_mag =  SunMalloc((grid->Nkmax+1)*sizeof(REAL),"OperatorQSlope");  
  REAL **dqdxi_face = vert->uf;

  if(prop->include_slope_terms) {
    // Compute the cell-centered gradients dqdx and dqdy
    // The first integer is the direction (0=x,1=y), the second zero means this is cell-centered (1=w-centered)    
    ComputeCellAveragedHorizontalGradientCell(vert->dqdx,0,x,0,grid,prop,phys,myproc);
    ComputeCellAveragedHorizontalGradientCell(vert->dqdy,1,x,0,grid,prop,phys,myproc);

    // Compute 1/J dq/dxi at the faces. This is not straightforward because it needs vertical
    // and horizontal interpolation
    Compute_dqdxi_Face(dqdxi_face,x,grid,phys,myproc,comm);
  }
  
  // Compute the Laplacian y = L(x)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    for(k=0;k<grid->Nk[i];k++) 
      y[i][k]=0;
    
    for(nf=0;nf<grid->nfaces[i];nf++) {
      // Only at non-boundary cells
      if((nc=grid->neigh[i*grid->maxfaces+nf])!=-1) {
	
        ne = grid->face[i*grid->maxfaces+nf];
	
        // create summed contributions for cell along the vertical
        for(k=grid->etop[ne];k<grid->Nke[ne];k++) 
          y[i][k]+=(x[nc][k]-x[i][k])*grid->dzf[ne][k]*grid->df[ne]/grid->dg[ne];
	// otherwise y[i][...] = 0 (boundary cell with no gradient in horiz.)
	
	if(prop->include_slope_terms) {
	  for(k=grid->etop[ne];k<grid->Nke[ne];k++) {
	    y[i][k]-=dqdxi_face[ne][k]*(vert->zc[nc][k]-vert->zc[i][k])*grid->dzf[ne][k]*grid->df[ne]/grid->dg[ne];
	  }
	}
      }
    }

    // Compute the half layer heights
    for(k=grid->ctop[i]+1;k<grid->Nk[i];k++)
      dzh[k]=0.5*(grid->dzz[i][k-1]+grid->dzz[i][k]);
    dzh[grid->ctop[i]]=grid->dzz[i][grid->ctop[i]];
    dzh[grid->Nk[i]]=grid->dzz[i][grid->Nk[i]-1];

    if(prop->include_slope_terms) {

      // Compute (dz/dx)^2 + (dz/dy)^2
      for(k=grid->ctop[i];k<=grid->Nk[i];k++) {
	slope_mag[k]=pow(vert->dzdx[i][k],2.0)+pow(vert->dzdy[i][k],2.0);
      }

      // Interpolate dq/dx and dqdy from cell centers to w faces "kmh" = k-1/2
      for(k=grid->ctop[i]+1;k<grid->Nk[i];k++) {
	dqdx_kmh = vert->dqdx[i][k]*grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])+
	  vert->dqdx[i][k-1]*grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1]);
	dqdy_kmh = vert->dqdy[i][k]*grid->dzz[i][k-1]/(grid->dzz[i][k]+grid->dzz[i][k-1])+
	  vert->dqdy[i][k-1]*grid->dzz[i][k]/(grid->dzz[i][k]+grid->dzz[i][k-1]);	
	Fq[k]=-vert->dzdx[i][k]*dqdx_kmh-vert->dzdy[i][k]*dqdy_kmh;
      }

      // Extrapolate to obtain dqdx and dqdy at z=eta (essentially AB2 on non-equispaced grid)
      k=grid->ctop[i];
      dqdx_kmh = vert->dqdx[i][k]*(2.0*grid->dzz[i][k]+grid->dzz[i][k+1])/(grid->dzz[i][k]+grid->dzz[i][k+1])-
	vert->dqdx[i][k+1]*grid->dzz[i][k+1]/(grid->dzz[i][k]+grid->dzz[i][k+1]);
      dqdy_kmh = vert->dqdy[i][k]*(2.0*grid->dzz[i][k]+grid->dzz[i][k+1])/(grid->dzz[i][k]+grid->dzz[i][k+1])-
	vert->dqdy[i][k+1]*grid->dzz[i][k+1]/(grid->dzz[i][k]+grid->dzz[i][k+1]);
      Fq[k]=-vert->dzdx[i][k]*dqdx_kmh-vert->dzdy[i][k]*dqdy_kmh;

      // Bottom boundary requires Fq=0 along with dq/dxi3=0
      k=grid->Nk[i];
      Fq[k]=0;
    } else {
      for(k=grid->ctop[i];k<grid->Nk[i]+1;k++) {
	Fq[k]=slope_mag[k]=0;
      }
    }
    
    // TURN OFF SLOPE TERMS WITH THIS
    /*
    for(k=0;k<=grid->Nkmax;k++)
      Fq[k]=slope_mag[k]=0;
    */

    int vert_slope_terms=1;
    if(!vert_slope_terms) {
      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
	y[i][k]+=grid->Ac[i]*(1.0/dzh[k]*(x[i][k-1]-x[i][k])-
			      1.0/dzh[k+1]*(x[i][k]-x[i][k+1]));
      
      if(grid->ctop[i]<grid->Nk[i]-1) {
	// Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
	k=grid->ctop[i];
	y[i][k]+=grid->Ac[i]*(-2.0/dzh[k]*x[i][k]-1.0/dzh[k+1]*(x[i][k]-x[i][k+1]));      
	
	// Bottom dq/dz = 0 so q[i][grid->Nk[i]]=q[i][grid->Nk[i]-1]
	k=grid->Nk[i]-1;
	y[i][k]+=grid->Ac[i]/dzh[k]*(x[i][k-1]-x[i][k]);
      } else {
	k=grid->ctop[i];
	y[i][grid->ctop[i]]-=2.0*grid->Ac[i]/dzh[grid->ctop[i]]*x[i][grid->ctop[i]];
      }
    } else {
      for(k=grid->ctop[i]+1;k<grid->Nk[i]-1;k++)
	y[i][k]+=grid->Ac[i]*((1.0+slope_mag[k])/dzh[k]*(x[i][k-1]-x[i][k])-
			      (1.0+slope_mag[k+1])/dzh[k+1]*(x[i][k]-x[i][k+1])+
			      (Fq[k]-Fq[k+1]));
      
      // apply top and bottom BCs in vertical, no flux on sides
      if(grid->ctop[i]<grid->Nk[i]-1) {
	// Top q=0 so q[i][grid->ctop[i]-1]=-q[i][grid->ctop[i]]
	k=grid->ctop[i];
	y[i][k]+=grid->Ac[i]*(-2.0*(1.0+slope_mag[k])/dzh[k]*x[i][k]-
			      (1.0+slope_mag[k+1])/dzh[k+1]*(x[i][k]-x[i][k+1])+
			      (Fq[k]-Fq[k+1]));
	
	// Bottom F + (1+S^2)dq/dz = 0 so set F and (1+S^2)dq/dz to 0
	k=grid->Nk[i]-1;
	y[i][k]+=grid->Ac[i]*((1.0+slope_mag[k])/dzh[k]*(x[i][k-1]-x[i][k])+Fq[k]);
      } else {
	k=grid->ctop[i];
	y[i][k]+=grid->Ac[i]*(-2.0*(1.0+slope_mag[k])/dzh[k]*x[i][k]+Fq[k]);
      }
    }
  }
 

  SunFree(dzh,(grid->Nkmax+1)*sizeof(REAL),"OperatorQSlope");
  SunFree(Fq,(grid->Nkmax+1)*sizeof(REAL),"OperatorQSlope");
  SunFree(slope_mag,(grid->Nkmax+1)*sizeof(REAL),"OperatorQSlope");    
}

/*
 * Function: InnerProduct3
 * Usage: InnerProduct3(x,y,grid,myproc,numprocs,comm);
 * ---------------------------------------------------
 * Compute the inner product of two two-dimensional arrays x and y.
 * Used for the CG method to solve for the nonhydrostatic pressure.
 *
 */
static REAL InnerProduct3(REAL **x, REAL **y, gridT *grid, int myproc, int numprocs, MPI_Comm comm) {

  int i, k, iptr;
  REAL sum, mysum=0;

  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];

    for(k=grid->ctop[i];k<grid->Nk[i];k++)
      mysum+=x[i][k]*y[i][k];
  }
  MPI_Reduce(&mysum,&(sum),1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Bcast(&sum,1,MPI_DOUBLE,0,comm);

  return sum;
}  
