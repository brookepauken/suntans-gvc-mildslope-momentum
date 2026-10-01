/*
 * File: partition.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains functions that use the ParMetis libraries
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "partition.h"
#include "memory.h"
#include "grid.h"

//Parmetis 2.0
//#include "parmetis.h"
//Parmetis 3.1
#include "parmetislib.h"

// Private function
static void GetGraph(GraphType *graph, gridT *grid, MPI_Comm comm);

/*
 * Function: GetPartitioning
 * Usage: GetPartitioning(maingrid,localgrid,myproc,numprocs,comm);
 * ----------------------------------------------------------------
 * This function uses the ParMetis libraries to compute the grid partitioning and places
 * the partition number into the maingrid->part array.
 *
 */
void GetPartitioning(gridT *maingrid, gridT **localgrid, int myproc, int numprocs, MPI_Comm comm) {
  int j, proc, numflag=0, wgtflag=0, options[5], edgecut;
  GraphType graph;
  MPI_Status status;

  //if doing 2d (x-z) locally, set to 1 to avoid the parmetis paramaterization
  const int STRIPE_PARTITION = 1;

  if(numprocs>1) {
    if(STRIPE_PARTITION) {
      if(myproc==0 && VERBOSE>1) printf("\tStripe partitioning by xv...\n");
      if(myproc==0) {
        int Nc = maingrid->Nc;
        int *idx = (int *)SunMalloc(Nc*sizeof(int),"GetPartitioning");
        for(j=0;j<Nc;j++) idx[j] = j;
        /* insertion sort by xv (Nc small enough for this ~10^4 cells; for
         * larger Nc swap to qsort with a comparator capturing xv). */
        for(int a=1;a<Nc;a++) {
          int ka = idx[a];
          REAL xa = maingrid->xv[ka];
          int b = a-1;
          while(b>=0 && maingrid->xv[idx[b]] > xa) { idx[b+1]=idx[b]; b--; }
          idx[b+1] = ka;
        }
        for(j=0;j<Nc;j++) {
          int rank = (int)((long long)j*numprocs/Nc);
          if(rank>=numprocs) rank = numprocs-1;
          maingrid->part[idx[j]] = rank;
        }
        SunFree(idx,Nc*sizeof(int),"GetPartitioning");
      }
      MPI_Bcast((void *)maingrid->part,maingrid->Nc,MPI_INT,0,comm);
    } else {
      options[0] = 0;
      wgtflag = 2;
      numflag = 0;

      GetGraph(&graph,maingrid,comm);

      (*localgrid)->part = SunMalloc(graph.nvtxs*sizeof(int),"GetPartitioning");

      if(myproc==0 && VERBOSE>2) printf("Partitioning with ParMETIS_PartKway...\n");
      ParMETIS_PartKway(graph.vtxdist,graph.xadj,graph.adjncy,graph.vwgt,NULL,
                        &wgtflag,&numflag,&numprocs,options,&edgecut,
                        (*localgrid)->part,&comm);

      if(myproc==0 && VERBOSE>2) printf("Redistributing the partition arrays...\n");
      if(myproc!=0)
        MPI_Send((void *)(*localgrid)->part,graph.nvtxs,MPI_INT,0,1,comm);
      else {
        for(j=0;j<graph.nvtxs;j++)
          maingrid->part[graph.vtxdist[myproc]+j]=(*localgrid)->part[j];
        for(proc=1;proc<numprocs;proc++)
          MPI_Recv((void *)(maingrid->part+graph.vtxdist[proc]),
                   graph.vtxdist[proc+1]-graph.vtxdist[proc],MPI_INT,proc,1,comm,&status);
      }
      MPI_Bcast((void *)maingrid->part,maingrid->Nc,MPI_INT,0,comm);

      SunFree((*localgrid)->part,graph.nvtxs*sizeof(int),"GetPartitioning");
    }
  } else {
    for(j=0;j<maingrid->Nc;j++)
      maingrid->part[j]=0;
  }
}

/*
 * Function: GetGraph
 * Usage GetGraph(graph,grid,comm);
 * --------------------------------
 * This code was adapted from the ParMetis-2.0 code in the ParMetis
 * distribution.
 *
 */
static void GetGraph(GraphType *graph, gridT *grid, MPI_Comm comm)
{
  int i, k, l, numprocs, myproc;
  int nvtxs, penum, snvtxs;
  idxtype *gxadj, *gadjncy, *gvwgt;  
  idxtype *vtxdist, *sxadj, *ssize, *svwgt;
  MPI_Status status;

  MPI_Comm_size(comm, &numprocs);
  MPI_Comm_rank(comm, &myproc);

  vtxdist = graph->vtxdist = idxsmalloc(numprocs+1, 0, "ReadGraph: vtxdist");

  if (myproc == 0) {
    ssize = idxsmalloc(numprocs, 0, "ReadGraph: ssize");

    nvtxs = grid->Nc;
    gxadj = grid->xadj;
    gadjncy = grid->adjncy;
    gvwgt = (idxtype *)grid->vwgt;

    /* Construct vtxdist and send it to all the processors */
    vtxdist[0] = 0;
    for (i=0,k=nvtxs; i<numprocs; i++) {
      l = k/(numprocs-i);
      vtxdist[i+1] = vtxdist[i]+l;
      k -= l;
    }
  }

  MPI_Bcast((void *)vtxdist, numprocs+1, IDX_DATATYPE, 0, comm);

  graph->gnvtxs = vtxdist[numprocs];
  graph->nvtxs = vtxdist[myproc+1]-vtxdist[myproc];
  graph->xadj = idxmalloc(graph->nvtxs+1, "ReadGraph: xadj");
  graph->vwgt = idxmalloc(graph->nvtxs, "ReadGraph: vwgt");

  if (myproc == 0) {
    for (penum=0; penum<numprocs; penum++) {
      snvtxs = vtxdist[penum+1]-vtxdist[penum];
      sxadj = idxmalloc(snvtxs+1, "ReadGraph: sxadj");
      svwgt = idxmalloc(snvtxs, "ReadGraph: svwgt");

      idxcopy(snvtxs+1, gxadj+vtxdist[penum], sxadj);
      idxcopy(snvtxs, gvwgt+vtxdist[penum], svwgt);

      if(VERBOSE>3) 
	for(k=0;k<snvtxs;k++)
	  printf("svwgt[%d]=%d\n",k,svwgt[k]);

      for (i=snvtxs; i>=0; i--)
        sxadj[i] -= sxadj[0];

      ssize[penum] = gxadj[vtxdist[penum+1]] - gxadj[vtxdist[penum]];

      if (penum == myproc) {
        idxcopy(snvtxs+1, sxadj, graph->xadj);
	idxcopy(snvtxs, svwgt, graph->vwgt);
      } else {
        MPI_Send((void *)sxadj, snvtxs+1, IDX_DATATYPE, penum, 1, comm); 
        MPI_Send((void *)svwgt, snvtxs, IDX_DATATYPE, penum, 1, comm); 
      }

      free(sxadj);
      free(svwgt);
    }
  }
  else {
    MPI_Recv((void *)graph->xadj, graph->nvtxs+1, IDX_DATATYPE, 0, 1, comm, &status);
    MPI_Recv((void *)graph->vwgt, graph->nvtxs, IDX_DATATYPE, 0, 1, comm, &status);
  }
  if(VERBOSE>3) {
    printf("Weights on each processor after MPI_Recv\n");
    for(k=0;k<graph->nvtxs;k++)
      printf("vwgt[%d]=%d\n",k,graph->vwgt[k]);
  }

  graph->nedges = graph->xadj[graph->nvtxs];
  graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: graph->adjncy");

  if (myproc == 0) {
    for (penum=0; penum<numprocs; penum++) {
      if (penum == myproc) 
        idxcopy(ssize[penum], gadjncy+gxadj[vtxdist[penum]], graph->adjncy);
      else
        MPI_Send((void *)(gadjncy+gxadj[vtxdist[penum]]), ssize[penum], IDX_DATATYPE, penum, 1, comm); 
    }

    free(ssize);
  }
  else 
    MPI_Recv((void *)graph->adjncy, graph->nedges, IDX_DATATYPE, 0, 1, comm, &status);

  if (myproc == 0) 
    GKfree(&gxadj, &gadjncy, LTERM);

  MALLOC_CHECK(NULL);
}



