/* Copyright (C) 2002-2012 The SSI Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
        #include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN32_H
        #include "lis_config_win32.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#ifdef USE_MPI
	#include <mpi.h>
#endif
#include "lislib.h"

char *lis_storagename2[]   = {"CRS", "CCS", "MSR", "DIA", "ELL", "JDS", "BSR", "BSC", "VBR", "COO", "DNS"};

#undef __FUNC__
#define __FUNC__ "main"
LIS_INT main(LIS_INT argc, char* argv[])
{
  LIS_MATRIX		A,A0;
  LIS_VECTOR		b,x,v;
  LIS_SCALAR		ntimes,nmflops,nnrm2;
  LIS_SCALAR		*value;

  LIS_INT		nprocs,my_rank;
  int    		int_nprocs,int_my_rank;
  LIS_INT		nthreads, maxthreads;
  LIS_INT		gn,nnz,mode;
  LIS_INT		i,j,jj,j0,j1,l,k,n,np,h,ih;
  LIS_INT		m,nn,ii;
  LIS_INT		block;
  LIS_INT		rn,rmin,rmax,rb;
  LIS_INT		is,ie,clsize,ci,*iw;
  LIS_INT		err,iter,matrix_type;
  LIS_INT	       	*ptr,*index;
  double		mem,val,ra,rs,ri,ria,ca,time,time2,convtime,val2,nnzs,nnzap,nnzt;
  double		commtime,comptime,flops;
  FILE			*file;
  char path[1024];

  LIS_DEBUG_FUNC_IN;
    
  lis_initialize(&argc, &argv);

#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD,&int_nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&int_my_rank);
  nprocs = int_nprocs;
  my_rank = int_my_rank;
#else
  nprocs  = 1;
  my_rank = 0;
#endif

  if( argc < 4 )
    {
      if( my_rank==0 ) printf("Usage: spmvtest5 matrix_filename matrix_type iter [block] \n");
      lis_finalize();
      exit(0);
    }

  file = fopen(argv[1], "r");
  if( file==NULL ) CHKERR(1);

  matrix_type  = atoi(argv[2]);
  iter = atoi(argv[3]);
  if (argv[4] == NULL) {
    block = 2;
  }
  else {
    block = atoi(argv[4]);
  }

  if( matrix_type<1 || matrix_type>11 )
    {
      if( my_rank==0 ) printf("matrix_type=%d <1 or matrix_type=%d >11\n",matrix_type,matrix_type);
      CHKERR(1);
    }
  if( iter<=0 )
    {
      if( my_rank==0 ) printf("iter=%d <= 0\n",iter);
      CHKERR(1);
    }

  if( my_rank==0 )
    {
      printf("\n");
      printf("number of processes = %d\n",nprocs);
    }

#ifdef _OPENMP
  if( my_rank==0 )
    {
      nthreads = omp_get_num_procs();
      maxthreads = omp_get_max_threads();
      printf("max number of threads = %d\n", nthreads);
      printf("number of threads = %d\n", maxthreads);
    }
#else
      nthreads = 1;
      maxthreads = 1;
#endif

  /* create matrix and vectors */
  lis_matrix_create(LIS_COMM_WORLD,&A0);
  err = lis_input(A0,NULL,NULL,argv[1]);
  CHKERR(err);

  n   = A0->n;
  gn  = A0->gn;
  nnz = A0->nnz;
  np  = A0->np-n;
#ifdef USE_MPI
  MPI_Allreduce(&nnz,&i,1,LIS_MPI_INT,MPI_SUM,A0->comm);
  nnzap = (double)i / (double)nprocs;
  nnzt  = ((double)nnz -nnzap)*((double)nnz -nnzap);
  nnz   = i;
  MPI_Allreduce(&nnzt,&nnzs,1,MPI_DOUBLE,MPI_SUM,A0->comm);
  nnzs  = (nnzs / (double)nprocs)/nnzap;
  MPI_Allreduce(&np,&i,1,LIS_MPI_INT,MPI_SUM,A0->comm);
  np = i;
#endif

  if( my_rank==0 ) 
    {
      printf("block size of BSR and BSC = %d x %d\n",block,block);
      printf("iteration count = %d\n\n",iter);
    }

  err = lis_vector_duplicate(A0,&x);
  if( err ) CHKERR(err);
  err = lis_vector_duplicate(A0,&b);
  if( err ) CHKERR(err);

  lis_matrix_get_range(A0,&is,&ie);
  for(i=0;i<n;i++)
    {
      err = lis_vector_set_value(LIS_INS_VALUE,i+is,1.0,x);
    }

  lis_matrix_duplicate(A0,&A);
  lis_matrix_set_type(A,matrix_type);
  err = lis_matrix_convert(A0,A);
  if( err ) CHKERR(err);
  if( A->matrix_type==LIS_MATRIX_BSR || A->matrix_type==LIS_MATRIX_BSC )
    {
      A->bnr = block;
      A->bnc = block;
    }
		    
  comptime = 0.0;
  commtime = 0.0;

  for(i=0;i<iter;i++)
    {
#ifdef USE_MPI
      MPI_Barrier(A->comm);
      time = lis_wtime();
      lis_send_recv(A->commtable,x->value);
      commtime += lis_wtime() - time;
#endif
      time2 = lis_wtime();
      lis_matvec(A,x,b);
      comptime += lis_wtime() - time2;
    }
  lis_vector_nrm2(b,&val);

  if( my_rank==0 )
    {
      flops = 2.0*nnz*iter*1.0e-6 / comptime;
      printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f %%, 2-norm = %e\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100,val);
    }

  lis_matrix_destroy(A);
  lis_matrix_destroy(A0);
  lis_vector_destroy(b);
  lis_vector_destroy(x);

  lis_finalize();

  LIS_DEBUG_FUNC_OUT;

  return 0;
}


