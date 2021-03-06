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
  LIS_INT		nthreads,maxthreads;
  LIS_INT		gn,nnz,np,mode;
  LIS_INT		i,ii,j,jj,k,kk,ctr;
  LIS_INT		l,m,n,nn;
  LIS_INT		rn,rmin,rmax,rb;
  LIS_INT		is,ie,clsize,ci,*iw;
  LIS_INT		err,iter,matrix_type,s,ss,se;
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

  if( argc < 5 )
    {
      if( my_rank==0 ) printf("Usage: spmvtest3 l m n iter [matrix_type]\n");
      CHKERR(1);
    }

  l  = atoi(argv[1]);
  m  = atoi(argv[2]);
  n  = atoi(argv[3]);
  iter = atoi(argv[4]);
  if (argv[5] == NULL) {
    s = 0;
  }
  else {
    s = atoi(argv[5]);
  }

  if( iter<=0 )
    {
      if( my_rank==0 ) printf("iter=%d <= 0\n",iter);
      CHKERR(1);
    }
  if( l<=0 || m<=0 || n<=0 )
    {
      if( my_rank==0 ) printf("l=%d <=0, m=%d <=0 or n=%d <=0\n",l,m,n);
      CHKERR(1);
    }
  if( s<0 || s>11 )
    {
      if( my_rank==0 ) printf("matrix_type=%d < 0 or matrix_type=%d > 11\n",s,s);
      CHKERR(1);
    }

  if( my_rank==0 )
    {
      printf("\n");
      printf("number of processes = %d\n",nprocs);
    }

#ifdef _OPENMP
  nthreads = omp_get_num_procs();
  maxthreads = omp_get_max_threads();
  if( my_rank==0 )
    {
      printf("max number of threads = %d\n", nthreads);
      printf("number of threads = %d\n", maxthreads);
    }
#else
  nthreads = 1;
  maxthreads = 1;
#endif

  /* create matrix and vectors */
  nn = l*m*n;
  err = lis_matrix_create(LIS_COMM_WORLD,&A0);
  err = lis_matrix_set_size(A0,0,nn);
  CHKERR(err);

  ptr   = (LIS_INT *)malloc((A0->n+1)*sizeof(LIS_INT));
  if( ptr==NULL ) CHKERR(1);
  index = (LIS_INT *)malloc(7*A0->n*sizeof(LIS_INT));
  if( index==NULL ) CHKERR(1);
  value = (LIS_SCALAR *)malloc(7*A0->n*sizeof(LIS_SCALAR));
  if( value==NULL ) CHKERR(1);

  lis_matrix_get_range(A0,&is,&ie);
  ctr = 0;
  for(ii=is;ii<ie;ii++)
    {
      i = ii/(m*n);
      jj = ii - i*(m*n);
      j = jj/n;
      k = jj - j*n;
      if( i>0 )   { kk = ii - m*n; index[ctr] = kk; value[ctr++] = -1.0;}
      if( i<l-1 ) { kk = ii + m*n; index[ctr] = kk; value[ctr++] = -1.0;}
      if( j>0 )   { kk = ii - n; index[ctr] = kk; value[ctr++] = -1.0;}
      if( j<m-1 ) { kk = ii + n; index[ctr] = kk; value[ctr++] = -1.0;}
      if( k>0 )   { kk = ii - 1; index[ctr] = kk; value[ctr++] = -1.0;}
      if( k<n-1 ) { kk = ii + 1; index[ctr] = kk; value[ctr++] = -1.0;}
      index[ctr] = ii; value[ctr++] = 6.0;
      ptr[ii-is+1] = ctr;
    }
  ptr[0] = 0;
  err = lis_matrix_set_crs(ptr[ie-is],ptr,index,value,A0);
  CHKERR(err);
  err = lis_matrix_assemble(A0);
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
#ifdef _LONGLONG
      printf("matrix size = %lld x %lld (%lld nonzero entries)\n",gn,gn,nnz);
#else
      printf("matrix size = %d x %d (%d nonzero entries)\n",gn,gn,nnz);
#endif
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
  for(i=0;i<n;i++)
    {
      lis_sort_id(A0->ptr[i],A0->ptr[i+1]-1,A0->index,A0->value);
    }
		
  /* 
     Parallel version of VBR is excluded temporally due to bug.
     DNS is also excluded to reduce memory usage.
  */

  if (s==0) 
    {
      ss = 1;
      se = 11;
    }
  else
    {
      ss = s;
      se = s+1;
    }
	
  for (matrix_type=ss;matrix_type<se;matrix_type++)
    {
      if ( (nprocs>1 || nthreads>1) && matrix_type==9 ) continue;
      lis_matrix_duplicate(A0,&A);
      lis_matrix_set_type(A,matrix_type);
      err = lis_matrix_convert(A0,A);
      if( err ) CHKERR(err);
		    
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
#ifdef USE_MPI
	  printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS, communication = %e sec, communication/computation = %3.3f %%, 2-norm = %e\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,commtime,commtime/comptime*100,val);
#else
	  printf("matrix_type = %2d (%s), computation = %e sec, %8.3f MFLOPS, 2-norm = %e\n",matrix_type,lis_storagename2[matrix_type-1],comptime,flops,val);
#endif
	}
      lis_matrix_destroy(A);
    }

  lis_matrix_destroy(A0);
  lis_vector_destroy(b);
  lis_vector_destroy(x);

  lis_finalize();

  LIS_DEBUG_FUNC_OUT;

  return 0;
}


