/*
 *  -- High Performance Computing Linpack Benchmark (HPL-GPU)
 *     HPL-GPU - 1.1 - 2011
 *
 *     David Rohr
 *     Matthias Kretz
 *     Matthias Bach
 *     Goethe Universität, Frankfurt am Main
 *     Frankfurt Institute for Advanced Studies
 *     (C) Copyright 2010 All Rights Reserved
 *
 *     Antoine P. Petitet
 *     University of Tennessee, Knoxville
 *     Innovative Computing Laboratory
 *     (C) Copyright 2000-2008 All Rights Reserved
 *
 *  -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 *  1. Redistributions  of  source  code  must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce  the above copyright
 *  notice, this list of conditions,  and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. All  advertising  materials  mentioning  features  or  use of this
 *  software must display the following acknowledgements:
 *  This  product  includes  software  developed  at  the  University  of
 *  Tennessee, Knoxville, Innovative Computing Laboratory.
 *  This product  includes software  developed at the Frankfurt Institute
 *  for Advanced Studies.
 *
 *  4. The name of the  University,  the name of the  Laboratory,  or the
 *  names  of  its  contributors  may  not  be used to endorse or promote
 *  products  derived   from   this  software  without  specific  written
 *  permission.
 *
 *  -- Disclaimer:
 *
 *  THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 *  OR  CONTRIBUTORS  BE  LIABLE FOR ANY  DIRECT,  INDIRECT,  INCIDENTAL,
 *  SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES  (INCLUDING,  BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,  OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * ======================================================================
 */

/*
 * Include files
 */
#include "hpl.h"
#include "util_timer.h"
#include "util_cal.h"
#ifdef HPL_GPU_TEMPERATURE_THRESHOLD
#include "util_adl.h"
#endif
/* 
 * Purpose
 * =======
 *
 * HPL_pdupdateTT broadcast - factorize and forward the panel PBCST and simultaneously
 * applies the row interchanges and updates part of the trailing (using
 * the panel PANEL) submatrix.
 *
 * Arguments
 * =========
 *
 * PBCST (local input/output) HPL_T_panel *
 * On entry, PBCST points to the data structure containing the
 * panel (to be broadcast) information.
 *
 * PANEL (local input/output) HPL_T_panel *
 * On entry, PANEL points to the data structure containing the
 * panel (to be updated) information.
 *
 * NN (local input) const int
 * On entry, NN specifies the local number of columns of the
 * trailing submatrix to be updated starting at the current
 * position. NN must be at least zero.
 *
 * ---------------------------------------------------------------------
 */ 
/* 
 * Purpose
 * =======
 *
 * HPL_pdgesv factors a N+1-by-N matrix using LU factorization with row
 * partial pivoting.  The main algorithm  is the "right looking" variant
 * with look-ahead.  The  lower  triangular factor is left unpivoted and
 * the pivots are not returned. The right hand side is the N+1 column of
 * the coefficient matrix.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPL_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * ALGO    (global input)                HPL_T_palg *
 *         On entry,  ALGO  points to  the data structure containing the
 *         algorithmic parameters.
 *
 * A       (local input/output)          HPL_T_pmat *
 *         On entry, A points to the data structure containing the local
 *         array information.
 *
 * ---------------------------------------------------------------------
 */ 

#ifdef HPL_CALL_CALDGEMM
extern volatile size_t HPL_CALDGEMM_swap_current_n;
extern int HPL_CALDGEMM_gpu_height;
volatile size_t HPL_CALDGEMM_swap_current_n;

HPL_T_grid* HPL_CALDGEMM_wrapper_grid = NULL;
HPL_T_panel* HPL_CALDGEMM_wrapper_panel = NULL;
HPL_T_panel* HPL_CALDGEMM_wrapper_panel_work = NULL;
int HPL_CALDGEMM_wrapper_icurcol = -1;
int HPL_CALDGEMM_wrapper_n = -1;
#endif

int* permU = NULL;
int dtrtri_(char *, char *, int *, double *, int *, int *);

void HPL_pdgesv_swap_prepare(HPL_T_grid* Grid, HPL_T_panel* panel, int n)
{
	fprintfctd(stderr, "Starting LASWP/DTRSM\n");
	
	if (panel->grid->nprow > 1)
	{
		HPL_ptimer_detail( HPL_TIMING_LASWP );
		permU = HPL_pdlaswp01T( panel, n );
		/*const size_t LDU = panel->grid->nprow == 1 ? panel->lda : (n + (8 - n % 8) % 8 + (((n + (8 - n % 8) % 8) % 16) == 0) * 8);
		double* Uptr = panel->grid->nprow == 1 ? panel->A : panel->U;
		HPL_dlaswp10N( n, panel->jb, Uptr, LDU, permU );*/
		HPL_ptimer_detail( HPL_TIMING_LASWP );
	}
	/*else
	{
		HPL_ptimer_detail( HPL_TIMING_LASWP );
		HPL_dlaswp00N( panel->jb, n, panel->A, panel->lda, panel->IWORK );
		HPL_ptimer_detail( HPL_TIMING_LASWP );
	}*/
}

void HPL_pdgesv_swap(HPL_T_grid* Grid, HPL_T_panel* panel, int n)
{
	int jb = panel->jb;
	size_t lda = panel->lda;
	const size_t LDU = panel->grid->nprow == 1 ? lda : (n + (8 - n % 8) % 8 + (((n + (8 - n % 8) % 8) % 16) == 0) * 8);
	double* Aptr = panel->A;
	double* L1ptr = panel->L1;
	double* Uptr = panel->grid->nprow == 1 ? panel->A : panel->U;
	int* ipiv = panel->IWORK;
	
#ifdef HPL_CALL_CALDGEMM
	size_t HPL_CALDGEMM_wrapper_laswp_stepsize = (HPL_CALDGEMM_gpu_height == 0 ? n : HPL_CALDGEMM_gpu_height);
#else
	const size_t HPL_CALDGEMM_wrapper_laswp_stepsize = n;
#endif

#ifndef CALDGEMM_TEST_DEBUG
	HPL_ptimer_detail( HPL_TIMING_DTRSM );
#endif
	int nremain = n;
	for (size_t i = 0;i < n;i += HPL_CALDGEMM_wrapper_laswp_stepsize)
	{
#ifdef HPL_CALL_CALDGEMM
		if (i) HPL_CALDGEMM_wrapper_laswp_stepsize *= 3;
		const int nn = Mmin(nremain, HPL_CALDGEMM_wrapper_laswp_stepsize);
#else
		const int nn = nremain;
#endif
		nremain -= nn;

#ifdef CALDGEMM_TEST_DEBUG
		HPL_ptimer_detail( HPL_TIMING_LASWP );
#endif
		if (panel->grid->nprow == 1)
		{
			HPL_dlaswp00N( jb, nn, Aptr + i * lda, lda, ipiv );
		}
		else
		{
			if (permU) HPL_dlaswp10N( nn, jb, Uptr + i, LDU, permU );
		}
#ifdef CALDGEMM_TEST_DEBUG
		HPL_ptimer_detail( HPL_TIMING_LASWP );
		HPL_ptimer_detail( HPL_TIMING_DTRSM );
#endif
		if (panel->grid->nprow == 1)
		{
#ifdef HPL_SLOW_CPU
			if (nn > 2 * jb)
			{
				int tmp_lda = jb;
				if (tmp_lda % 8) tmp_lda += (8 - tmp_lda % 8);
				if (((tmp_lda / 8) & 1) == 0) tmp_lda += 8;
				int tmp_ldb = jb;
				if (tmp_ldb % 8) tmp_ldb += (8 - tmp_ldb % 8);
				if (((tmp_ldb / 8) & 1) == 0) tmp_ldb += 8;
				
				double* tmpa = (double*) malloc(sizeof(double) * jb * tmp_lda);
				double* tmpb = (double*) malloc(sizeof(double) * nn * tmp_ldb);
				HPL_dlacpy(jb, jb, L1ptr, jb, tmpa, tmp_lda, 1);
				HPL_dlacpy(jb, nn, Uptr + i * LDU, LDU, tmpb, tmp_ldb, 1);
			
				int ii,jj;
			
				for (ii = 0;ii < jb;ii++)
				{
					tmpa[ii * (tmp_lda + 1)] = 1.;
					for (jj = ii + 1;jj < jb;jj++)
					{
						tmpa[ii * tmp_lda + jj] = 0;
					}
				}
				int tmp;
				dtrtri_("U", "U", &jb, tmpa, &tmp_lda, &tmp);
				HPL_gpu_dgemm(HplColumnMajor, HplTrans, HplNoTrans, jb, nn, jb, 1.0, tmpa, tmp_lda, tmpb, tmp_ldb, 0.0, Uptr + i * LDU, LDU, 0 );
				free(tmpa);
				free(tmpb);
			}
			else
#endif
			{
				HPL_dtrsm( HplColumnMajor, HplLeft, HplUpper, HplTrans, HplUnit, jb, nn, HPL_rone, L1ptr, jb, Uptr + i * LDU, LDU );
			}
			
		}
		else
		{
			HPL_dtrsm( HplColumnMajor, HplRight, HplUpper, HplNoTrans, HplUnit, nn, jb, HPL_rone, L1ptr, jb, Uptr + i, LDU );
		}
#ifdef CALDGEMM_TEST_DEBUG
		HPL_ptimer_detail( HPL_TIMING_DTRSM );
#endif
#ifdef HPL_CALL_CALDGEMM
		HPL_CALDGEMM_swap_current_n = i + nn;
#endif
		
		//fprintf(stderr, "Done at %lld\n", (size_t) i + nn);
	}

#ifndef CALDGEMM_TEST_DEBUG
	HPL_ptimer_detail( HPL_TIMING_DTRSM );
#endif
	
	fprintfctd(stderr, "LASWP/DTRSM finished\n");
}


void HPL_pdgesv_factorize(HPL_T_grid* Grid, HPL_T_panel* panel, int icurcol)
{
	int mycol = Grid->mycol;
	fprintfctd(stderr, "Running Factorize\n");
	if(mycol == icurcol)
	{
		HPL_pdfact(panel);    //factor current panel
	}
	fprintfctd(stderr, "Factorize Ended\n");
}

void HPL_pdgesv_broadcast(HPL_T_grid* Grid, HPL_T_panel* panel, int icurcol)
{
#ifdef HPL_DETAILED_TIMING
   struct timeval tp;
   long start, startu;
   double time, throughput;
#endif

	fprintfctd(stderr, "Starting Broadcast\n");
	HPL_ptimer_detail(HPL_TIMING_BCAST);
	HPL_binit(panel);

#ifdef HPL_DETAILED_TIMING
   if (panel->grid->mycol == panel->pcol)
   {
	(void) gettimeofday( &tp, NULL );
	start  = tp.tv_sec;
	startu = tp.tv_usec;
   }
#endif

	HPL_bcast(panel);

#ifdef HPL_DETAILED_TIMING
   if (panel->grid->mycol == panel->pcol)
   {
	(void) gettimeofday( &tp, NULL );
	time = (double)( tp.tv_sec - start ) + ( (double)( tp.tv_usec-startu ) / 1000000.0 );
	throughput = (double) panel->len * sizeof(double) / time / 1000000.;
	fprintf(stderr, "MPI Broadcast: size=%f MB - throughput=%f MB/s\n", (double) panel->len * (double) sizeof(double) / 1024. / 1024., throughput);
   }
#endif

	HPL_ptimer_detail(HPL_TIMING_BCAST);
	fprintfctd(stderr, "Broadcast Ended\n");
}

#ifdef HPL_CALL_CALDGEMM

void HPL_CALDGEMM_wrapper_factorize()
{
	HPL_pdgesv_factorize(HPL_CALDGEMM_wrapper_grid, HPL_CALDGEMM_wrapper_panel, HPL_CALDGEMM_wrapper_icurcol);
}
void HPL_CALDGEMM_wrapper_broadcast()
{
	HPL_pdgesv_broadcast(HPL_CALDGEMM_wrapper_grid, HPL_CALDGEMM_wrapper_panel, HPL_CALDGEMM_wrapper_icurcol);
}
void HPL_CALDGEMM_wrapper_swap()
{
	HPL_pdgesv_swap(HPL_CALDGEMM_wrapper_grid, HPL_CALDGEMM_wrapper_panel_work, HPL_CALDGEMM_wrapper_n);
}

#endif

void HPL_pdupdateTT(HPL_T_grid* Grid, HPL_T_panel* PBCST, HPL_T_panel* PANEL, const int NN, int factorize, int depth2)
{
	//.. Local Variables ..
	double * Aptr, * L1ptr, * L2ptr, * Uptr, * dpiv;
	int * ipiv;
	int curr, i, iroff, jb, lda, ldl2, mp, n, nb;
#ifdef HPL_GPU_TEMPERATURE_THRESHOLD
	double temperature;
#endif
	//.. Executable Statements ..
	fprintfctd(stderr, "Running pdupdateTT\n");
	HPL_ptimer_detail( HPL_TIMING_UPDATE );
	nb = PANEL->nb;
	jb = PANEL->jb;
	n = PANEL->nq;
	lda = PANEL->lda;
	if( NN >= 0 ) n = Mmin( NN, n );

	const int LDU = PANEL->grid->nprow == 1 ? lda : (n + (8 - n % 8) % 8 + (((n + (8 - n % 8) % 8) % 16) == 0) * 8);

	Aptr = PANEL->A;
	L2ptr = PANEL->L2;
	L1ptr = PANEL->L1;
	ldl2 = PANEL->ldl2;
	curr = ( PANEL->grid->myrow == PANEL->prow ? 1 : 0 );
	Uptr = PANEL->grid->nprow == 1 ? PANEL->A : PANEL->U;
	dpiv = PANEL->DPIV;
	ipiv = PANEL->IWORK;
	iroff = PANEL->ii;
	mp = PANEL->mp - ( PANEL->grid->nprow == 1 || curr != 0 ? jb : 0 );

	if( PANEL->grid->nprow == 1 ) for( i = 0; i < jb; i++ ) { ipiv[i] = (int)(dpiv[i]) - iroff; }

	if (n)
	{
#ifdef HPL_CALL_CALDGEMM
		HPL_CALDGEMM_wrapper_n = n;
		HPL_CALDGEMM_wrapper_grid = Grid;
		HPL_CALDGEMM_wrapper_panel = PBCST;
		HPL_CALDGEMM_wrapper_panel_work = PANEL;
		HPL_CALDGEMM_wrapper_icurcol = factorize;
#endif
		HPL_pdgesv_swap_prepare(Grid, PANEL, n);
#ifdef HPL_CALL_CALDGEMM
		if (depth2
#ifdef HPL_LOOKAHEAD2_TURNOFF
		 && n >= HPL_LOOKAHEAD2_TURNOFF
 #endif
		 )
		{
		    CALDGEMM_enable_async_laswp(1);
		}
		else
		{
		    HPL_CALDGEMM_gpu_height = 0;
		    CALDGEMM_enable_async_laswp(0);
#endif
		    HPL_pdgesv_swap(Grid, PANEL, n);
#ifdef HPL_CALL_CALDGEMM
		}
#endif
	
		HPL_ptimer_detail( HPL_TIMING_DGEMM );
#ifdef HPL_CALL_CALDGEMM
		int caldgemm_linpack_mode = (factorize != -1) ? (Grid->mycol == HPL_CALDGEMM_wrapper_icurcol ? 2 : 1) : 0;
		//caldgemm_linpack_mode = 0;
		HPL_gpu_dgemm( HplColumnMajor, HplNoTrans, PANEL->grid->nprow == 1 ? HplNoTrans : HplTrans, mp, n, jb, -HPL_rone, L2ptr, ldl2, Uptr, LDU, HPL_rone, (PANEL->grid->nprow == 1 || curr != 0) ? Mptr( Aptr, jb, 0, lda ) : Aptr, lda, caldgemm_linpack_mode );
#ifdef HPL_GPU_TEMPERATURE_THRESHOLD
		if (adl_temperature_check_run(&temperature, 1))
		{
			fprintf(stderr, "Error running temperature check\n");
			exit(1);
		}
		if (temperature >= HPL_GPU_TEMPERATURE_THRESHOLD)
		{
			fprintf(stderr, "GPU Temperature Threshold of %f exceeded, GPU temperature is %f\n", (double) HPL_GPU_TEMPERATURE_THRESHOLD, temperature);
			exit(1);
		}
#endif
#else
		HPL_dgemm( HplColumnMajor, HplNoTrans, PANEL->grid->nprow == 1 ? HplNoTrans : HplTrans, mp, n, jb, -HPL_rone, L2ptr, ldl2, Uptr, LDU, HPL_rone, (PANEL->grid->nprow == 1 || curr != 0) ? Mptr( Aptr, jb, 0, lda ) : Aptr, lda );
#endif
		HPL_ptimer_detail( HPL_TIMING_DGEMM );
		
#ifndef HPL_CALL_CALDGEMM
		if (factorize != -1)
		{
			HPL_pdgesv_factorize(Grid, PBCST, factorize);
			HPL_pdgesv_broadcast(Grid, PBCST, factorize);
		}
#endif

		if (PANEL->grid->nprow != 1 && curr != 0)
		{
			HPL_ptimer_detail( HPL_TIMING_DLATCPY );
			HPL_dlatcpy( jb, n, Uptr, LDU, Aptr, lda );
			HPL_ptimer_detail( HPL_TIMING_DLATCPY );
		}
	}
	else if (factorize != -1)
	{
		HPL_pdgesv_factorize(Grid, PBCST, factorize);
		HPL_pdgesv_broadcast(Grid, PBCST, factorize);
	}

	HPL_ptimer_detail( HPL_TIMING_UPDATE );

	fprintfctd(stderr, "pdupdateTT ended\n");
}

void PrintMatrix(HPL_T_grid* GRID, HPL_T_pmat* A)
{
	for (int i = 0;i < A->n;i++)
	{
		for (int j = 0;j < A->n + 1;j++)
		{
			int tmp[6];
			HPL_infog2l(i, j, A->nb, A->nb, 0, 0, GRID->myrow, GRID->mycol, GRID->nprow, GRID->npcol, tmp, tmp+1, tmp+2, tmp+3, GRID);
			if (GRID->myrow == 0 && GRID->mycol == 0)
			{
				if (tmp[2] || tmp[3])
				{
					int rank;
					if (GRID->order == HPL_ROW_MAJOR)
					{
						rank = tmp[2] * GRID->npcol + tmp[3];
					}
					else
					{
						rank = tmp[3] * GRID->nprow + tmp[2];
					}

					MPI_Recv(tmp, 6, MPI_INT, rank, i * A->n + j, GRID->all_comm, MPI_STATUS_IGNORE);
					//fprintf(stderr, "Guessing Row %d Col %d Rank %d\n", tmp[2], tmp[3], rank);
				}
				else
				{
					*((double*) (tmp+4)) = A->A[tmp[1] * A->ld + tmp[0]];
				}
			}
			else if (tmp[2] == GRID->myrow && tmp[3] == GRID->mycol)
			{
				//fprintf(stderr, "Row %d Col %d Rank %d\n", tmp[2], tmp[3], GRID->iam);
				*((double*) (tmp+4)) = A->A[tmp[1] * A->ld + tmp[0]];
				MPI_Send(tmp, 6, MPI_INT, 0, i * A->n + j, GRID->all_comm);
			}
			if (GRID->iam == 0)
			{
				//fprintf(stderr, "i %d j %d Row %d Col %d II %d JJ %d entry %lf\n", i, j, tmp[2], tmp[3], tmp[0], tmp[1], *((double*) (tmp+4)));
				fprintf(stderr, "%lf\t", *((double*) (tmp+4)));
			}
		}
		if (GRID->iam == 0) fprintf(stderr, "\n");
	}
	if (GRID->iam == 0) fprintf(stderr, "\n");
}

void PrintVector(HPL_T_grid* GRID, HPL_T_pmat* A)
{
	//Currently only works with P=1 (probably)
	if (GRID->iam == 0) fprintf(stderr, "\n");
	double* buffer = malloc(A->n * sizeof(double));
	for (int i = 0;i < GRID->npcol;i++)
	{
		if (i == 0)
		{
			for (int j = 0;j < A->n;j++) buffer[j] = A->X[j];
		}
		else
		{
			if (GRID->iam == 0)
			{
				MPI_Recv(buffer, A->n, MPI_DOUBLE, i, 0, GRID->all_comm, MPI_STATUS_IGNORE);
			}
			else if (GRID->mycol == i)
			{
				MPI_Send(A->X, A->n, MPI_DOUBLE, 0, 0, GRID->all_comm);
			}
		}
		if (GRID->iam == 0) fprintf(stderr, "Rank %d:\t", i);
		for (int j = 0;j < A->n;j++)
		{
			if (GRID->iam == 0) fprintf(stderr, "%lf\t", buffer[j]);
		}
		if (GRID->iam == 0) fprintf(stderr, "\n");
	}
	free(buffer);
	if (GRID->iam == 0) fprintf(stderr, "\n");
}

void HPL_pdgesv(HPL_T_grid* GRID, HPL_T_palg* ALGO, HPL_T_pmat* A)
{
	//.. Local Variables ..
	HPL_T_panel *p, **panel = NULL;
	int N, depth1, depth2, icurcol, j, jb, mycol, n, nb, nn, npcol, nq, tag=MSGID_BEGIN_FACT;
#ifdef HPL_PRINT_INTERMEDIATE
	uint64_t total_gflop;
	uint64_t time_start;
#endif
	//.. Executable Statements ..
	if( A->n <= 0 ) return;
	A->info = 0;
	
	mycol = GRID->mycol;
	npcol = GRID->npcol;
	N = A->n;
	nb = A->nb;
	depth1 = (ALGO->depth >= 1);
	depth2 = (ALGO->depth >= 2);
	

#ifdef HPL_PRINT_INTERMEDIATE
	total_gflop = 2 * (uint64_t) N * N * N / 3 / 1e9;
	time_start = util_getTimestamp();
#endif /* HPL_PRINT_INTERMEDIATE */

	//Allocate the panel list
	panel = (HPL_T_panel **)malloc((size_t)(depth1+1) * sizeof(HPL_T_panel *));
	if(panel == NULL) HPL_pabort(__LINE__, "HPL_pdgesvK2", "Memory allocation failed");

	//Create and initialize the lookahead panel
	nq = HPL_numcol(N+1, nb, mycol, GRID);
	nn = N;

	if (depth1)
	{
		HPL_pdpanel_new(GRID, ALGO, nn, nn+1, Mmin(nn, nb), A, 0, 0, tag, &panel[0]);
	}

	//Create main panel
	HPL_pdpanel_new(GRID, ALGO, nn, nn+1, Mmin(nn, nb), A, 0, 0, tag, &panel[depth1]);
	tag = MNxtMgid(tag, MSGID_BEGIN_FACT, MSGID_END_FACT);
	
#ifdef HPL_START_PERCENTAGE
	double fullwork = N;
	fullwork = fullwork * fullwork * fullwork;
	fullwork *= 1.0 - (double) HPL_START_PERCENTAGE / 100.0;
	fullwork = pow(fullwork, 1.0 / 3.0);
	int startrow = N - fullwork;
	if (startrow < 0) startrow = 0;
	startrow -= startrow % (npcol * nb);
	if( GRID->myrow == 0 && GRID->mycol == 0 )
	{
	    fprintf(stderr, "Starting at col %d which corresponds to approx %2.1lf %% of execution time\n", startrow, 100.0 * (double) (N - startrow) * (double) (N - startrow) * (double) (N - startrow) / (double) N / (double) N / (double) N);
	}
#else
	const int startrow = 0;
#endif
	
	//Main loop over the columns of A
	for(j = startrow; j < N; j += nb)
	{
		icurcol = MColToPCol(j, nb, npcol, GRID);
		n = N - j;
#ifdef HPL_HALF_BLOCKING
		if (n <= HPL_HALF_BLOCKING)// && n + nb > HPL_HALF_BLOCKING && npcol == 1 && GRID->nprow == 1 && depth1 == 0)
		{
			nb = A->nb / 2;
		}
#endif
#ifdef HPL_CUSTOM_PARAMETER_CHANGE
		HPL_CUSTOM_PARAMETER_CHANGE
#endif
		jb = Mmin(n, nb);
#ifdef HPL_DETAILED_TIMING
		fprintfct(stderr, "Iteration j=%d N=%d n=%d jb=%d Totaltime=%2.3lf\n", j, N, n, jb, HPL_ptimer_inquire( HPL_WALL_PTIME, HPL_TIMING_ITERATION ));
#else
		fprintfct(stderr, "Iteration j=%d N=%d n=%d jb=%d\n", j, N, n, jb);
#endif

#if defined(HPL_PRINT_INTERMEDIATE) & !defined(HPL_PAUSE)
		// there are still n rows to compute
		// j rows have already been computed
		if( GRID->myrow == 0 && GRID->mycol == 0 )
		{
			uint64_t todo_gflop = 2 * (uint64_t) n * n * n / 3 / 1e9;
			uint64_t gFlop = total_gflop - todo_gflop;
			float ratio = (float) gFlop / total_gflop;

			float modifier = ratio * ratio * ratio * ratio;
			modifier = 1. - modifier;
			modifier = 1. + 0.1 * modifier;

			uint64_t time_now = util_getTimestamp();
			float seconds = (float) util_getTimeDifference( time_start, time_now ) / 1e6;
			if( seconds != 0 && j != 0 )
			{
				float flops = (float) gFlop / (float) seconds / modifier;
				uint64_t eta = (seconds / ratio - seconds) * modifier;
				uint64_t wall_now = util_getWalltime() / 1e3;
				//printf( "%f %% of factorization (%.2f GFlop) done in %ld s at approx. %.2f Gflops\n", ratio * 100, (float) gFlop, seconds, flops );
				printf( "[%ld] %.f %% (j = %d/%d) of factorization at approx. %.2f Gflops, assuming to finish in %ld s.\n", wall_now, ratio * 100, j, N, flops, eta );
			}
		}
#endif /* HPL_PRINT_INTERMEDIATE */

		//Initialize current panel
		HPL_ptimer_detail( HPL_TIMING_ITERATION );

		if (j == startrow || depth1 == 0)
		{
			HPL_pdpanel_free(panel[depth1]);
			HPL_pdpanel_init(GRID, ALGO, n, n + 1, jb, A, j, j, tag, panel[depth1]);
			HPL_pdgesv_factorize(GRID, panel[depth1], icurcol);
			HPL_pdgesv_broadcast(GRID, panel[depth1], icurcol);
		}
		
		tag = MNxtMgid(tag, MSGID_BEGIN_FACT, MSGID_END_FACT);
		
		if (depth1 && j + nb < N)
		{
			HPL_pdpanel_free(panel[0]);
			HPL_pdpanel_init(GRID, ALGO, n - nb, n - nb + 1, Mmin(n - nb, nb), A, j + nb, j + nb, tag, panel[0]);
		}
		
		nn = (mycol == icurcol) ? HPL_numcolI(jb, j, nb, mycol, GRID) : 0;

		//Finish the latest update and broadcast the current panel

		int olddepth1 = depth1;
#ifdef HPL_DISABLE_LOOKAHEAD
		if (n <= HPL_DISABLE_LOOKAHEAD + nb + 1)
		{
			depth1 = 0;
		}
#endif
		HPL_pdupdateTT(GRID, panel[0], panel[olddepth1], nq-nn, (depth1 && j + nb < N) ? MColToPCol(j + nb, nb, npcol, GRID) : -1, depth2);

		HPL_ptimer_detail( HPL_TIMING_ITERATION );
		//Switch panel pointers
		if (depth1)
		{
			p = panel[0];
			panel[0] = panel[1];
			panel[1] = p;
		}

		if(mycol == icurcol)
		{
			nq -= jb;
		}

#ifdef HPL_PAUSE
		HPL_ptimer( 0 );
		printf("HPL_PAUSE: Sleeping for %d usec\n", HPL_PAUSE * n);
		usleep(HPL_PAUSE * n);
		HPL_ptimer( 0 );
#endif
	}
	//Clean-up: Release panels and panel list
	if(depth1)
	{
		HPL_pdpanel_disp( &panel[0]);
	}
	HPL_pdpanel_disp(&panel[depth1]);
	if(panel) free(panel);
	
	//Solve upper triangular system
	if( A->info == 0 ) HPL_pdtrsv( GRID, A );
}
