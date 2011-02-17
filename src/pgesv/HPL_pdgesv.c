/*
 *  -- High Performance Computing Linpack Benchmark (HPL-GPU)
 *     HPL-GPU - 1.0 - 2010
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
volatile size_t HPL_CALDGEMM_swap_current_n;

HPL_T_grid* HPL_CALDGEMM_wrapper_grid = NULL;
HPL_T_panel* HPL_CALDGEMM_wrapper_panel = NULL;
HPL_T_panel* HPL_CALDGEMM_wrapper_panel_work = NULL;
int HPL_CALDGEMM_wrapper_icurcol = -1;
int HPL_CALDGEMM_wrapper_n = -1;
size_t HPL_CALDGEMM_wrapper_laswp_stepsize;
#endif

int* permU = NULL;

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
	
#ifndef HPL_CALL_CALDGEMM
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
			HPL_dtrsm( HplColumnMajor, HplLeft, HplUpper, HplTrans,	HplUnit, jb, nn, HPL_rone, L1ptr, jb, Uptr + i * LDU, LDU );
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
	int test = HPL_KEEP_TESTING;
	fprintfctd(stderr, "Starting Broadcast\n");
	HPL_ptimer_detail(HPL_TIMING_BCAST);
	HPL_binit(panel);
	do
	{
		HPL_bcast(panel, &test);
	}
	while(test != HPL_SUCCESS);
	HPL_bwait(panel);
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
		if (depth2 && n >= 56 * 1024)
		{
		    HPL_CALDGEMM_wrapper_laswp_stepsize = 5120;
		    CALDGEMM_enable_async_laswp(1);
		}
		else
		{
		    HPL_CALDGEMM_wrapper_laswp_stepsize = n;
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
	nq = HPL_numcol(N+1, nb, nb, mycol, npcol, GRID);
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
		icurcol = MColToPCol(j, nb, npcol);
		n = N - j;
		jb = Mmin(n, nb);
#ifdef HPL_DETAILED_TIMING
		fprintfct(stderr, "Iteration j=%d N=%d n=%d jb=%d Totaltime=%2.3lf\n", j, N, n, jb, HPL_ptimer_inquire( HPL_WALL_PTIME, HPL_TIMING_ITERATION ));
#else
		fprintfct(stderr, "Iteration j=%d N=%d n=%d jb=%d\n", j, N, n, jb);
#endif

#ifdef HPL_PRINT_INTERMEDIATE
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
			uint64_t seconds = util_getTimeDifference( time_start, time_now ) / 1e6;
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
		

		nn = (mycol == icurcol) ? HPL_numcolI(jb, j, nb, nb, mycol, npcol, GRID) : 0;

		//Finish the latest update and broadcast the current panel
		HPL_pdupdateTT(GRID, panel[0], panel[depth1], nq-nn, (depth1 && j + nb < N) ? MColToPCol(j + nb, nb, npcol) : -1, depth2);

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
