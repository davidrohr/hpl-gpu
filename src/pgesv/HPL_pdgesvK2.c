/* 
 * This is a modified version of the High Performance Computing Linpack
 * Benchmark (HPL). All code not contained in the original HPL version
 * 2.0 is property of the Frankfurt Institute for Advanced Studies
 * (FIAS). None of the material may be copied, reproduced, distributed,
 * republished, downloaded, displayed, posted or transmitted in any form
 * or by any means, including, but not limited to, electronic,
 * mechanical, photocopying, recording, or otherwise, without the prior
 * written permission of FIAS. For those parts contained in the
 * unmodified version of the HPL the below copyright notice applies.
 * 
 * Authors:
 * David Rohr (drohr@jwdt.org)
 * Matthias Bach (bach@compeng.uni-frankfurt.de)
 * Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 * 
 * -- High Performance Computing Linpack Benchmark (HPL)                
 *    HPL - 2.0 - September 10, 2008                          
 *    Antoine P. Petitet                                                
 *    University of Tennessee, Knoxville                                
 *    Innovative Computing Laboratory                                 
 *    (C) Copyright 2000-2008 All Rights Reserved                       
 *                                                                      
 * -- Copyright notice and Licensing terms:                             
 *                                                                      
 * Redistribution  and  use in  source and binary forms, with or without
 * modification, are  permitted provided  that the following  conditions
 * are met:                                                             
 *                                                                      
 * 1. Redistributions  of  source  code  must retain the above copyright
 * notice, this list of conditions and the following disclaimer.        
 *                                                                      
 * 2. Redistributions in binary form must reproduce  the above copyright
 * notice, this list of conditions,  and the following disclaimer in the
 * documentation and/or other materials provided with the distribution. 
 *                                                                      
 * 3. All  advertising  materials  mentioning  features  or  use of this
 * software must display the following acknowledgement:                 
 * This  product  includes  software  developed  at  the  University  of
 * Tennessee, Knoxville, Innovative Computing Laboratory.             
 *                                                                      
 * 4. The name of the  University,  the name of the  Laboratory,  or the
 * names  of  its  contributors  may  not  be used to endorse or promote
 * products  derived   from   this  software  without  specific  written
 * permission.                                                          
 *                                                                      
 * -- Disclaimer:                                                       
 *                                                                      
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR  CONTRIBUTORS  BE  LIABLE FOR ANY  DIRECT,  INDIRECT,  INCIDENTAL,
 * SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES  (INCLUDING,  BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,  OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 * ---------------------------------------------------------------------
 */ 
/*
 * Include files
 */
#include "hpl.h"

void HPL_pdgesvK2(HPL_T_grid* GRID, HPL_T_palg* ALGO, HPL_T_pmat* A)
{
/* 
 * Purpose
 * =======
 *
 * HPL_pdgesvK2 factors a N+1-by-N matrix using LU factorization with row
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
	//.. Local Variables ..
	HPL_T_panel                * p, * * panel = NULL;
	int                        N, depth, icurcol=0, j, jb, jj=0, jstart,
		k, mycol, n, nb, nn, npcol, nq,
		tag=MSGID_BEGIN_FACT, test=HPL_KEEP_TESTING;
	//.. Executable Statements ..
	mycol = GRID->mycol;
	npcol = GRID->npcol;
	depth = ALGO->depth;
	N = A->n;
	nb = A->nb;

	if( N <= 0 ) return;

	//Allocate a panel list of length depth + 1 (depth >= 1)
	panel = (HPL_T_panel **)malloc( (size_t)(depth+1) * sizeof( HPL_T_panel *) );
	if( panel == NULL ) HPL_pabort( __LINE__, "HPL_pdgesvK2", "Memory allocation failed" );

	//Create and initialize the lookahead panel
	nq = HPL_numroc( N+1, nb, nb, mycol, 0, npcol );
	nn = N;
	jstart = 0;

	if (depth)
	{
		HPL_pdpanel_new( GRID, ALGO, nn, nn+1, Mmin( nn, nb ), A, jstart, jstart, tag, &panel[0] );
	}

	//Create main panel
	HPL_pdpanel_new( GRID, ALGO, nn, nn+1, Mmin( nn, nb ), A, jstart, jstart, tag, &panel[depth] );
	tag = MNxtMgid( tag, MSGID_BEGIN_FACT, MSGID_END_FACT );

	//Main loop over the remaining columns of A
	for( j = jstart; j < N; j += nb )
	{
		n = N - j;
		jb = Mmin( n, nb );

		//Initialize current panel
		(void) HPL_pdpanel_free( panel[depth] );
		HPL_pdpanel_init( GRID, ALGO, n, n+1, jb, A, j, j, tag, panel[depth] );

		if( mycol == icurcol )
		{
			nn = HPL_numrocI( jb, j, nb, nb, mycol, 0, npcol );
			HPL_pdfact( panel[depth] );    //factor current panel
		}
		else
		{
			nn = 0;
		}

		HPL_ptimer_detail( HPL_TIMING_BCAST );
		(void) HPL_binit(         panel[depth] );
		do
		{
			(void) HPL_bcast(       panel[depth], &test );
		}
		while( test != HPL_SUCCESS );
		(void) HPL_bwait(         panel[depth] );
		HPL_ptimer_detail( HPL_TIMING_BCAST );

		//Finish the latest update and broadcast the current panel
		HPL_pdupdateTT( NULL, NULL, panel[depth], nq-nn );

		//Circular of the panel pointers: * xtmp = x[0]; for( k=0; k < depth; k++ ) x[k] = x[k+1]; x[d] = xtmp;
		//Go to next process row and column - update the message ids for broadcast
		if (depth)
		{
			p = panel[0];
			panel[0] = panel[1];
			panel[1] = p;
		}

		if( mycol == icurcol ) { jj += jb; nq -= jb; }
		icurcol = MModAdd1( icurcol, npcol );
		tag = MNxtMgid( tag, MSGID_BEGIN_FACT, MSGID_END_FACT );
	}

	//Clean-up: Finish updates - release panels and panel list
	nn = HPL_numrocI( 1, N, nb, nb, mycol, 0, npcol );
	if(depth)
	{
		(void) HPL_pdpanel_disp(  &panel[0] );
	}
	(void) HPL_pdpanel_disp( &panel[depth] );

	if( panel ) free( panel );
}
