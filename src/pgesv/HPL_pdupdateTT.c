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
 * HPL - 2.0 - September 10, 2008 
 * Antoine P. Petitet 
 * University of Tennessee, Knoxville 
 * Innovative Computing Laboratory 
 * (C) Copyright 2000-2008 All Rights Reserved 
 * 
 * -- Copyright notice and Licensing terms: 
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met: 
 * 
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer. 
 * 
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions, and the following disclaimer in the
 * documentation and/or other materials provided with the distribution. 
 * 
 * 3. All advertising materials mentioning features or use of this
 * software must display the following acknowledgement: 
 * This product includes software developed at the University of
 * Tennessee, Knoxville, Innovative Computing Laboratory. 
 * 
 * 4. The name of the University, the name of the Laboratory, or the
 * names of its contributors may not be used to endorse or promote
 * products derived from this software without specific written
 * permission. 
 * 
 * -- Disclaimer: 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 * ---------------------------------------------------------------------
 */ 
/*
 * Include files
 */
#include "hpl.h"

void HPL_pdupdateTT(HPL_T_panel* PBCST, int* IFLAG, HPL_T_panel* PANEL, const int NN)
{
/* 
 * Purpose
 * =======
 *
 * HPL_pdupdateTT broadcast - forward the panel PBCST and simultaneously
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
 * IFLAG (local output) int *
 * On exit, IFLAG indicates whether or not the broadcast has
 * been completed when PBCST is not NULL on entry. In that case,
 * IFLAG is left unchanged.
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
	//.. Local Variables ..
	double * Aptr, * L1ptr, * L2ptr, * Uptr, * dpiv;
	int * ipiv;
	int curr, i, iroff, jb, lda, ldl2, mp, n, nb, test;
	//.. Executable Statements ..
	fprintfct(stderr, "Running pdupdateTT\n");
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
		HPL_ptimer_detail( HPL_TIMING_LASWP );
		if (PANEL->grid->nprow == 1) HPL_dlaswp00N( jb, n, Aptr, lda, ipiv );
		else HPL_pdlaswp01T( PBCST, &test, PANEL, n );
		HPL_ptimer_detail( HPL_TIMING_LASWP );

		HPL_ptimer_detail( HPL_TIMING_DTRSM );
		if (PANEL->grid->nprow == 1) HPL_dtrsm( HplColumnMajor, HplLeft, HplUpper, HplTrans,	HplUnit, jb, n, HPL_rone, L1ptr, jb, Uptr, LDU );
		else                         HPL_dtrsm( HplColumnMajor, HplRight, HplUpper, HplNoTrans, HplUnit, n, jb, HPL_rone, L1ptr, jb, Uptr, LDU );
		HPL_ptimer_detail( HPL_TIMING_DTRSM );

		HPL_ptimer_detail( HPL_TIMING_DGEMM );
		HPL_dgemm( HplColumnMajor, HplNoTrans, PANEL->grid->nprow == 1 ? HplNoTrans : HplTrans, mp, n, jb, -HPL_rone, L2ptr, ldl2, Uptr, LDU, HPL_rone, (PANEL->grid->nprow == 1 || curr != 0) ? Mptr( Aptr, jb, 0, lda ) : Aptr, lda );
		HPL_ptimer_detail( HPL_TIMING_DGEMM );

		if (PANEL->grid->nprow != 1 && curr != 0)
		{
			HPL_ptimer_detail( HPL_TIMING_DLATCPY );
			HPL_dlatcpy( jb, n, Uptr, LDU, Aptr, lda );
			HPL_ptimer_detail( HPL_TIMING_DLATCPY );
		}
	}

	HPL_ptimer_detail( HPL_TIMING_UPDATE );

	fprintfct(stderr, "pdupdateTT ended\n");
}
