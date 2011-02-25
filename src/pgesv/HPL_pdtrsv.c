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

void HPL_pdtrsv(HPL_T_grid* GRID, HPL_T_pmat* AMAT)
{
/* 
 * Purpose
 * =======
 *
 * HPL_pdtrsv solves an upper triangular system of linear equations.
 *  
 * The rhs is the last column of the N by N+1 matrix A. The solve starts
 * in the process  column owning the  Nth  column of A, so the rhs b may
 * need to be moved one process column to the left at the beginning. The
 * routine therefore needs  a column  vector in every process column but
 * the one owning  b. The result is  replicated in all process rows, and
 * returned in XR, i.e. XR is of size nq = LOCq( N ) in all processes.
 *  
 * The algorithm uses decreasing one-ring broadcast in process rows  and
 * columns  implemented  in terms of  synchronous communication point to
 * point primitives.  The  lookahead of depth 1 is used to minimize  the
 * critical path. This entire operation is essentially ``latency'' bound
 * and an estimate of its running time is given by:
 *  
 *    (move rhs) lat + N / ( P bdwth ) +            
 *    (solve)    ((N / NB)-1) 2 (lat + NB / bdwth) +
 *               gam2 N^2 / ( P Q ),                
 *  
 * where  gam2   is an estimate of the   Level 2 BLAS rate of execution.
 * There are  N / NB  diagonal blocks. One must exchange  2  messages of
 * length NB to compute the next  NB  entries of the vector solution, as
 * well as performing a total of N^2 floating point operations.
 *
 * Arguments
 * =========
 *
 * GRID    (local input)                 HPL_T_grid *
 *         On entry,  GRID  points  to the data structure containing the
 *         process grid information.
 *
 * AMAT    (local input/output)          HPL_T_pmat *
 *         On entry,  AMAT  points  to the data structure containing the
 *         local array information.
 *
 * ---------------------------------------------------------------------
 */ 
//Local Variables
	MPI_Comm Ccomm, Rcomm;
	double *A=NULL, *Aprev=NULL, *Aptr, *XC=NULL, *XR=NULL, *Xd=NULL, *Xdprev=NULL, *W=NULL;
	int Alcol_matrix, Alcol_process, Alrow, Anpprev, Anp, Anq, Bcol, Cmsgid, GridIsNotPx1, GridIsNot1xQ, Rmsgid,
		Wfr=0, colprev, kb, kbprev, lda, mycol,	myrow, n, n1, n1p, n1pprev=0, nb, npcol, nprow, rowprev, tmp1, tmp2;

//Executable Statements
	HPL_ptimer_detail( HPL_TIMING_PTRSV );
	if ((n = AMAT->n) <= 0) return;
	nb = AMAT->nb;
	lda = AMAT->ld;
	A = AMAT->A;
	XR = AMAT->X;

	(void) HPL_grid_info(GRID, &nprow, &npcol, &myrow, &mycol);
	//if (mycol >= 2) return;
	//npcol = 2;
	Rcomm = GRID->row_comm;
	Rmsgid = MSGID_BEGIN_PTRSV;
	Ccomm = GRID->col_comm;
	Cmsgid = MSGID_BEGIN_PTRSV + 1;
	GridIsNot1xQ = (nprow > 1);
	GridIsNotPx1 = (npcol > 1);

//Move the rhs in the process column owning the last column of A.
	Mnumrow(Anp, n, nb, myrow, nprow);
	Mnumcol(Anq, n, nb, mycol, GRID);
	fprintfqt(stderr, "Rank %d Anq %d Anp %d\n", GRID->iam, Anq, Anp);

	tmp1 = (n - 1) / nb;
	Alrow = tmp1 % nprow;
	Alcol_matrix = tmp1;
	Alcol_process = MColBlockToPCol(Alcol_matrix, 0, GRID);
	kb = n - tmp1 * nb;

	Aptr = (double *) (A);
	XC = Mptr(Aptr, 0, Anq, lda);
	Mindxg2p_col(n, nb, nb, Bcol, 0, GRID);

	if ((Anp > 0) && (Alcol_process != Bcol))
	{
		if(mycol == Bcol)
		{
			(void) HPL_send(XC, Anp, Alcol_process, Rmsgid, Rcomm);
			fprintfqt(stderr, "Send 1: %d to %d\n", GRID->iam, Alcol_process);
		}
		else if(mycol == Alcol_process)
		{
			(void) HPL_recv(XC, Anp, Bcol, Rmsgid, Rcomm);
			fprintfqt(stderr, "Recv 1: %d from %d\n", GRID->iam, Bcol);
		}
	}
	Rmsgid = (Rmsgid + 2 > MSGID_END_PTRSV ? MSGID_BEGIN_PTRSV : Rmsgid + 2);
	if (mycol != Alcol_process)
	{
		for(tmp1 = 0; tmp1 < Anp; tmp1++)
		{
			XC[tmp1] = HPL_rzero;
		}
	}

//Set up lookahead
	n1 = (npcol - 1) * nb;
	n1 = Mmax(n1, nb);
	if (Anp > 0)
	{
		W = (double*) malloc((size_t) (Mmin(n1, Anp)) * sizeof(double));
		if (W == NULL)
		{
			HPL_pabort(__LINE__, "HPL_pdtrsv", "Memory allocation failed");
		}
		Wfr = 1;
	}

	Anpprev = Anp;
	Xdprev = XR;
	Aprev = Aptr = Mptr(Aptr, 0, Anq, lda);
	tmp1 = n - kb;
	tmp1 -= (tmp2 = Mmin(tmp1, n1));
	MnumrowI(n1pprev, tmp2, Mmax(0, tmp1), nb, myrow, nprow);

	if (myrow == Alrow)
	{
		Anpprev = (Anp -= kb);
	}
	if (mycol == Alcol_process)
	{
		Aprev = (Aptr -= lda * kb);
		Anq -= kb;
		Xdprev = (Xd = XR + Anq);
		if (myrow == Alrow)
		{
			HPL_dtrsv(HplColumnMajor, HplUpper, HplNoTrans, HplNonUnit, kb, Aptr+Anp, lda, XC+Anp, 1);
			HPL_dcopy(kb, XC+Anp, 1, Xd, 1);
		}
	}

	rowprev = Alrow;
	Alrow = MModSub1(Alrow, nprow);
	colprev = Alcol_process;
	Alcol_matrix--;
	Alcol_process = MColBlockToPCol(Alcol_matrix, 0, GRID);
	kbprev  = kb;
	n -= kb;
	tmp1 = n - (kb = nb);
	tmp1 -= (tmp2 = Mmin(tmp1, n1));
	MnumrowI(n1p, tmp2, Mmax(0, tmp1), nb, myrow, nprow);

// Start the operations
	while(n > 0)
	{
		if(mycol == Alcol_process)
		{
			Aptr -= lda * kb;
			Anq -= kb;
			Xd = XR + Anq;
		}
		if(myrow == Alrow)
		{
			Anp -= kb;
		}
/*
 * Broadcast  (decreasing-ring)  of  previous solution block in previous
 * process column,  compute  partial update of current block and send it
 * to current process column.
 */
		if (mycol == colprev)
		{
			//Send previous solution block in process row above
			if (myrow == rowprev)
			{
				if (GridIsNot1xQ)
				{
					(void) HPL_send(Xdprev, kbprev, MModSub1(myrow, nprow), Cmsgid, Ccomm);
					fprintfqt(stderr, "Send 2: %d to %d\n", GRID->iam, MModSub1(myrow, nprow));
				}
			}
			else
			{
				(void) HPL_recv(Xdprev, kbprev, MModAdd1(myrow, nprow), Cmsgid, Ccomm);
				fprintfqt(stderr, "Recv 2: %d from %d\n", GRID->iam, MModAdd1(myrow, nprow));
			} 

			//Compute partial update of previous solution block and send it to current column
			if(n1pprev > 0)
			{
				tmp1 = Anpprev - n1pprev;
				HPL_dgemv(HplColumnMajor, HplNoTrans, n1pprev, kbprev, -HPL_rone, Aprev+tmp1, lda, Xdprev, 1, HPL_rone, XC+tmp1, 1 );
				if(GridIsNotPx1)
				{
					(void) HPL_send(XC+tmp1, n1pprev, Alcol_process, Rmsgid, Rcomm);
					fprintfqt(stderr, "Send 3: %d to %d\n", GRID->iam, Alcol_process);
				}
			}

			//Finish the (decreasing-ring) broadcast of the solution block in previous process column
			if((myrow != rowprev) && (myrow != MModAdd1(rowprev, nprow)))
			{
				(void) HPL_send(Xdprev, kbprev, MModSub1(myrow, nprow), Cmsgid, Ccomm);
				fprintfqt(stderr, "Send 4: %d to %d\n", GRID->iam, MModSub1(myrow, nprow));
			}
		}
		else if (mycol == Alcol_process)
		{
			//Current column receives and accumulates partial update of previous solution block
			if (n1pprev > 0)
			{
				(void) HPL_recv(W, n1pprev, colprev, Rmsgid, Rcomm);
				fprintfqt(stderr, "Recv 5: %d from %d\n", GRID->iam, colprev);
				HPL_daxpy(n1pprev, HPL_rone, W, 1, XC+Anpprev-n1pprev, 1);
			}
		}
		
		//Solve current diagonal block 
		if((mycol == Alcol_process) && (myrow == Alrow))
		{
			HPL_dtrsv(HplColumnMajor, HplUpper, HplNoTrans, HplNonUnit, kb, Aptr+Anp, lda, XC+Anp, 1);
			HPL_dcopy(kb, XC+Anp, 1, XR+Anq, 1);
		}

		//Finish previous update
		if((mycol == colprev) && ((tmp1 = Anpprev - n1pprev ) > 0))
		{
			HPL_dgemv(HplColumnMajor, HplNoTrans, tmp1, kbprev, -HPL_rone, Aprev, lda, Xdprev, 1, HPL_rone, XC, 1);
		}

		//Save info of current step and update info for the next step
		if (mycol == Alcol_process)
		{
			Xdprev = Xd;
			Aprev = Aptr;
		}
		if (myrow == Alrow)
		{
			Anpprev -= kb;
		}
		rowprev = Alrow;
		colprev = Alcol_process;
		n1pprev = n1p;
		kbprev  = kb; n -= kb;
		Alrow = MModSub1(Alrow, nprow);
		Alcol_matrix--;
		Alcol_process = MColBlockToPCol(Alcol_matrix, 0, GRID);
		tmp1 = n - (kb = nb);
		tmp1 -= (tmp2 = Mmin(tmp1, n1));
		MnumrowI(n1p, tmp2, Mmax(0, tmp1), nb, myrow, nprow);

		Rmsgid = (Rmsgid+2 > MSGID_END_PTRSV ? MSGID_BEGIN_PTRSV : Rmsgid+2);
		Cmsgid = (Cmsgid+2 > MSGID_END_PTRSV ? MSGID_BEGIN_PTRSV+1 : Cmsgid+2);
	}
	//Replicate last solution block
	if (mycol == colprev)
	{
		(void) HPL_broadcast((void *) (XR), kbprev, HPL_DOUBLE, rowprev, Ccomm);
		fprintfqt(stderr, "Bcase 6: %d to %d\n", rowprev, GRID->iam);
	}

	if (Wfr) free(W);
	HPL_ptimer_detail(HPL_TIMING_PTRSV);

	//End of HPL_pdtrsv
}
