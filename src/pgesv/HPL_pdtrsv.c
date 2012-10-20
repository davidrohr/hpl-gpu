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

int HPL_n1(int matrix_col, int nb, const HPL_T_grid* GRID)
{
	int cols = 0;
	if (matrix_col)
	{
		int process_col = GRID->col_mapping[matrix_col];
		matrix_col--;
		while (matrix_col >= 0 && GRID->col_mapping[matrix_col] != process_col)
		{
			cols++;
			matrix_col--;
		}
	}
	if (cols == 0) cols = 1;
	return(cols * nb);
}

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
		colprev, kb, kbprev, lda, mycol, myrow, n, n1, n1p, n1pprev=0, nb, npcol, nprow, rowprev, tmp1, tmp2, Wsize;
	int sendcol_matrix = -1;

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

	tmp1 = (n - 1) / nb;
	Alrow = tmp1 % nprow;
	Alcol_matrix = tmp1;
	Alcol_process = MColBlockToPCol(Alcol_matrix, GRID);
	kb = n - tmp1 * nb;

	Aptr = (double *) (A);
	XC = Mptr(Aptr, 0, Anq, lda);
	Mindxg2p_col(n, nb, nb, Bcol, 0, GRID);

	if ((Anp > 0) && (Alcol_process != Bcol))
	{
		if(mycol == Bcol)
		{
			(void) HPL_send(XC, Anp, Alcol_process, Rmsgid, Rcomm);
		}
		else if(mycol == Alcol_process)
		{
			(void) HPL_recv(XC, Anp, Bcol, Rmsgid, Rcomm);
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
	//n1 = (npcol - 1) * nb;
	//n1 = Mmax(n1, nb);
	n1 = HPL_n1(Alcol_matrix, nb, GRID);
	Wsize = Mmin((npcol - 1) * nb, Anp);
	if (Wsize > 0)
	{
		W = (double*) malloc((size_t) Wsize * sizeof(double));
		if (W == NULL)
		{
			HPL_pabort(__LINE__, "HPL_pdtrsv", "Memory allocation failed");
		}
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
			//fprintfqt(STD_OUT, "Process %d: dtrsv, offset %d\n", GRID->iam, Anp);
			HPL_dcopy(kb, XC+Anp, 1, Xd, 1);
		}
	}

	n -= kb;

// Start the operations
	while(n > 0)
	{
		rowprev = Alrow;
		Alrow = MModSub1(Alrow, nprow);
		colprev = Alcol_process;
		Alcol_matrix--;
		Alcol_process = MColBlockToPCol(Alcol_matrix, GRID);
		kbprev = kb;
		n1 = HPL_n1(Alcol_matrix, nb, GRID);
		tmp1 = n - (kb = nb);
		tmp1 -= (tmp2 = Mmin(tmp1, n1));
		MnumrowI(n1p, tmp2, Mmax(0, tmp1), nb, myrow, nprow);
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
				}
			}
			else
			{
				(void) HPL_recv(Xdprev, kbprev, MModAdd1(myrow, nprow), Cmsgid, Ccomm);
			}
		}

		if (Alcol_process < colprev ? (mycol <= colprev && mycol > Alcol_process) : (mycol <= colprev || mycol > Alcol_process))
		{
			if (mycol == colprev)
			{
				//Compute partial update of previous solution block and send it to current column
				tmp1 = Anpprev - n1pprev;
				HPL_dgemv(HplColumnMajor, HplNoTrans, n1pprev, kbprev, -HPL_rone, Aprev+tmp1, lda, Xdprev, 1, HPL_rone, XC+tmp1, 1 );
				//fprintfqt(STD_OUT, "Process %d: dgemv %d rows starting from %d\n", GRID->iam, n1pprev, tmp1);
				sendcol_matrix = Alcol_matrix;
			}
			
			if(GridIsNotPx1)
			{
				if (sendcol_matrix != -1)
				{
					tmp2 = 1;
					while (sendcol_matrix >= tmp2 && !(GRID->col_mapping[sendcol_matrix - tmp2 + 1] > GRID->col_mapping[sendcol_matrix - tmp2] ?
						(GRID->col_mapping[sendcol_matrix - tmp2 + 1] >= mycol && GRID->col_mapping[sendcol_matrix - tmp2] <= mycol) :
						(GRID->col_mapping[sendcol_matrix - tmp2 + 1] >= mycol || GRID->col_mapping[sendcol_matrix - tmp2] <= mycol)))
					{
						tmp2++;
					}
					sendcol_matrix -= tmp2;
					tmp2 *= nb;
					MnumrowI(tmp1, tmp2, n - tmp2, nb, myrow, nprow);
					tmp2 = Anpprev - tmp1;
				}
				else
				{
					tmp2 = 0;
					tmp1 = 0;
				}
				
				//fprintfqt(STD_OUT, "Process %d: sending to %d (%d bytes starting from %d, partial update)\n", GRID->iam, Alcol_process, tmp1, tmp2);
				(void) MPI_Send(XC+tmp2, tmp1, MPI_DOUBLE, Alcol_process, Rmsgid, Rcomm);
			}
		}

		if (mycol == colprev)
		{
			//Finish the (decreasing-ring) broadcast of the solution block in previous process column
			if((myrow != rowprev) && (myrow != MModAdd1(rowprev, nprow)))
			{
				//(void) HPL_send(Xdprev, kbprev, MModSub1(myrow, nprow), Cmsgid, Ccomm);
				MPI_Send(Xdprev, kbprev, MPI_DOUBLE, MModSub1(myrow, nprow), Cmsgid, Ccomm);
			}
		}
		else if (mycol == Alcol_process)
		{
			//Current column receives and accumulates partial update of previous solution block
			for (int i = colprev;(i - mycol) % npcol != 0;i = (i + npcol - 1) % npcol)
			{
				MPI_Status tmpstatus;
				int recvsize;
				//fprintfqt(STD_OUT, "Process %d starting receive from %d (buffer %d)\n", GRID->iam, i, Wsize);
				MPI_Recv(W, Wsize, MPI_DOUBLE, i, Rmsgid, Rcomm, &tmpstatus);
				MPI_Get_count(&tmpstatus, MPI_DOUBLE, &recvsize);
				//fprintfqt(STD_OUT, "Process %d: received from %d (%d bytes starting from %d)\n", GRID->iam, i, recvsize, Anpprev - recvsize);
				HPL_daxpy(recvsize, HPL_rone, W, 1, XC+Anpprev-recvsize, 1);
			}
		}
		
		//Solve current diagonal block 
		if((mycol == Alcol_process) && (myrow == Alrow))
		{
			HPL_dtrsv(HplColumnMajor, HplUpper, HplNoTrans, HplNonUnit, kb, Aptr+Anp, lda, XC+Anp, 1);
			//fprintfqt(STD_OUT, "Process %d: dtrsv, offset %d\n", GRID->iam, Anp);
			HPL_dcopy(kb, XC+Anp, 1, XR+Anq, 1);
		}

		//Finish previous update
		if((mycol == colprev) && ((tmp1 = Anpprev - n1pprev ) > 0))
		{
			HPL_dgemv(HplColumnMajor, HplNoTrans, tmp1, kbprev, -HPL_rone, Aprev, lda, Xdprev, 1, HPL_rone, XC, 1);
			//fprintfqt(STD_OUT, "Process %d: dgemv (%d rows starting from %d, finishing)\n", GRID->iam, tmp1, 0);
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
		n1pprev = n1p;
		n -= kb;
		Rmsgid = (Rmsgid+2 > MSGID_END_PTRSV ? MSGID_BEGIN_PTRSV : Rmsgid+2);
		Cmsgid = (Cmsgid+2 > MSGID_END_PTRSV ? MSGID_BEGIN_PTRSV+1 : Cmsgid+2);
	}
	rowprev = Alrow;
	colprev = Alcol_process;
	kbprev = kb;

	//Replicate last solution block
	if (mycol == colprev)
	{
		(void) HPL_broadcast((void *) (XR), kbprev, HPL_DOUBLE, rowprev, Ccomm);
	}

	if (Wsize) free(W);
	HPL_ptimer_detail(HPL_TIMING_PTRSV);

	//End of HPL_pdtrsv
}
