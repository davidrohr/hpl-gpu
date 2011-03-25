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
#include "util_trace.h"

#ifndef HPL_dtrsm

void HPL_dtrsm(const enum HPL_ORDER ORDER, const enum HPL_SIDE SIDE, const enum HPL_UPLO UPLO, const enum HPL_TRANS TRANS, const enum HPL_DIAG DIAG, const int M, const int N,
   const double ALPHA, const double *A, const int LDA, double *B, const int LDB)
{
/* 
 * Purpose
 * =======
 *
 * HPL_dtrsm solves one of the matrix equations
 *  
 *    op( A ) * X = alpha * B,   or  X * op( A ) = alpha * B,
 *  
 * where alpha is a scalar, X and B are m by n matrices, A is a unit, or
 * non-unit, upper or lower triangular matrix and op(A) is one of
 *  
 *    op( A ) = A   or   op( A ) = A^T.
 *  
 * The matrix X is overwritten on B.
 *  
 * No test for  singularity  or  near-singularity  is included  in  this
 * routine. Such tests must be performed before calling this routine.
 *
 * Arguments
 * =========
 *
 * ORDER   (local input)                 const enum HPL_ORDER
 *         On entry, ORDER  specifies the storage format of the operands
 *         as follows:                                                  
 *            ORDER = HplRowMajor,                                      
 *            ORDER = HplColumnMajor.                                   
 *
 * SIDE    (local input)                 const enum HPL_SIDE
 *         On entry, SIDE  specifies  whether  op(A) appears on the left
 *         or right of X as follows:
 *            SIDE==HplLeft    op( A ) * X = alpha * B,
 *            SIDE==HplRight   X * op( A ) = alpha * B.
 *
 * UPLO    (local input)                 const enum HPL_UPLO
 *         On  entry,   UPLO   specifies  whether  the  upper  or  lower
 *         triangular  part  of the array  A  is to be referenced.  When
 *         UPLO==HplUpper, only  the upper triangular part of A is to be
 *         referenced, otherwise only the lower triangular part of A is 
 *         to be referenced. 
 *
 * TRANS   (local input)                 const enum HPL_TRANS
 *         On entry, TRANSA  specifies the form of  op(A)  to be used in
 *         the matrix-matrix operation follows:                         
 *            TRANSA==HplNoTrans    : op( A ) = A,                     
 *            TRANSA==HplTrans      : op( A ) = A^T,                   
 *            TRANSA==HplConjTrans  : op( A ) = A^T.                   
 *
 * DIAG    (local input)                 const enum HPL_DIAG
 *         On entry,  DIAG  specifies  whether  A  is unit triangular or
 *         not. When DIAG==HplUnit,  A is assumed to be unit triangular,
 *         and otherwise, A is not assumed to be unit triangular.
 *
 * M       (local input)                 const int
 *         On entry,  M  specifies  the number of rows of the  matrix B.
 *         M must be at least zero.
 *
 * N       (local input)                 const int
 *         On entry, N  specifies the number of columns of the matrix B.
 *         N must be at least zero.
 *
 * ALPHA   (local input)                 const double
 *         On entry, ALPHA specifies the scalar alpha.   When  ALPHA  is
 *         supplied  as  zero then the elements of the matrix B need not
 *         be set on input.
 *
 * A       (local input)                 const double *
 *         On entry,  A  points  to an array of size equal to or greater
 *         than LDA * k,  where  k is m  when  SIDE==HplLeft  and  is  n
 *         otherwise.  Before  entry  with  UPLO==HplUpper,  the leading
 *         k by k upper triangular  part of the array A must contain the
 *         upper triangular  matrix and the  strictly  lower  triangular
 *         part of A is not referenced.  When  UPLO==HplLower on  entry,
 *         the  leading k by k lower triangular part of the array A must
 *         contain the lower triangular matrix  and  the  strictly upper
 *         triangular part of A is not referenced.
 *          
 *         Note that  when  DIAG==HplUnit,  the  diagonal elements of  A
 *         not referenced  either,  but are assumed to be unity.
 *
 * LDA     (local input)                 const int
 *         On entry,  LDA  specifies  the  leading  dimension  of  A  as
 *         declared  in  the  calling  (sub) program.  LDA  must  be  at
 *         least MAX(1,m) when SIDE==HplLeft, and MAX(1,n) otherwise.
 *
 * B       (local input/output)          double *
 *         On entry,  B  points  to an array of size equal to or greater
 *         than LDB * n.  Before entry, the leading  m by n  part of the
 *         array B must contain the matrix  B, except when beta is zero,
 *         in which case B need not be set on entry.  On exit, the array
 *         B is overwritten by the m by n solution matrix.
 *
 * LDB     (local input)                 const int
 *         On entry,  LDB  specifies  the  leading  dimension  of  B  as
 *         declared  in  the  calling  (sub) program.  LDB  must  be  at
 *         least MAX(1,m).
 *
 * ---------------------------------------------------------------------
 */ 
START_TRACE( DTRSM )

   cblas_dtrsm( ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB );

END_TRACE
/*
 * End of HPL_dtrsm
 */
}

#endif
