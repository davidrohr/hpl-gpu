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

#ifndef HPL_BLAS_H
#define HPL_BLAS_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_misc.h"
/*
 * ---------------------------------------------------------------------
 * typedef definitions
 * ---------------------------------------------------------------------
 */
enum HPL_ORDER
{  HplRowMajor = 101,  HplColumnMajor  = 102 };
enum HPL_TRANS
{  HplNoTrans  = 111,  HplTrans        = 112,  HplConjTrans    = 113 };
enum HPL_UPLO
{  HplUpper    = 121,  HplLower        = 122 };
enum HPL_DIAG
{  HplNonUnit  = 131,  HplUnit         = 132 };
enum HPL_SIDE
{  HplLeft     = 141,  HplRight        = 142 }; 

/*
 * ---------------------------------------------------------------------
 * The C interface of the BLAS is available ...
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define    CBLAS_INDEX         int
 
#define    CBLAS_ORDER         HPL_ORDER
#define    CblasRowMajor       HplRowMajor
#define    CblasColMajor       HplColMajor
 
#define    CBLAS_TRANSPOSE     HPL_TRANS
#define    CblasNoTrans        HplNoTrans
#define    CblasTrans          HplTrans
#define    CblasConjTrans      HplConjTrans
 
#define    CBLAS_UPLO          HPL_UPLO
#define    CblasUpper          HplUpper
#define    CblasLower          HplLower
 
#define    CBLAS_DIAG          HPL_DIAG
#define    CblasNonUnit        HplNonUnit
#define    CblasUnit           HplUnit
 
#define    CBLAS_SIDE          HPL_SIDE
#define    CblasLeft           HplLeft
#define    CblasRight          HplRight
/*
 * ---------------------------------------------------------------------
 * CBLAS Function prototypes
 * ---------------------------------------------------------------------
 */
CBLAS_INDEX cblas_idamax(  const int,       const double *,  const int );
void cblas_dswap (  const int,       double *,        const int,       double *,   const int );
void cblas_dcopy (  const int,       const double *,  const int,       double *,   const int );
void cblas_daxpy (  const int,       const double,    const double *,  const int,   double *,        const int );
void cblas_dscal (  const int,       const double,    double *,        const int );

void cblas_dgemv (  const enum CBLAS_ORDER,           const enum CBLAS_TRANSPOSE, const int,       const int,       const double,    const double *,
   const int,       const double *,  const int,       const double,   double *,        const int );

void cblas_dger (  const enum CBLAS_ORDER,           const int,       const int,   const double,    const double *,  const int,       const double *,
   const int,       double *,        const int );
void cblas_dtrsv (  const enum CBLAS_ORDER,           const enum CBLAS_UPLO,   const enum CBLAS_TRANSPOSE,       const enum CBLAS_DIAG,
   const int,       const double *,  const int,       double *,   const int );

void cblas_dgemm (  const enum CBLAS_ORDER,           const enum CBLAS_TRANSPOSE,   const enum CBLAS_TRANSPOSE,       const int,       const int,
   const int,       const double,    const double *,  const int,   const double *,  const int,       const double,    double *,   const int );
void cblas_dtrsm (  const enum CBLAS_ORDER,           const enum CBLAS_SIDE,   const enum CBLAS_UPLO,            const enum CBLAS_TRANSPOSE,
   const enum CBLAS_DIAG,            const int,       const int,    const double,    const double *,  const int,       double *,   const int );
/*
 * ---------------------------------------------------------------------
 * HPL C BLAS macro definition
 * ---------------------------------------------------------------------
 */
#ifdef NO_TRACE_CALLS

#ifdef HPL_CALDGEMM_CBLAS_WRAPPER

#define blasint int
#include "caldgemm_cblas_wrapper.h"
#undef blasint

#define    HPL_dswap           cblas_dswap
#define    HPL_dcopy           cblas_dcopy
#define    HPL_daxpy           cblas_daxpya
#define    HPL_dscal           cblas_dscala
#define    HPL_idamax          cblas_idamax

#define    HPL_dgemv           cblas_dgemva
#define    HPL_dtrsv           cblas_dtrsv
#define    HPL_dger            cblas_dger

#define    HPL_dgemm           cblas_dgemma
#ifdef HPL_CALL_CALDGEMM
#define    HPL_gpu_dgemm       CALDGEMM_dgemm
#else
#define    HPL_gpu_dgemm       cblas_dgemma
#endif
#define    HPL_dtrsm           cblas_dtrsma
#else
#define    HPL_dswap           cblas_dswap
#define    HPL_dcopy           cblas_dcopy
#define    HPL_daxpy           cblas_daxpy
#define    HPL_dscal           cblas_dscal
#define    HPL_idamax          cblas_idamax

#define    HPL_dgemv           cblas_dgemv
#define    HPL_dtrsv           cblas_dtrsv
#define    HPL_dger            cblas_dger

#define    HPL_dgemm           cblas_dgemm
#ifdef HPL_CALL_CALDGEMM
#define    HPL_gpu_dgemm       CALDGEMM_dgemm
#else
#define    HPL_gpu_dgemm       cblas_dgemm
#endif
#define    HPL_dtrsm           cblas_dtrsm
#endif

#ifdef HPL_CALDGEMM_ASYNC_FACT_DGEMM
#undef     HPL_dgemm
#define    HPL_dgemm           CALDGEMM_async_dgemm
#endif
#ifdef HPL_CALDGEMM_ASYNC_DTRSM
#undef     HPL_dtrsm
#define    HPL_dtrsm           CALDGEMM_async_dtrsm
#endif
#endif


/*
 * ---------------------------------------------------------------------
 * HPL BLAS Function prototypes
 * ---------------------------------------------------------------------
 */
#if defined(TRACE_CALLS)

int HPL_idamax(   const int,   const double *,   const int);
void HPL_daxpy(   const int,   const double,   const double *,   const int,   double *,   const int );
void HPL_dcopy(   const int,   const double *,   const int,   double *,   const int );
void HPL_dscal(   const int,   const double,   double *,   const int );
void HPL_dswap(   const int,   double *,   const int,   double *,   const int );
void HPL_dgemv(   const enum HPL_ORDER,   const enum HPL_TRANS,   const int,   const int,   const double,   const double *,   const int,   const double *,   const int,
   const double,   double *,   const int );
void HPL_dger(   const enum HPL_ORDER,   const int,   const int,   const double,   const double *,   const int,   double *,   const int,   double *,   const int );
void HPL_dtrsv(   const enum HPL_ORDER,   const enum HPL_UPLO,   const enum HPL_TRANS,   const enum HPL_DIAG,   const int,   const double *,   const int,   double *,
   const int );
void HPL_dgemm(   const enum HPL_ORDER,   const enum HPL_TRANS,   const enum HPL_TRANS,   const int,   const int,   const int,   const double,   const double *,
   const int,   const double *,   const int,   const double,   double *,   const int );
void HPL_gpu_dgemm(   const enum HPL_ORDER,   const enum HPL_TRANS,   const enum HPL_TRANS,   const int,   const int,   const int,   const double,   const double *,
   const int,   const double *,   const int,   const double,   double *,   const int, int );
void HPL_dtrsm(   const enum HPL_ORDER,   const enum HPL_SIDE,   const enum HPL_UPLO,   const enum HPL_TRANS,   const enum HPL_DIAG,   const int,   const int,
   const double,   const double *,   const int,   double *,   const int );

#endif

#endif
/*
 * hpl_blas.h
 */
