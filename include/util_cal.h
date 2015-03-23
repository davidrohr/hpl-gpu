/**
 * A wrapper for the C++ CALDGEMM
 *
 * Copyright 2010:
 *  - David Rohr (drohr@jwdt.org)
 *  - Matthias Bach (bach@compeng.uni-frankfurt.de)
 *  - Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 *
 * This file is part of HPL-GPU.
 *
 * HPL-GPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HPL-GPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HPL-GPU.  If not, see <http://www.gnu.org/licenses/>.
 *
 * In addition to the rules layed out by the GNU General Public License
 * the following exception is granted:
 *
 * Use with the Original BSD License.
 *
 * Notwithstanding any other provision of the GNU General Public License
 * Version 3, you have permission to link or combine any covered work with
 * a work licensed under the 4-clause BSD license into a single combined
 * work, and to convey the resulting work.  The terms of this License will
 * continue to apply to the part which is the covered work, but the special
 * requirements of the 4-clause BSD license, clause 3, concerning the
 * requirement of acknowledgement in advertising materials will apply to
 * the combination as such.
 */

#ifndef UTIL_CAL_H
#define UTIL_CAL_H

#ifdef __cplusplus
extern "C"
{
#endif

void* CALDGEMM_GetObject();
int CALDGEMM_Init();
void CALDGEMM_Shutdown();
void CALDGEMM_dgemm( const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE,
                     const enum CBLAS_TRANSPOSE, const int, const int,
   const int,       const double,    const double *,  const int,
   const double *,  const int,       const double,    double *,
   const int, int LinpackCallbacks );
void CALDGEMM_async_dgemm( const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE,
                     const enum CBLAS_TRANSPOSE, const int, const int,
   const int,       const double,    const double *,  const int,
   const double *,  const int,       const double,    double *,
   const int );
void CALDGEMM_async_dtrsm(const enum HPL_ORDER ORDER, const enum HPL_SIDE SIDE, const enum HPL_UPLO UPLO, const enum HPL_TRANS TRANS, const enum HPL_DIAG DIAG, const int M, const int N,
   const double ALPHA, const double *A, const int LDA, double *B, const int LDB);
void CALDGEMM_async_dtrsm2(const enum HPL_ORDER ORDER, const enum HPL_SIDE SIDE, const enum HPL_UPLO UPLO, const enum HPL_TRANS TRANS, const enum HPL_DIAG DIAG, const int M, const int N,
   const double ALPHA, const double *A, const int LDA, double *B, const int LDB);
void* CALDGEMM_alloc(size_t size, int interleave);
void CALDGEMM_free(void* ptr);
void CALDGEMM_set_num_nodes(int num, int rank);
void CALDGEMM_enable_async_laswp(int enable);
void CALDGEMM_reset();

#ifdef __cplusplus
}
#endif

#endif


