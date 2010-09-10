/* 
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

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include "permutationhelper.h"
//#include <xmmintrin.h>
//#include <ammintrin.h>
#include <mm3dnow.h>
#include <iostream>
#include <cstdlib>

namespace
{
    int HPL_init_dlaswp00N()
    {
        const char *num_threads_string = getenv("LASWP_NUM_THREADS");
        int num_threads = tbb::task_scheduler_init::default_num_threads();
        if (num_threads_string) {
            num_threads = atoi(num_threads_string);
        }
        static tbb::task_scheduler_init init(num_threads);
        return 0;
    }

    int _HPL_init_dlaswp00N = HPL_init_dlaswp00N();

    template<typename T>
        static inline void swap(__restrict__ T &a, __restrict__ T &b)
        {
            register T tmp = a;
            a = b;
            b = tmp;
        }

    template<typename T> static inline T max(T a, T b) { return a > b ? a : b; }

    class HPL_dlaswp00N_impl
    {
        private:
            double *__restrict__ const A;
            const PermutationHelper &__restrict__ perm;
            const int LDA, permSize;
        public:
            inline HPL_dlaswp00N_impl(double *_A, const int _LDA, const int _permSize, const PermutationHelper &_perm)
                : A(_A), perm(_perm), LDA(_LDA), permSize(_permSize)
            {}

            inline void operator()(const tbb::blocked_range<size_t> &range) const {
                for (size_t colIndex = range.begin(); colIndex < range.end(); ++colIndex) {
                    double *__restrict__ col = &A[colIndex * LDA];
                    int i = 0;
                    for (; i < permSize - 8; i += 2) {
                        //_m_prefetchw(&col[p[8].a]);
                        //_m_prefetchw(&col[p[8].b]);
                        swap(col[perm[i + 0].a], col[perm[i + 0].b]);
                        swap(col[perm[i + 1].a], col[perm[i + 1].b]);
                    }
                    for (; i < permSize; ++i) {
                        const int rowIndex = perm[i].a;
                        const int otherRow = perm[i].b;
                        swap(col[rowIndex], col[otherRow]);
                    }
                }
            }
    };
}

/**
 * Performs a series of local row interchanges on a matrix A.
 * One row interchange is initiated for rows 0 through M-1 of A.
 *
 * \param M
 *         On entry, M specifies the number of rows of the array A to be
 *         interchanged. M must be at least zero.
 *
 * \param N
 *         On entry, N  specifies  the number of columns of the array A.
 *         N must be at least zero.
 *
 * \param A
 *         On entry, A  points to an array of dimension (LDA,N) to which
 *         the row interchanges will be  applied.  On exit, the permuted
 *         matrix.
 *
 * \param LDA
 *         On entry, LDA specifies the leading dimension of the array A.
 *         LDA must be at least MAX(1,M).
 *
 * \param IPIV
 *         On entry,  IPIV  is  an  array of size  M  that  contains the
 *         pivoting  information.  For  k  in [0..M),  IPIV[k]=IROFF + l
 *         implies that local rows k and l are to be interchanged.
 */
extern "C" void HPL_dlaswp00N(const int M, const int N, double *__restrict__ A, const int LDA, const int *__restrict__ IPIV)
{
    // A is stored as
    // r0c0 r1c0 r2c0 ... rM-1c0 ... rLDA-1c0 r0c1
    // M   : #rows
    // LDA : offset to go to the next column (same row)
    // N   : #columns
    if ((M <= 0) || (N <= 0)) {
        return;
    }
    if (M * N < 128 * 128) {
        for (int rowIndex = 0; rowIndex < M; ++rowIndex) {
            const int otherRow = IPIV[rowIndex];
            double *__restrict__ col = A;
            if (otherRow != rowIndex) {
                int colIndex;
                for (colIndex = 0; colIndex < N - 15; colIndex += 16) {
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                }
                for (; colIndex < N; ++colIndex) {
                    swap(col[rowIndex], col[otherRow]); col += LDA;
                }
            }
        }
        return;
    }

    PermutationHelper &perm = PermutationHelper::instance();
    perm.ensureSize(M);
    int permSize = 0;
    for (int rowIndex = 0; rowIndex < M; ++rowIndex) {
        const int otherRow = IPIV[rowIndex];
        if (otherRow != rowIndex) {
            perm[permSize].a = rowIndex;
            perm[permSize].b = otherRow;
            ++permSize;
        }
    }

    // number of columns to process per thread: should result in worksets of 64 kB == L1d
    const unsigned int chunksize = max(1lu, 64 * 1024 / sizeof(double) / permSize);

    tbb::parallel_for (tbb::blocked_range<size_t>(0, N, chunksize),
            HPL_dlaswp00N_impl(A, LDA, permSize, perm)
            );
}
