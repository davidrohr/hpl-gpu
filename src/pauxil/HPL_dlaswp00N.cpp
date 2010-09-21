/*
    Copyright (C) 2010 Matthias Kretz <kretz@kde.org>
    Copyright (C) 2010 Frankfurt Institute for Advanced Studies (FIAS)

    The source code is property of the Frankfurt Institute for Advanced Studies
    (FIAS). None of the material may be copied, reproduced, distributed,
    republished, downloaded, displayed, posted or transmitted in any form or by
    any means, including, but not limited to, electronic, mechanical,
    photocopying, recording, or otherwise, without the prior written permission
    of FIAS.

*/

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include "permutationhelper.h"
#include <mm3dnow.h>
#include <iostream>
#include <cstdlib>

#include "util_timer.h"
#include "util_trace.h"

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

extern "C" void HPL_dlaswp00N(const int M, const int N, double *__restrict__ A, const int LDA, const int *__restrict__ IPIV)
{
#ifdef TRACE_CALLS
   uint64_t tr_start, tr_end, tr_diff;
   tr_start = util_getTimestamp();
#endif /* TRACE_CALLS */

#ifdef USE_ORIGINAL_LASWP
#include "HPL_dlaswp00N.c"
#else
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
#endif

#ifdef TRACE_CALLS
   tr_end = util_getTimestamp();
   tr_diff = util_getTimeDifference( tr_start, tr_end );

   fprintf( trace_dgemm, "DLASWP00N,M=%i,N=%i,LDA=%i,TIME=%lu\n", M, N, LDA, tr_diff );
#endif /* TRACE_CALLS */
}
