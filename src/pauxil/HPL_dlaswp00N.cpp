/*
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

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include "permutationhelper.h"
#include <mm3dnow.h>
#include <iostream>
#include <cstdlib>
#include "helpers.h"

#include "util_timer.h"
#include "util_trace.h"

namespace
{
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
START_TRACE( DLASWP00N )

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
    } else {
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
#endif

END_TRACE
}
