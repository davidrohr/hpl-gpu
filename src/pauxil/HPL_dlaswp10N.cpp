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

#include "util_timer.h"
#include "util_trace.h"
#include "helpers.h"
#include <tbb/parallel_for.h>

typedef MyRange<16, 96> Range10N;

class dlaswp10N_impl
{
    const size_t N, LDA;
    double *__restrict__ const A;
    const int *__restrict__ const IPIV;
    public:
        dlaswp10N_impl(size_t _N, double *_A, size_t _LDA, const int *_IPIV)
            : N(_N), LDA(_LDA), A(_A), IPIV(_IPIV)
        {}

        void operator()(const Range10N &range) const
        {
            const size_t begin = range.begin();
            const size_t M = range.N();

            for (size_t i = 0; i < N; ++i ) {
                const size_t ii = IPIV[i];
                if (__builtin_expect(i != ii, 1)) {
                    double *__restrict__ a0 = A + i  * LDA + begin;
                    double *__restrict__ a1 = A + ii * LDA + begin;

                    _m_prefetchw(&a0[ 0]);
                    _m_prefetchw(&a1[ 0]);
                    _m_prefetchw(&a0[ 8]);
                    _m_prefetchw(&a1[ 8]);
                    _m_prefetchw(&a0[16]);
                    _m_prefetchw(&a1[16]);
                    _m_prefetchw(&a0[24]);
                    _m_prefetchw(&a1[24]);

                    for (size_t j = M; j; j -= 16, a0 += 16, a1 += 16 ) {
                        _m_prefetchw(&a0[32]);
                        _m_prefetchw(&a1[32]);
#ifdef HPL_LASWP_AVX
                        swapAVX(&a0[ 0], &a1[ 0]);
                        swapAVX(&a0[ 4], &a1[ 4]);
#else
                        swapSSE(a0[ 0], a1[ 0]);
                        swapSSE(a0[ 2], a1[ 2]);
                        swapSSE(a0[ 4], a1[ 4]);
                        swapSSE(a0[ 6], a1[ 6]);
#endif
                        _m_prefetchw(&a0[40]);
                        _m_prefetchw(&a1[40]);
#ifdef HPL_LASWP_AVX
                        swapAVX(&a0[ 8], &a1[ 8]);
                        swapAVX(&a0[12], &a1[12]);
#else
                        swapSSE(a0[ 8], a1[ 8]);
                        swapSSE(a0[10], a1[10]);
                        swapSSE(a0[12], a1[12]);
                        swapSSE(a0[14], a1[14]);
#endif
                    }
                }
            }
#ifndef HPL_LASWP_AVX
            _mm_sfence();
#endif
        }
};

extern "C" void HPL_dlaswp10N(const int _M, const int N, double *A,
        const int LDA, const int *IPIV)
{
START_TRACE( DLASWP10N )

#ifdef USE_ORIGINAL_LASWP
const int M = _M;
#include "HPL_dlaswp10N.c"
#else
    // we require that LDA is even!

    size_t M = _M;

    if (__builtin_expect(M <= 0 || N <= 0, 0)) {
        return;
    }

    if (__builtin_expect((A - static_cast<double *>(0)) & 1, 0)) {
        // A does not start on a 16 Byte boundary

        // TODO: should be started as a separate task
        for (int i = 0; i < N; ++i) {
            const int ii = IPIV[i];
            if (__builtin_expect(i != ii, 1)) {
                double *__restrict__ a0 = A + i * LDA;
                double *__restrict__ a1 = A + ii * LDA;
                swap(*a0, *a1);
            }
        }
        A += 1; M -= 1;
    }

    if (M & 15) {
        // TODO: should be started as a separate task
        for (int i = 0; i < N; ++i) {
            const int ii = IPIV[i];
            if (__builtin_expect(i != ii, 1)) {
                double *__restrict__ a0 = A + (M & ~15) + i * LDA;
                double *__restrict__ a1 = A + (M & ~15) + ii * LDA;
                size_t j;
#ifdef HPL_LASWP_AVX
                for (j = 0; j < (M & 14); j += 4) {
                    swapAVX(&a0[j], &a1[j]);
                }
#else
                for (j = 0; j < (M & 14); j += 2) {
                    swapSSE(a0[j], a1[j]);
                }
#endif
                if (M & 1) {
                    swap(a0[j], a1[j]);
                }
            }
        }
    }

    tbb::parallel_for (Range10N(0, M & ~15), dlaswp10N_impl(N, A, LDA, IPIV));
#endif

END_TRACE
#ifdef TRACE_LASWP
   char filename[256];
   snprintf(filename, 256, "dlaswp10N.%04d.%05d.%05d.%7.4fs.dat", M, N, LDA, laswp_time);
   FILE *permdata = fopen(filename, "w");
   fwrite(IPIV, sizeof(IPIV[0]), N, permdata);
   fclose(permdata);
#endif
}
