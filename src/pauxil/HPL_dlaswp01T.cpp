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

#include <cstddef>
#include "util_timer.h"
#include "util_trace.h"

#ifndef USE_ORIGINAL_LASWP
#include <tbb/parallel_for.h>
#include "helpers.h"
#endif

#ifndef USE_ORIGINAL_LASWP
class dlaswp01T_impl
{
    private:
        const size_t M, LDA, LDU;
        double *const __restrict__ A;
        double *const __restrict__ U;
        const int *const __restrict__ LINDXA;
        const int *const __restrict__ LINDXAU;

    public:
        dlaswp01T_impl(const int _M, double *_A, const unsigned int _LDA,
                double *_U, const unsigned int _LDU, const int *_LINDXA, const int *_LINDXAU)
            : M(_M), LDA(_LDA), LDU(_LDU), A(_A), U(_U), LINDXA(_LINDXA), LINDXAU(_LINDXAU)
        {}

        void operator()(const MyRange<16, 64> &range) const
        {
            const size_t begin = range.begin();
            size_t N = range.N();
            double *__restrict__ AA = &A[begin * LDA];
            double *__restrict__ UU = &U[begin];
            for (; N >= 16; AA += 16 * LDA, UU += 16, N -= 16) {
                for (size_t i = 0; i < M; ++i) {
                    const double *__restrict__ Ar = &AA[LINDXA[i]];
                    const ptrdiff_t ArNext = &AA[LINDXA[(i + 1 == M) ? 0 : (i + 1)]] - Ar;
                    if (LINDXAU[i] >= 0) {
                        const size_t rowUw = LINDXAU[i];
                        double *__restrict__ Uw = &UU[rowUw * LDU];
                        _m_prefetchw(Uw + 0);
                        _m_prefetchw(Uw + 8);
                        _m_prefetchw(Uw + 15);
                        // from the LINDXAU data I've seen it can be expected that:
                        // rowUw == LINDXAU[i] => rowUw + 1 == LINDXAU[i + 1]
                        _m_prefetchw(Uw + LDU + 0);
                        _m_prefetchw(Uw + LDU + 8);
                        _m_prefetchw(Uw + LDU + 15);
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 0], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 1], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 2], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 3], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 4], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 5], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 6], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 7], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 8], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[ 9], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[10], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[11], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[12], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[13], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[14], Ar); Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[15], Ar);
                    } else {
                        const size_t rowAw = -LINDXAU[i];
                        double *__restrict__ Aw = &AA[rowAw];
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar);
                    }
                }
            }
        }
};
#endif /* USE_ORIGINAL_LASWP */

extern "C" void HPL_dlaswp01T(const int M, const int N, double *A, const int LDA,
        double *U, const int LDU, const int *LINDXA, const int *LINDXAU)
{
START_TRACE( DLASWP01T )

#ifdef USE_ORIGINAL_LASWP
#include "HPL_dlaswp01T.c"
#else
    const size_t largeN = N & ~15;
    const size_t smallN = N - largeN;
    tbb::parallel_for (MyRange<16, 64>(0, largeN), dlaswp01T_impl(M, A, LDA, U, LDU, LINDXA, LINDXAU), tbb::auto_partitioner());

    if (__builtin_expect(smallN > 0, 1)) {
        for (int i = 0; i < M - 1; ++i) {
            const double *__restrict__ Ar = &A[LINDXA[i] + largeN * LDA];
            const ptrdiff_t ArNext = &A[LINDXA[i + 1] + largeN * LDA] - Ar;
            if (LINDXAU[i] >= 0) {
                double *__restrict__ Uw = &U[LINDXAU[i] * static_cast<size_t>(LDU) + largeN];
                for (size_t col = 0; col < smallN; ++col) {
                    _mm_prefetch(Ar + ArNext, _MM_HINT_T0);
                    copy(&Uw[col], Ar);
                    Ar += LDA;
                }
            } else {
                double *__restrict__ Aw = &A[-LINDXAU[i] + largeN * LDA];
                for (size_t col = 0; col < smallN; ++col) {
                    _mm_prefetch(Ar + ArNext, _MM_HINT_T0);
                    streamingCopy(Aw, Ar);
                    Ar += LDA;
                    Aw += LDA;
                }
            }
        }
        const double *__restrict__ Ar = &A[LINDXA[M - 1] + largeN * LDA];
        if (LINDXAU[M - 1] >= 0) {
            double *__restrict__ Uw = &U[LINDXAU[M - 1] * static_cast<size_t>(LDU) + largeN];
            for (size_t col = 0; col < smallN; ++col) {
                copy(&Uw[col], Ar);
                Ar += LDA;
            }
        } else {
            double *__restrict__ Aw = &A[-LINDXAU[M - 1] + largeN * LDA];
            for (size_t col = 0; col < smallN; ++col) {
                streamingCopy(Aw, Ar);
                Ar += LDA;
                Aw += LDA;
            }
        }
    }
    _mm_sfence();

#endif

END_TRACE
#ifdef TRACE_CALLS
#ifdef TRACE_PERMDATA
   char filename[256];
   snprintf(filename, 256, "dlaswp01T.%04d.%05d.%05d.%05d.dat", M, N, LDA, LDU);
   FILE *permdata = fopen(filename, "w");
   fwrite(LINDXA, sizeof(LINDXA[0]), M, permdata);
   fwrite(LINDXAU, sizeof(LINDXAU[0]), M, permdata);
   fclose(permdata);
#endif
#endif /* TRACE_CALLS */
} 
