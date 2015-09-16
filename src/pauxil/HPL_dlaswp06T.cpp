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
class dlaswp06T_impl
{
    private:
        const size_t M;
        const size_t LDA;
        const size_t LDU;
        double *__restrict__ const A;
        double *__restrict__ const U;
        const int *__restrict__ const LINDXA;
    public:
        dlaswp06T_impl(size_t _M, double *_A, size_t _LDA, double *_U, size_t _LDU,
                const int *_LINDXA)
            : M(_M), LDA(_LDA), LDU(_LDU), A(_A), U(_U), LINDXA(_LINDXA) {}
        void operator()(const tbb::blocked_range<size_t> &range) const
        {
            const size_t begin = range.begin();
            const size_t columns = range.end() - begin;
            double *u = &U[begin];
            double *uNext = u + LDU;
            double *A2 = &A[begin * LDA];
            for (size_t i = 0; i < M; ++i) {
                double *a = &A2[LINDXA[i]];
                ptrdiff_t aNext = &A2[LINDXA[i + 1]] - a;
                size_t j = 7;
                for (; j < columns; j += 8) {
                    _m_prefetchw(&uNext[j - 7]);
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 7]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 6]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 5]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 4]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 3]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 2]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 1]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j    ]); a += LDA;
                }
                _m_prefetchw(&uNext[j - 7]);
                double *u2 = &u[columns];
                j = (j - columns) * 0x15;
                long tmp1;
#ifdef HPL_HAVE_PREFETCHW
#define HPL_PREFETCHW "prefetchw"
#else
#define HPL_PREFETCHW "prefetchnta"
#endif
                asm volatile(
                    "lea 0x6(%%rip), %[tmp1]\n"
                    "lea (%[tmp1], %[tmp0], 1), %[tmp1]\n"
                    "jmpq *%[tmp1]\n"

                    HPL_PREFETCHW " (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x38(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x38(%[u])\n"
                    "add %[LDA], %[a]\n"

                    HPL_PREFETCHW " (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x30(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x30(%[u])\n"
                    "add %[LDA], %[a]\n"

                    HPL_PREFETCHW " (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x28(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x28(%[u])\n"
                    "add %[LDA], %[a]\n"

                    HPL_PREFETCHW " (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x20(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x20(%[u])\n"
                    "add %[LDA], %[a]\n"

                    HPL_PREFETCHW " (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x18(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x18(%[u])\n"
                    "add %[LDA], %[a]\n"

                    HPL_PREFETCHW " (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x10(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x10(%[u])\n"
                    "add %[LDA], %[a]\n"

                    HPL_PREFETCHW " (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x8(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x8(%[u])\n"
                    "add %[LDA], %[a]\n"
                    : [a]"+d"(a),
                      [tmp0]"+S"(j),
                      [tmp1]"=a"(tmp1)
                    : [u]"D"(u2),
                      [aNext]"c"(aNext * 8),
                      [LDA]"b"(LDA * 8)
                   );
                u = uNext;
                uNext += LDU;
            }
#ifndef HPL_LASWP_AVX
            _mm_mfence();
#endif
        }
};
#endif /* USE_ORIGINAL_LASWP */

extern "C" void HPL_dlaswp06T(const int M, const int N, double *A,
        const int LDA, double *U, const int LDU, const int *LINDXA)
{
START_TRACE( DLASWP06T )

#ifdef USE_ORIGINAL_LASWP
#include "HPL_dlaswp06T.c"
#else
    if (M <= 0 || N <= 0) {
        return;
    }

    tbb::parallel_for (tbb::blocked_range<size_t>(0, N, 32),
            dlaswp06T_impl(M, A, LDA, U, LDU, LINDXA));
#endif

END_TRACE
#ifdef TRACE_LASWP
   char filename[256];
   snprintf(filename, 256, "dlaswp06T.%04d.%05d.%05d.%05d.%7.4fs.dat", M, N, LDA, LDU, laswp_time);
   FILE *permdata = fopen(filename, "w");
   fwrite(LINDXA, sizeof(LINDXA[0]), M, permdata);
   fclose(permdata);
#endif
}
