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
                asm volatile(
                    "lea 0x6(%%rip), %[tmp1]\n"
                    "lea (%[tmp1], %[tmp0], 1), %[tmp1]\n"
                    "jmpq *%[tmp1]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x38(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x38(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x30(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x30(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x28(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x28(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x20(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x20(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x18(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x18(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x10(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x10(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
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
        }
};
#endif /* USE_ORIGINAL_LASWP */

extern "C" void HPL_dlaswp06T(const int M, const int N, double *A,
        const int LDA, double *U, const int LDU, const int *LINDXA)
{
#ifdef TRACE_CALLS
   uint64_t tr_start, tr_end, tr_diff;
   tr_start = util_getTimestamp();
#endif /* TRACE_CALLS */

#ifdef USE_ORIGINAL_LASWP
#include "HPL_dlaswp06T.c"
#else
    if (M <= 0 || N <= 0) {
        return;
    }

    tbb::parallel_for (tbb::blocked_range<size_t>(0, N, 32),
            dlaswp06T_impl(M, A, LDA, U, LDU, LINDXA));
#endif

#ifdef TRACE_CALLS
   tr_end = util_getTimestamp();
   tr_diff = util_getTimeDifference( tr_start, tr_end );

   fprintf( trace_dgemm, "DLASWP06T,M=%i,N=%i,LDA=%i,LDU=%i,TIME=%lu,THRPT=%.2fGB/s\n", M, N, LDA, LDU, tr_diff,
           0.004 * sizeof(double) * M * N / tr_diff );
#ifdef TRACE_PERMDATA
   char filename[256];
   snprintf(filename, 256, "dlaswp06T.%04d.%05d.%05d.%05d.dat", M, N, LDA, LDU);
   FILE *permdata = fopen(filename, "w");
   fwrite(LINDXA, sizeof(LINDXA[0]), M, permdata);
   fclose(permdata);
#endif
#endif /* TRACE_CALLS */
}
