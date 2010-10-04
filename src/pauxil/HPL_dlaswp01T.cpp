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

        void operator()(const tbb::blocked_range<size_t> &range) const
        {
            const size_t N = range.end() - range.begin();

            for (size_t i = 0; i < M - 1; ++i) {
                const double *__restrict__ Ar = &A[LINDXA[i]];
                const ptrdiff_t ArNext = &A[LINDXA[i + 1]] - Ar;
                if (LINDXAU[i] >= 0) {
                    const size_t rowUw = LINDXAU[i];
                    double *__restrict__ Uw = &U[rowUw * LDU];
                    double *__restrict__ const UwEnd = Uw + range.end();
                    Uw += range.begin();
                    // from the LINDXAU data I've seen it can be expected that:
                    // rowUw == LINDXAU[i] => rowUw + 1 == LINDXAU[i + 1]
                    for (size_t col = 0; col < N; col += 8) {
                        _m_prefetchw(Uw + col);
                        _m_prefetchw(Uw + LDU + col);
                    }
                    Ar += range.begin() * LDA;
                    do {
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
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); copy(&Uw[15], Ar); Ar += LDA;
                        Uw += 16;
                    } while ((UwEnd - Uw) >= 16);
                    for (;Uw < UwEnd; ++Uw) {
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1);
                        copy(Uw, Ar);
                        Ar += LDA;
                    }
                } else {
                    const size_t rowAw = -LINDXAU[i];
                    double *__restrict__ Aw = &A[rowAw];
                    size_t col = range.begin();
                    Aw += col * LDA;
                    Ar += col * LDA;
                    col = range.end() - col;
                    do {
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
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                        col -= 16;
                    } while (col >= 16);
                    for (; col; --col) {
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T1); streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    }
                }
            }
            const double *__restrict__ Ar = &A[LINDXA[M - 1]];
            if (LINDXAU[M - 1] >= 0) {
                const size_t rowUw = LINDXAU[M - 1];
                double *__restrict__ Uw = &U[rowUw * LDU];
                double *__restrict__ const UwEnd = Uw + range.end();
                Uw += range.begin();
                _m_prefetchw(Uw);
                _m_prefetchw(Uw + 8);
                _m_prefetchw(Uw + 16);
                Ar += range.begin() * LDA;
                do {
                    copy(&Uw[ 0], Ar); Ar += LDA;
                    copy(&Uw[ 1], Ar); Ar += LDA;
                    copy(&Uw[ 2], Ar); Ar += LDA;
                    copy(&Uw[ 3], Ar); Ar += LDA;
                    copy(&Uw[ 4], Ar); Ar += LDA;
                    copy(&Uw[ 5], Ar); Ar += LDA;
                    copy(&Uw[ 6], Ar); Ar += LDA;
                    copy(&Uw[ 7], Ar); Ar += LDA;
                    copy(&Uw[ 8], Ar); Ar += LDA;
                    copy(&Uw[ 9], Ar); Ar += LDA;
                    copy(&Uw[10], Ar); Ar += LDA;
                    copy(&Uw[11], Ar); Ar += LDA;
                    copy(&Uw[12], Ar); Ar += LDA;
                    copy(&Uw[13], Ar); Ar += LDA;
                    copy(&Uw[14], Ar); Ar += LDA;
                    copy(&Uw[15], Ar); Ar += LDA;
                    Uw += 16;
                } while ((UwEnd - Uw) >= 16);
                for (;Uw < UwEnd; ++Uw) {
                    copy(Uw, Ar); Ar += LDA;
                }
            } else {
                const size_t rowAw = -LINDXAU[M - 1];
                double *__restrict__ Aw = &A[rowAw];
                size_t col = range.begin();
                Aw += col * LDA;
                Ar += col * LDA;
                col = range.end() - col;
                do {
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                    col -= 16;
                } while (col >= 16);
                for (; col; --col) {
                    streamingCopy(Aw, Ar); Aw += LDA; Ar += LDA;
                }
            }
        }
};
#endif /* USE_ORIGINAL_LASWP */

extern "C" void HPL_dlaswp01T(const int M, const int N, double *A, const int LDA,
        double *U, const int LDU, const int *LINDXA, const int *LINDXAU)
{
#ifdef TRACE_CALLS
   uint64_t tr_start, tr_end, tr_diff;
   tr_start = util_getTimestamp();
   struct timespec cputime0, cputime1;
   clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cputime0);
#endif /* TRACE_CALLS */

#ifdef USE_ORIGINAL_LASWP
#include "HPL_dlaswp01T.c"
#else
    if(__builtin_expect(N < 16, 0)) {
        if (__builtin_expect(N <= 0, 0)) {
            return;
        }

        for (int i = 0; i < M; ++i) {
            const double *Ar = &A[LINDXA[i]];
            if (LINDXAU[i] >= 0) {
                double *Uw = &U[LINDXAU[i] * static_cast<size_t>(LDU)];
                for (int col = 0; col < N; ++col) {
                    Uw[col] = *Ar;
                    Ar += LDA;
                }
            } else {
                double *Aw = &A[-LINDXAU[i]];
                for (int col = 0; col < N; ++col) {
                    *Aw = *Ar;
                    Ar += LDA;
                    Aw += LDA;
                }
            }
        }
    } else {
        tbb::parallel_for (tbb::blocked_range<size_t>(0, N, 48),
                dlaswp01T_impl(M, A, LDA, U, LDU, LINDXA, LINDXAU),
                tbb::simple_partitioner());
    }
#endif

#ifdef TRACE_CALLS
   tr_end = util_getTimestamp();
   clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cputime1);
   const uint64_t cputime = cputime1.tv_sec * 10000000ull + cputime1.tv_nsec / 1000ull - cputime0.tv_sec * 10000000ull - cputime0.tv_nsec / 1000ull;
   tr_diff = util_getTimeDifference( tr_start, tr_end );

   fprintf( trace_dgemm, "DLASWP01T,M=%i,N=%i,LDA=%i,LDU=%i,TIME=%lu,CPU=%lu,THRPT=%.2fGB/s\n", M, N, LDA, LDU, tr_diff, cputime,
           0.002 * sizeof(double) * M * N / tr_diff );
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
