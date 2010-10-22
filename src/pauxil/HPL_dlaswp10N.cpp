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
                        swapSSE(a0[ 0], a1[ 0]);
                        swapSSE(a0[ 2], a1[ 2]);
                        swapSSE(a0[ 4], a1[ 4]);
                        swapSSE(a0[ 6], a1[ 6]);
                        _m_prefetchw(&a0[40]);
                        _m_prefetchw(&a1[40]);
                        swapSSE(a0[ 8], a1[ 8]);
                        swapSSE(a0[10], a1[10]);
                        swapSSE(a0[12], a1[12]);
                        swapSSE(a0[14], a1[14]);
                    }
                }
            }
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
                for (j = 0; j < (M & 14); j += 2) {
                    swapSSE(a0[j], a1[j]);
                }
                if (M & 1) {
                    swap(a0[j], a1[j]);
                }
            }
        }
    }

    tbb::parallel_for (Range10N(0, M & ~15), dlaswp10N_impl(N, A, LDA, IPIV));
#endif

END_TRACE
#ifdef TRACE_CALLS
#ifdef TRACE_PERMDATA
   char filename[256];
   snprintf(filename, 256, "dlaswp10N.%04d.%05d.%05d.dat", M, N, LDA);
   FILE *permdata = fopen(filename, "w");
   fwrite(IPIV, sizeof(IPIV[0]), N, permdata);
   fclose(permdata);
#endif
#endif /* TRACE_CALLS */
}
