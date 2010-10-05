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

extern "C" void HPL_dlaswp10N(const int M, const int N, double *A,
        const int LDA, const int *IPIV)
{
START_TRACE( HPL_DLASWP10N )

#ifdef USE_ORIGINAL_LASWP
#include "HPL_dlaswp10N.c"
#else
   if (M <= 0) {
       return;
   }

   const int mu = M & ~15u;
   const int mr = M - mu;

   for(int j = 0; j < N; j++ ) {
       const int jp = IPIV[j];
       if( j != jp ) {
           double *__restrict__ a0 = A + j * LDA;
           double *__restrict__ a1 = A + jp * LDA;

           _mm_prefetch(&a0[ 0], _MM_HINT_NTA);
           _mm_prefetch(&a1[ 0], _MM_HINT_NTA);
           _mm_prefetch(&a0[ 8], _MM_HINT_NTA);
           _mm_prefetch(&a1[ 8], _MM_HINT_NTA);
           _mm_prefetch(&a0[16], _MM_HINT_NTA);
           _mm_prefetch(&a1[16], _MM_HINT_NTA);
           _mm_prefetch(&a0[24], _MM_HINT_NTA);
           _mm_prefetch(&a1[24], _MM_HINT_NTA);
           if ((a0 - static_cast<double *>(0)) & 1) {
               swap(a0[0], a1[0]);
               ++a0;
               ++a1;
           }
           // If LDA is odd life sucks
           for(int i = 0; i < mu; i += 16, a0 += 16, a1 += 16 ) {
               _mm_prefetch(&a0[32], _MM_HINT_NTA);
               _mm_prefetch(&a1[32], _MM_HINT_NTA);
               swapSSE(a0[ 0], a1[ 0]);
               swapSSE(a0[ 2], a1[ 2]);
               swapSSE(a0[ 4], a1[ 4]);
               swapSSE(a0[ 6], a1[ 6]);
               _mm_prefetch(&a0[40], _MM_HINT_NTA);
               _mm_prefetch(&a1[40], _MM_HINT_NTA);
               swapSSE(a0[ 8], a1[ 8]);
               swapSSE(a0[10], a1[10]);
               swapSSE(a0[12], a1[12]);
               swapSSE(a0[14], a1[14]);
           }

           for(int i = 0; i < mr; i++ ) {
               swap(a0[i], a1[i]);
           }
       }
   }
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
