/* 
 * This is a modified version of the High Performance Computing Linpack
 * Benchmark (HPL). All code not contained in the original HPL version
 * 2.0 is property of the Frankfurt Institute for Advanced Studies
 * (FIAS). None of the material may be copied, reproduced, distributed,
 * republished, downloaded, displayed, posted or transmitted in any form
 * or by any means, including, but not limited to, electronic,
 * mechanical, photocopying, recording, or otherwise, without the prior
 * written permission of FIAS. For those parts contained in the
 * unmodified version of the HPL the below copyright notice applies.
 * 
 * Authors:
 * David Rohr (drohr@jwdt.org)
 * Matthias Bach (bach@compeng.uni-frankfurt.de)
 * Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 * 
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

#include "util_timer.h"
#include "util_trace.h"
#include "helpers.h"

extern "C" void HPL_dlaswp10N(const int M, const int N, double *A,
        const int LDA, const int *IPIV)
{
#ifdef TRACE_CALLS
   uint64_t tr_start, tr_end, tr_diff;
   tr_start = util_getTimestamp();
   int realN = 0;
#endif /* TRACE_CALLS */

   if( ( M <= 0 ) || ( N <= 0 ) ) return;

   const int mu = M & ~15u;
   const int mr = M - mu;

   for(int j = 0; j < N; j++ ) {
       const int jp = IPIV[j];
       if( j != jp ) {
#ifdef TRACE_CALLS
           ++realN;
#endif
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

#ifdef TRACE_CALLS
   tr_end = util_getTimestamp();
   tr_diff = util_getTimeDifference( tr_start, tr_end );

   fprintf( trace_dgemm, "DLASWP10N,M=%i,N=%i,LDA=%i,TIME=%lu,THRPT=%.2fGB/s\n", M, N, LDA, tr_diff,
           0.004 * sizeof(double) * M * realN / tr_diff);
#ifdef TRACE_PERMDATA
   char filename[256];
   snprintf(filename, 256, "dlaswp10N.%04d.%05d.%05d.dat", M, N, LDA);
   FILE *permdata = fopen(filename, "w");
   fwrite(IPIV, sizeof(IPIV[0]), N, permdata);
   fclose(permdata);
#endif
#endif /* TRACE_CALLS */
}
