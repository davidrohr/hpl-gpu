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
/*
 * Include files
 */
#include "util_trace.h"

#include "tsc.h"
#include "../pauxil/helpers.h"
#include <tbb/parallel_for.h>

typedef MyRange<8, 16> Range;

class HPL_dlatcpy_impl
{
    private:
        const size_t N;
        const size_t LDA;
        const size_t LDB;
        const double *__restrict__ const A;
        double *__restrict__ const B;

    public:
        HPL_dlatcpy_impl(size_t _N, const double *_A, size_t _LDA,
                double *_B, size_t _LDB)
            : N(_N), LDA(_LDA), LDB(_LDB),
            A(_A), B(_B)
        {
        }

        void operator()(const Range &range) const
        {
            const size_t end = range.begin() + range.N();
            for ( size_t i = range.begin(); i < end; i += 8 )
            {
                double *__restrict__ B_ij = &B[ i ];
                const double *__restrict__ A_ji = &A[ i * LDA ];
                size_t j = 0;
                for ( ; j < (N & ~7); j += 8 )
                {
                    _mm_prefetch( &A_ji[ 8 + 0 * LDA ], _MM_HINT_NTA );
                    _mm_prefetch( &A_ji[ 8 + 1 * LDA ], _MM_HINT_NTA );
#if 0
#define prefetchw(addr) _m_prefetchw(addr)
#else
#define prefetchw(addr) _mm_prefetch(addr, _MM_HINT_NTA )
#endif
                    prefetchw( &B_ij[ 8 + 0 * LDB ] );
                    prefetchw( &B_ij[ 8 + 1 * LDB ] );
                    prefetchw( &B_ij[ 8 + 2 * LDB ] );
                    prefetchw( &B_ij[ 8 + 3 * LDB ] );
                    prefetchw( &B_ij[ 8 + 4 * LDB ] );
                    prefetchw( &B_ij[ 8 + 5 * LDB ] );
                    prefetchw( &B_ij[ 8 + 6 * LDB ] );
                    prefetchw( &B_ij[ 8 + 7 * LDB ] );
#undef prefetchw

                    _mm_prefetch( &A_ji[ 8 + 2 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 3 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 4 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 5 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 6 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 7 * LDA ], _MM_HINT_T1 );
                    for ( size_t i2 = 0; i2 < 8; i2 += 2 )
                    {
                       const __m128d tmp0 = _mm_load_pd( &A_ji[ (i2 + 0) * LDA     ] );
                       const __m128d tmp1 = _mm_load_pd( &A_ji[ (i2 + 1) * LDA     ] );
                       const __m128d tmp2 = _mm_load_pd( &A_ji[ (i2 + 0) * LDA + 2 ] );
                       const __m128d tmp3 = _mm_load_pd( &A_ji[ (i2 + 1) * LDA + 2 ] );
                       const __m128d tmp4 = _mm_load_pd( &A_ji[ (i2 + 0) * LDA + 4 ] );
                       const __m128d tmp5 = _mm_load_pd( &A_ji[ (i2 + 1) * LDA + 4 ] );
                       const __m128d tmp6 = _mm_load_pd( &A_ji[ (i2 + 0) * LDA + 6 ] );
                       const __m128d tmp7 = _mm_load_pd( &A_ji[ (i2 + 1) * LDA + 6 ] );
                       _mm_prefetch( &A_ji[ (i2 + 2) * LDA ], _MM_HINT_T0 );
                       _mm_prefetch( &A_ji[ (i2 + 3) * LDA ], _MM_HINT_T0 );
                       _mm_store_pd( &B_ij[ i2           ], _mm_unpacklo_pd( tmp0, tmp1 ) );
                       _mm_store_pd( &B_ij[ i2 + 1 * LDB ], _mm_unpackhi_pd( tmp0, tmp1 ) );
                       _mm_store_pd( &B_ij[ i2 + 2 * LDB ], _mm_unpacklo_pd( tmp2, tmp3 ) );
                       _mm_store_pd( &B_ij[ i2 + 3 * LDB ], _mm_unpackhi_pd( tmp2, tmp3 ) );
                       _mm_store_pd( &B_ij[ i2 + 4 * LDB ], _mm_unpacklo_pd( tmp4, tmp5 ) );
                       _mm_store_pd( &B_ij[ i2 + 5 * LDB ], _mm_unpackhi_pd( tmp4, tmp5 ) );
                       _mm_store_pd( &B_ij[ i2 + 6 * LDB ], _mm_unpacklo_pd( tmp6, tmp7 ) );
                       _mm_store_pd( &B_ij[ i2 + 7 * LDB ], _mm_unpackhi_pd( tmp6, tmp7 ) );
                    }
                    B_ij += 8 * LDB;
                    A_ji += 8;
                }
                for ( ; j < N; ++j )
                {
                   streamingCopy( &B_ij[ 0 ], &A_ji[ 0 * LDA ] );
                   streamingCopy( &B_ij[ 1 ], &A_ji[ 1 * LDA ] );
                   streamingCopy( &B_ij[ 2 ], &A_ji[ 2 * LDA ] );
                   streamingCopy( &B_ij[ 3 ], &A_ji[ 3 * LDA ] );
                   streamingCopy( &B_ij[ 4 ], &A_ji[ 4 * LDA ] );
                   streamingCopy( &B_ij[ 5 ], &A_ji[ 5 * LDA ] );
                   streamingCopy( &B_ij[ 6 ], &A_ji[ 6 * LDA ] );
                   streamingCopy( &B_ij[ 7 ], &A_ji[ 7 * LDA ] );

                   B_ij += LDB;
                   ++A_ji;
                }
            }
        }
};

class HPL_dlatcpy_impl2
{
    private:
        const size_t N;
        const size_t LDA;
        const size_t LDB;
        const double *__restrict__ const A;
        double *__restrict__ const B;

    public:
        HPL_dlatcpy_impl2(size_t _N, const double *_A, size_t _LDA,
                double *_B, size_t _LDB)
            : N(_N), LDA(_LDA), LDB(_LDB),
            A(_A), B(_B)
        {
        }

        void operator()(const Range &range) const
        {
            const size_t end = range.begin() + range.N();
            for (size_t i = range.begin(); i < end; i += 8) {
                double *__restrict__ B_ij = &B[ i ];
                const double *__restrict__ A_ji = &A[ i * LDA ];
                size_t j = 0;
                for ( ; j < (N & ~7); j += 8 )
                {
                    _mm_prefetch( &A_ji[ 8 + 0 * LDA ], _MM_HINT_NTA );
                    _mm_prefetch( &A_ji[ 8 + 1 * LDA ], _MM_HINT_NTA );
                    _mm_prefetch( &A_ji[ 8 + 2 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 3 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 4 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 5 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 6 * LDA ], _MM_HINT_T1 );
                    _mm_prefetch( &A_ji[ 8 + 7 * LDA ], _MM_HINT_T1 );
                    for ( size_t j2 = 0; j2 < 8; j2 += 2 )
                    {
                        for ( size_t i2 = 0; i2 < 8; i2 += 2 )
                        {
                            const __m128d tmp0 = _mm_load_pd( &A_ji[ i2 * LDA ] );
                            const __m128d tmp1 = _mm_load_pd( &A_ji[ i2 * LDA + LDA ] );
                            _mm_stream_pd( &B_ij[ i2       ], _mm_unpacklo_pd( tmp0, tmp1 ) );
                            _mm_stream_pd( &B_ij[ i2 + LDB ], _mm_unpackhi_pd( tmp0, tmp1 ) );
                        }
                        B_ij += 2 * LDB;
                        A_ji += 2;
                    }
                }
                for ( ; j < N; ++j )
                {
                    for ( size_t i2 = 0; i2 < 8; ++i2 )
                    {
                        streamingCopy( &B_ij[ i2 ], &A_ji[ i2 * LDA ] );
                    }

                    B_ij += LDB;
                    ++A_ji;
                }
            }
        }
};

/**
 * Purpose
 * =======
 *
 * HPL_dlatcpy copies the transpose of an array A into an array B.
 * 
 *
 * Arguments
 * =========
 *
 * M       (local input)                 const int
 *         On entry,  M specifies the number of  rows of the array B and
 *         the number of columns of A. M must be at least zero.
 *
 * N       (local input)                 const int
 *         On entry,  N specifies the number of  rows of the array A and
 *         the number of columns of B. N must be at least zero.
 *
 * A       (local input)                 const double *
 *         On entry, A points to an array of dimension (LDA,M).
 *
 * LDA     (local input)                 const int
 *         On entry, LDA specifies the leading dimension of the array A.
 *         LDA must be at least MAX(1,N).
 *
 * B       (local output)                double *
 *         On entry, B points to an array of dimension (LDB,N). On exit,
 *         B is overwritten with the transpose of A.
 *
 * LDB     (local input)                 const int
 *         On entry, LDB specifies the leading dimension of the array B.
 *         LDB must be at least MAX(1,M).
 *
 * ---------------------------------------------------------------------
 */
extern "C" void HPL_dlatcpy(const int _M, const int _N, const double *A, const int _LDA, double *B, const int _LDB)
{
   TimeStampCounter tsc;
   START_TRACE( DLATCPY )
   tsc.start();

   if ( _M <= 0 || _N <= 0 ) {
      return;
   }

   const size_t M = _M;
   const size_t MM = M & ~7;
   const size_t N = _N;
   const size_t LDA = _LDA;
   const size_t LDB = _LDB;

   // B_ij = A_ji

   size_t BS = 16;
   if (N < 1024) {
      // BS * N ~ 512kB ~ 60k doubles => BS ~ 60000 / N
      BS = 60000 / N;
      BS = BS > MM ? BS / 2 : BS;
      tbb::parallel_for( Range(0, MM, BS), HPL_dlatcpy_impl( N, A, LDA, B, LDB ),
            tbb::simple_partitioner() );
   } else {
      tbb::parallel_for( Range(0, MM), HPL_dlatcpy_impl2( N, A, LDA, B, LDB ),
            tbb::simple_partitioner() );
   }

   if ( M & 7 )
   {
      for ( size_t j = 0; j < N; ++j )
      {
         for ( size_t i = MM; i < M; ++i )
         {
            streamingCopy( &B[ i + j * LDB ], &A[ j + i * LDA ] );
         }
      }
   }

   // make sure the streaming stores are visible to subsequent loads
   // and stores to B
   _mm_mfence();

   tsc.stop();
   END_TRACE
   fprintf( stderr, "%s(M = %4d, N = %5d, LDA = %5d, LDB = %5d) => Bs = %4ld, %.2f GB/s\n", __func__,
         _M, _N, _LDA, _LDB, BS, 2.1 * 2 * M * N * sizeof( double ) / tsc.cycles() );
}
