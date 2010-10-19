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
   const size_t N = _N;
   const size_t LDA = _LDA;
   const size_t LDB = _LDB;

   // B_ij = A_ji

   size_t i = 0;
   for ( ; i < M - 7; i += 8 )
   {
      double *__restrict__ B_ij = &B[ i ];
      const double *__restrict__ A_ji = &A[ i * LDA ];
      size_t j = 0;
      for ( ; j < N - 7; j += DoublesInCacheline )
      {
         _mm_prefetch( &A_ji[ 8 + 0 * LDA ], _MM_HINT_NTA );
         _mm_prefetch( &A_ji[ 8 + 1 * LDA ], _MM_HINT_NTA );
         _mm_prefetch( &A_ji[ 8 + 2 * LDA ], _MM_HINT_NTA );
         _mm_prefetch( &A_ji[ 8 + 3 * LDA ], _MM_HINT_NTA );
         _mm_prefetch( &A_ji[ 8 + 4 * LDA ], _MM_HINT_NTA );
         _mm_prefetch( &A_ji[ 8 + 5 * LDA ], _MM_HINT_NTA );
         _mm_prefetch( &A_ji[ 8 + 6 * LDA ], _MM_HINT_NTA );
         _mm_prefetch( &A_ji[ 8 + 7 * LDA ], _MM_HINT_NTA );
         for ( size_t j2 = 0; j2 < DoublesInCacheline; ++j2 )
         {
            // collect one cacheline in a store buffer
            // reading 8 streams from A linearly
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
   for ( size_t j = 0; j < N; ++j )
   {
      for ( ; i < M; ++i )
      {
         streamingCopy( &B[ i + j * LDB ], &A[ j + i * LDA ] );
      }
   }

   // make sure the streaming stores are visible to subsequent loads and stores to B
   _mm_mfence();

#if 0
   const double               * A0 = A,              * A1 = A + 1;
   double                     * B0 = B,              * B1 = B +     LDB;
   const int                  incA = -M * LDA + (1 << 1),
                              incB = ( (unsigned int)(LDB) << 1 ) - M,
                              incA0 = -M * LDA + 1, incB0 = LDB - M;
   int                        mu, nu;


   mu = (int)( ( (unsigned int)(M) >> 2 ) << 2 );
   nu = (int)( ( (unsigned int)(N) >> 1 ) << 1 );

   for( int j = 0; j < nu; j += 2 ) {
      for( int i = 0; i < mu; i += 4 ) {
         B0[ 0] = *A0; A0 += LDA; B1[ 0] = *A1; A1 += LDA;
         B0[ 1] = *A0; A0 += LDA; B1[ 1] = *A1; A1 += LDA;
         B0[ 2] = *A0; A0 += LDA; B1[ 2] = *A1; A1 += LDA;
         B0[ 3] = *A0; A0 += LDA; B1[ 3] = *A1; A1 += LDA;
         B0 += 4; B1 += 4;
      }

      for( int i = mu; i < M; i++ ) {
         *B0 = *A0; B0++; A0 += LDA; *B1 = *A1; B1++; A1 += LDA;
      }

      A0 += incA; A1 += incA; B0 += incB; B1 += incB;
   }

   for( int j = nu; j < N; j++, B0 += incB0, A0 += incA0 ) {
      for( int i = 0; i < mu; i += 4, B0 += 4 ) {
         B0[ 0]=*A0; A0 += LDA;
         B0[ 1]=*A0; A0 += LDA;
         B0[ 2]=*A0; A0 += LDA; B0[ 3]=*A0; A0 += LDA;
      }

      for( int i = mu; i < M; i++, B0++, A0 += LDA ) {
         *B0 = *A0;
      }
   }
#endif

   tsc.stop();
   END_TRACE
   fprintf( stderr, "=> %f GB/s\n", M * N * sizeof( double ) / tsc.cycles() );
}
