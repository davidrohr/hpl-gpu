/*
 *  -- High Performance Computing Linpack Benchmark (HPL-GPU)
 *     HPL-GPU - 1.1 - 2011
 *
 *     David Rohr
 *     Matthias Kretz
 *     Matthias Bach
 *     Goethe Universität, Frankfurt am Main
 *     Frankfurt Institute for Advanced Studies
 *     (C) Copyright 2010 All Rights Reserved
 *
 *     Antoine P. Petitet
 *     University of Tennessee, Knoxville
 *     Innovative Computing Laboratory
 *     (C) Copyright 2000-2008 All Rights Reserved
 *
 *  -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 *  1. Redistributions  of  source  code  must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce  the above copyright
 *  notice, this list of conditions,  and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. All  advertising  materials  mentioning  features  or  use of this
 *  software must display the following acknowledgements:
 *  This  product  includes  software  developed  at  the  University  of
 *  Tennessee, Knoxville, Innovative Computing Laboratory.
 *  This product  includes software  developed at the Frankfurt Institute
 *  for Advanced Studies.
 *
 *  4. The name of the  University,  the name of the  Laboratory,  or the
 *  names  of  its  contributors  may  not  be used to endorse or promote
 *  products  derived   from   this  software  without  specific  written
 *  permission.
 *
 *  -- Disclaimer:
 *
 *  THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 *  OR  CONTRIBUTORS  BE  LIABLE FOR ANY  DIRECT,  INDIRECT,  INCIDENTAL,
 *  SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES  (INCLUDING,  BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,  OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * ======================================================================
 */

/*
 * Define default value for unrolling factor
 */
#ifndef HPL_LASWP01T_DEPTH
#define    HPL_LASWP01T_DEPTH       32
#define    HPL_LASWP01T_LOG2_DEPTH   5
#endif

/*
 * .. Local Variables ..
 */
   double                     * a0, * a1;
   const int                  incA = (int)( (unsigned int)(LDA) <<
                                            HPL_LASWP01T_LOG2_DEPTH ),
                              incU = ( 1 << HPL_LASWP01T_LOG2_DEPTH );
   int                        nu, nr;
   register int               i, j;
/* ..
 * .. Executable Statements ..
 */
   if( ( M <= 0 ) || ( N <= 0 ) ) return;

   nr = N - ( nu = (int)( ( (unsigned int)(N) >> HPL_LASWP01T_LOG2_DEPTH ) <<
                            HPL_LASWP01T_LOG2_DEPTH ) );

   for( j = 0; j < nu; j += HPL_LASWP01T_DEPTH, A += incA, U += incU )
   {
      for( i = 0; i < M; i++ )
      {
         a0 = A + (size_t)(LINDXA[i]);

         if( LINDXAU[i] >= 0 )
         {
            a1 = U + (size_t)(LINDXAU[i]) * (size_t)(LDU);

            a1[ 0] = *a0; a0 += LDA;
#if ( HPL_LASWP01T_DEPTH >  1 )
            a1[ 1] = *a0; a0 += LDA;
#endif
#if ( HPL_LASWP01T_DEPTH >  2 )
            a1[ 2] = *a0; a0 += LDA; a1[ 3] = *a0; a0 += LDA;
#endif
#if ( HPL_LASWP01T_DEPTH >  4 )
            a1[ 4] = *a0; a0 += LDA; a1[ 5] = *a0; a0 += LDA;
            a1[ 6] = *a0; a0 += LDA; a1[ 7] = *a0; a0 += LDA;
#endif
#if ( HPL_LASWP01T_DEPTH >  8 )
            a1[ 8] = *a0; a0 += LDA; a1[ 9] = *a0; a0 += LDA;
            a1[10] = *a0; a0 += LDA; a1[11] = *a0; a0 += LDA;
            a1[12] = *a0; a0 += LDA; a1[13] = *a0; a0 += LDA;
            a1[14] = *a0; a0 += LDA; a1[15] = *a0; a0 += LDA;
#endif
#if ( HPL_LASWP01T_DEPTH > 16 )
            a1[16] = *a0; a0 += LDA; a1[17] = *a0; a0 += LDA;
            a1[18] = *a0; a0 += LDA; a1[19] = *a0; a0 += LDA;
            a1[20] = *a0; a0 += LDA; a1[21] = *a0; a0 += LDA;
            a1[22] = *a0; a0 += LDA; a1[23] = *a0; a0 += LDA;
            a1[24] = *a0; a0 += LDA; a1[25] = *a0; a0 += LDA;
            a1[26] = *a0; a0 += LDA; a1[27] = *a0; a0 += LDA;
            a1[28] = *a0; a0 += LDA; a1[29] = *a0; a0 += LDA;
            a1[30] = *a0; a0 += LDA; a1[31] = *a0; a0 += LDA;
#endif
         }
         else
         {
            a1 = A - (size_t)(LINDXAU[i]);

            *a1 = *a0; a1 += LDA; a0 += LDA;
#if ( HPL_LASWP01T_DEPTH >  1 )
            *a1 = *a0; a1 += LDA; a0 += LDA;
#endif
#if ( HPL_LASWP01T_DEPTH >  2 )
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
#endif
#if ( HPL_LASWP01T_DEPTH >  4 )
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
#endif
#if ( HPL_LASWP01T_DEPTH >  8 )
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
#endif
#if ( HPL_LASWP01T_DEPTH > 16 )
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
            *a1 = *a0; a1 += LDA; a0 += LDA; *a1 = *a0; a1 += LDA; a0 += LDA;
#endif
         }
      }
   }

   if( nr > 0 )
   {
      for( i = 0; i < M; i++ )
      {
         a0 = A + (size_t)(LINDXA[i]);

         if( LINDXAU[i] >= 0 )
         {
            a1 = U + (size_t)(LINDXAU[i]) * (size_t)(LDU);
            for( j = 0; j < nr; j++, a0 += LDA ) { a1[j] = *a0; }
         }
         else
         {
            a1 = A - (size_t)(LINDXAU[i]);
            for( j = 0; j < nr; j++, a1 += LDA, a0 += LDA ) { *a1 = *a0; }
         }
      }
   }
