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
 *	HPL - 2.0 - September 10, 2008						  
 *	Antoine P. Petitet												
 *	University of Tennessee, Knoxville								
 *	Innovative Computing Laboratory								 
 *	(C) Copyright 2000-2008 All Rights Reserved					   
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

#include "../pauxil/helpers.h"
#include <tbb/parallel_for.h>

typedef MyRange<8, 16> Range;

inline void dlacpy_worker(const double* __restrict__ A, double* __restrict__ B, size_t N, size_t begin, size_t end, size_t LDA, size_t LDB)
{
	for (size_t j = 0 ; j < N; j ++ )
	{
		double *__restrict__ B_ij = &B[ j * LDB ];
		const double *__restrict__ A_ij = &A[ j * LDA ];
		for (size_t i = begin; i < end; i += 8)
		{
			//_mm_prefetch( &A_ij[i + LDA], _MM_HINT_NTA);
			_mm_stream_pd( &B_ij[i], _mm_load_pd( &A_ij[i]));
			_mm_stream_pd( &B_ij[i + 2], _mm_load_pd( &A_ij[i + 2]));
			_mm_stream_pd( &B_ij[i + 4], _mm_load_pd( &A_ij[i + 4]));
			_mm_stream_pd( &B_ij[i + 6], _mm_load_pd( &A_ij[i + 6]));
		}
	}
}

class HPL_dlacpy_impl
{
	private:
		const size_t N;
		const size_t LDA;
		const size_t LDB;
		const double *__restrict__ const A;
		double *__restrict__ const B;

	public:
		HPL_dlacpy_impl(size_t _N, const double *_A, size_t _LDA, double *_B, size_t _LDB)
			: N(_N), LDA(_LDA), LDB(_LDB), A(_A), B(_B)
		{
		}

		void operator()(const Range &range) const
		{
			const size_t end = range.begin() + range.N();
			dlacpy_worker(A, B, N, range.begin(), end, LDA, LDB);
		}
};

/* 
 * Purpose
 * =======
 *
 * HPL_dlacpy copies an array A into an array B.
 * 
 *
 * Arguments
 * =========
 *
 * M	   (local input)				 const int
 *		 On entry,  M specifies the number of rows of the arrays A and
 *		 B. M must be at least zero.
 *
 * N	   (local input)				 const int
 *		 On entry,  N specifies  the number of columns of the arrays A
 *		 and B. N must be at least zero.
 *
 * A	   (local input)				 const double *
 *		 On entry, A points to an array of dimension (LDA,N).
 *
 * LDA	 (local input)				 const int
 *		 On entry, LDA specifies the leading dimension of the array A.
 *		 LDA must be at least MAX(1,M).
 *
 * B	   (local output)				double *
 *		 On entry, B points to an array of dimension (LDB,N). On exit,
 *		 B is overwritten with A.
 *
 * LDB	 (local input)				 const int
 *		 On entry, LDB specifies the leading dimension of the array B.
 *		 LDB must be at least MAX(1,M).
 *
 * ---------------------------------------------------------------------
 */
extern "C" void HPL_dlacpy(const int _M, const int _N, const double *A, const int _LDA, double *B, const int _LDB, int multithread)
{
   START_TRACE( DLACPY )

   if ( _M <= 0 || _N <= 0 ) {
	  return;
   }

   const size_t M = _M;
   const size_t MM = M & ~7;
   const size_t N = _N;
   const size_t LDA = _LDA;
   const size_t LDB = _LDB;

   // B_ij = A_ji

   if (multithread)
   {
    tbb::parallel_for( Range(0, MM), HPL_dlacpy_impl( N, A, LDA, B, LDB ), tbb::auto_partitioner() );
   }
   else
   {
    dlacpy_worker(A, B, N, 0, MM, LDA, LDB);
   }

   if ( M & 7 )
   {
	  for ( size_t j = 0; j < N; ++j )
	  {
		 for ( size_t i = MM; i < M; ++i )
		 {
			streamingCopy( &B[ i + j * LDB ], &A[ i + j * LDA ] );
		 }
	  }
   }

   // make sure the streaming stores are visible to subsequent loads
   // and stores to B
   _mm_mfence();

   END_TRACE
}
