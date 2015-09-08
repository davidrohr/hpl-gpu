/*
 *  -- High Performance Computing Linpack Benchmark (HPL-GPU)
 *     HPL-GPU - 1.0 - 2010
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
 * Include files
 */
#ifndef STD_OUT
#define STD_OUT stdout
#endif
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
#ifdef HPL_LASWP_AVX
			_mm256_stream_pd( &B_ij[i + 0], _mm256_load_pd( &A_ij[i + 0]));
			_mm256_stream_pd( &B_ij[i + 4], _mm256_load_pd( &A_ij[i + 4]));
#else
			_mm_stream_pd( &B_ij[i + 0], _mm_load_pd( &A_ij[i + 0]));
			_mm_stream_pd( &B_ij[i + 2], _mm_load_pd( &A_ij[i + 2]));
			_mm_stream_pd( &B_ij[i + 4], _mm_load_pd( &A_ij[i + 4]));
			_mm_stream_pd( &B_ij[i + 6], _mm_load_pd( &A_ij[i + 6]));
#endif
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
   
   if ((_LDA & 1) || (_LDB & 1))
   {
      fprintf(STD_OUT, "ERROR: Uneven leading dimension not supported by dlacpy\n");
      exit(1);
   }

   const size_t M = _M;
   size_t MM = M & ~7;
   const size_t N = _N;
   const size_t LDA = _LDA;
   const size_t LDB = _LDB;

   // B_ij = A_ji
   
   if (LDA & 1 || LDB & 1 || ((size_t) A) & 15 || ((size_t) B) & 15)
   {
    MM = 0;
    fprintf(STD_OUT, "WARNING, parameters unaligned, using slow copy routing LDA %lld LDB %lld A %lld B %lld\n", (long long int) LDA, (long long int) LDB, (long long int) (size_t) A, (long long int) (size_t) B);
    goto Unaligned;
   }
   else if (multithread)
   {
    tbb::parallel_for( Range(0, MM), HPL_dlacpy_impl( N, A, LDA, B, LDB ), tbb::auto_partitioner() );
   }
   else
   {
    dlacpy_worker(A, B, N, 0, MM, LDA, LDB);
   }
   if ( M & 7 )
   {
Unaligned:
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
#ifdef TRACE_LASWP
   char filename[256];
   snprintf(filename, 256, "dlacpy.%04d.%05d.%05d.%05d.%7.4fs.dat", M, N, LDA, LDB, laswp_time);
   FILE *permdata = fopen(filename, "w");
   fclose(permdata);
#endif
}
