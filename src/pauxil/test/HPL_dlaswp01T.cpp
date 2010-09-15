#include "../HPL_dlaswp00N.cpp"
#include "../permutationhelper.h"
#include "../permutationhelper.cpp"
#include "matrix.h"
#include <mm3dnow.h>
#include <emmintrin.h>
#include <tbb/parallel_for.h>

#ifndef HPL_LASWP01T_LOG2_DEPTH
#define    HPL_LASWP01T_LOG2_DEPTH   5
#endif
#ifndef HPL_LASWP01T_DEPTH
#define    HPL_LASWP01T_DEPTH        (1 << HPL_LASWP01T_LOG2_DEPTH)
#endif

extern "C" void HPL_dlaswp01T
(
   const int                        M,
   const int                        N,
   double *                         A,
   const int                        LDA,
   double *                         U,
   const int                        LDU,
   const int *                      LINDXA,
   const int *                      LINDXAU
)
{
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
/*
 * End of HPL_dlaswp01T
 */

/*
 * Purpose
 * =======
 *
 * HPL_dlaswp01T copies  scattered rows  of  A  into itself  and into an
 * array U.  The row offsets in  A  of the source rows  are specified by
 * LINDXA.  The  destination of those rows are specified by  LINDXAU.  A
 * positive value of LINDXAU indicates that the array  destination is U,
 * and A otherwise. Rows of A are stored as columns in U.
 *
 * Arguments
 * =========
 *
 * M       (local input)                 const int
 *         On entry, M  specifies the number of rows of A that should be
 *         moved within A or copied into U. M must be at least zero.
 *
 * N       (local input)                 const int
 *         On entry, N  specifies the length of rows of A that should be
 *         moved within A or copied into U. N must be at least zero.
 *
 * A       (local input/output)          double *
 *         On entry, A points to an array of dimension (LDA,N). The rows
 *         of this array specified by LINDXA should be moved within A or
 *         copied into U.
 *
 * LDA     (local input)                 const int
 *         On entry, LDA specifies the leading dimension of the array A.
 *         LDA must be at least MAX(1,M).
 *
 * U       (local input/output)          double *
 *         On entry, U points to an array of dimension (LDU,M). The rows
 *         of A specified by  LINDXA  are copied within this array  U at
 *         the  positions indicated by positive values of LINDXAU.  The
 *         rows of A are stored as columns in U.
 *
 * LDU     (local input)                 const int
 *         On entry, LDU specifies the leading dimension of the array U.
 *         LDU must be at least MAX(1,N).
 *
 * LINDXA  (local input)                 const int *
 *         On entry, LINDXA is an array of dimension M that contains the
 *         local  row indexes  of  A  that should be moved within  A  or
 *         or copied into U.
 *
 * LINDXAU (local input)                 const int *
 *         On entry, LINDXAU  is an array of dimension  M that  contains
 *         the local  row indexes of  U  where the rows of  A  should be
 *         copied at. This array also contains the  local row offsets in
 *         A where some of the rows of A should be moved to.  A positive
 *         value of  LINDXAU[i]  indicates that the row  LINDXA[i]  of A
 *         should be copied into U at the position LINDXAU[i]; otherwise
 *         the row  LINDXA[i]  of  A  should be moved  at  the  position
 *         -LINDXAU[i] within A.
 *
 * ---------------------------------------------------------------------
 */

class dlaswp01T_task_U : public tbb::task
{
    private:
    public:
};

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
//X             if (range.end() - range.begin() < 16) {
//X                 fprintf(stderr, "range too small: %ld\n", range.end() - range.begin());
//X                 abort();
//X             }
//X             if (range.end() - range.begin() > 32) {
//X                 fprintf(stderr, "range too large: %ld\n", range.end() - range.begin());
//X                 abort();
//X             }
            for (size_t i = 0; i < M - 1; ++i) {
                const double *Ar = &A[LINDXA[i]];
                const ptrdiff_t ArNext = &A[LINDXA[i + 1]] - Ar;
                if (LINDXAU[i] >= 0) {
                    const size_t rowUw = LINDXAU[i];
                    double *Uw = &U[rowUw * LDU];
                    double *const UwEnd = Uw + range.end();
                    Uw += range.begin();
                    _m_prefetchw(Uw);
                    _m_prefetchw(Uw + 8);
                    _m_prefetchw(Uw + 16);
                    Ar += range.begin() * LDA;
                    do {
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 0] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 1] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 2] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 3] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 4] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 5] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 6] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 7] = *Ar; Ar += LDA;
                        _mm_clflush(&Uw[0]);
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 8] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 9] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[10] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[11] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[12] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[13] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[14] = *Ar; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[15] = *Ar; Ar += LDA;
                        _mm_clflush(&Uw[8]);
                        Uw += 16;
                    } while ((UwEnd - Uw) >= 16);
                    for (;Uw < UwEnd; ++Uw) {
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); Uw[ 0] = *Ar; Ar += LDA;
                    }
                    _mm_clflush(&Uw[-1]);
                } else {
                    const size_t rowAw = -LINDXAU[i];
                    double *Aw = &A[rowAw];
                    size_t col = range.begin();
                    Aw += col * LDA;
                    Ar += col * LDA;
                    col = range.end() - col;
                    do {
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                        col -= 16;
                    } while (col >= 16);
                    for (; col; --col) {
                        _mm_prefetch(Ar + ArNext, _MM_HINT_T0); *Aw = *Ar; Aw += LDA; Ar += LDA;
                    }
                }
            }
            const double *Ar = &A[LINDXA[M - 1]];
            if (LINDXAU[M - 1] >= 0) {
                const size_t rowUw = LINDXAU[M - 1];
                double *Uw = &U[rowUw * LDU];
                double *const UwEnd = Uw + range.end();
                Uw += range.begin();
                _m_prefetchw(Uw);
                _m_prefetchw(Uw + 8);
                _m_prefetchw(Uw + 16);
                Ar += range.begin() * LDA;
                do {
                    Uw[ 0] = *Ar; Ar += LDA;
                    Uw[ 1] = *Ar; Ar += LDA;
                    Uw[ 2] = *Ar; Ar += LDA;
                    Uw[ 3] = *Ar; Ar += LDA;
                    Uw[ 4] = *Ar; Ar += LDA;
                    Uw[ 5] = *Ar; Ar += LDA;
                    Uw[ 6] = *Ar; Ar += LDA;
                    Uw[ 7] = *Ar; Ar += LDA;
                    _mm_clflush(&Uw[0]);
                    Uw[ 8] = *Ar; Ar += LDA;
                    Uw[ 9] = *Ar; Ar += LDA;
                    Uw[10] = *Ar; Ar += LDA;
                    Uw[11] = *Ar; Ar += LDA;
                    Uw[12] = *Ar; Ar += LDA;
                    Uw[13] = *Ar; Ar += LDA;
                    Uw[14] = *Ar; Ar += LDA;
                    Uw[15] = *Ar; Ar += LDA;
                    _mm_clflush(&Uw[8]);
                    Uw += 16;
                } while ((UwEnd - Uw) >= 16);
                for (;Uw < UwEnd; ++Uw) {
                    Uw[ 0] = *Ar; Ar += LDA;
                }
                _mm_clflush(&Uw[-1]);
            } else {
                const size_t rowAw = -LINDXAU[M - 1];
                double *Aw = &A[rowAw];
                size_t col = range.begin();
                Aw += col * LDA;
                Ar += col * LDA;
                for (; col < range.end(); ++col) {
                    *Aw = *Ar;
                    Aw += LDA;
                    Ar += LDA;
                }
            }
        }
};
extern "C" void dlaswp01T(const int M, const int N, double *A, const unsigned int LDA,
        double *U, const unsigned int LDU, const int *LINDXA, const int *LINDXAU)
{
    if( N <= 0 ) {
        return;
    }

#if 0
    // The idea is to create a list of tasks that operate on an area of memory such that no cache
    // line (64B == 8 doubles) which is written will be read/written by another task. If possible the granularity
    // should even be on page size (4096B == 512 doubles) to optimize the TLB misses.
    PermutationHelper &perm = PermutationHelper::instance();
    perm.ensureSize(M);
    int pi = 0;
    int pir = M - 1;
    for (int i = 0; i < M; ++i) {
        if (LINDXAU[i] >= 0) {
            perm[pi].a = i;
            perm[pi].b = LINDXAU[i];
            ++pi;
        } else {
            perm[pir].a = i;
            perm[pir].b = -LINDXAU[i];
            --pir;
        }
    }
#endif
    tbb::parallel_for (tbb::blocked_range<size_t>(0, N, 48),
            dlaswp01T_impl(M, A, LDA, U, LDU, LINDXA, LINDXAU),
            tbb::simple_partitioner());
    return;
#if 0
    const size_t prefetchStride = 90 * LDA;

    for (int i = 0; i < M; i++) {
        const int rowAr = LINDXA[i];
        double *Ar = &A[rowAr];
        if (LINDXAU[i] >= 0) {
            // good: copy one row from A into the transposed U matrix
            const size_t rowUw = LINDXAU[i];
            double *Uw = &U[rowUw * LDU];
            int NN = N;
#ifdef HPL_LASWP01T_STREAMING
            const double *null = 0;

            long offset = (Uw - null) & 7;
            if (N < offset) {
                offset = N;
            }
            int NN -= offset;
            offset = (7 - offset) * 0x12;
            double tmp;
            // ensure cache line alignment of Uw
            asm(
                "lea 0x6(%%rip), %[rip]\n"
                "lea (%[rip],%[tmp],1), %[rip]\n"
                "jmpq *%[rip]\n"
                "prefetchnta (%[Ar],%[prefetch],8)\n"
                "mov (%[Ar]), %[tmp]\n"
                "mov %[tmp], (%[Uw])\n"
                "add $0x8, %[Uw]\n"
                "lea (%[Ar], %[LDA], 8), %[Ar]\n"
                "prefetchnta (%[Ar],%[prefetch],8)\n"
                "mov (%[Ar]), %[tmp]\n"
                "mov %[tmp], (%[Uw])\n"
                "add $0x8, %[Uw]\n"
                "lea (%[Ar], %[LDA], 8), %[Ar]\n"
                "prefetchnta (%[Ar],%[prefetch],8)\n"
                "mov (%[Ar]), %[tmp]\n"
                "mov %[tmp], (%[Uw])\n"
                "add $0x8, %[Uw]\n"
                "lea (%[Ar], %[LDA], 8), %[Ar]\n"
                "prefetchnta (%[Ar],%[prefetch],8)\n"
                "mov (%[Ar]), %[tmp]\n"
                "mov %[tmp], (%[Uw])\n"
                "add $0x8, %[Uw]\n"
                "lea (%[Ar], %[LDA], 8), %[Ar]\n"
                "prefetchnta (%[Ar],%[prefetch],8)\n"
                "mov (%[Ar]), %[tmp]\n"
                "mov %[tmp], (%[Uw])\n"
                "add $0x8, %[Uw]\n"
                "lea (%[Ar], %[LDA], 8), %[Ar]\n"
                "prefetchnta (%[Ar],%[prefetch],8)\n"
                "mov (%[Ar]), %[tmp]\n"
                "mov %[tmp], (%[Uw])\n"
                "add $0x8, %[Uw]\n"
                "lea (%[Ar], %[LDA], 8), %[Ar]\n"
                "prefetchnta (%[Ar],%[prefetch],8)\n"
                "mov (%[Ar]), %[tmp]\n"
                "mov %[tmp], (%[Uw])\n"
                "add $0x8, %[Uw]\n"
                "lea (%[Ar], %[LDA], 8), %[Ar]\n"
                : [Ar]"+S"(Ar),
                  [Uw]"+D"(Uw),
                  [tmp]"+a"(offset), [rip]"=b"(tmp)
                : [prefetch]"c"(prefetchStride),
                  [LDA]"d"(static_cast<size_t>(LDA))
               );
#endif
            for (int col = NN / 16; col; --col) {
#ifdef HPL_LASWP01T_STREAMING
                double *Ar3 = Ar + 3 * LDA;
                double *prefetch0 = Ar + prefetchStride;
                double *prefetch3 = Ar3 + prefetchStride;
                asm(
                    "movsd (%[Ar0]), %%xmm0\n"
                    "prefetchnta (%[prefetch0])\n"
                    "movsd (%[Ar0],%[LDA8],2), %%xmm1\n"
                    "prefetchnta (%[prefetch0],%[LDA8],2)\n"
                    "movsd (%[Ar0],%[LDA8],4), %%xmm2\n"
                    "prefetchnta (%[prefetch0],%[LDA8],4)\n"
                    "movsd (%[Ar0],%[LDA8],8), %%xmm4\n"
                    "prefetchnta (%[prefetch0],%[LDA8],8)\n"
                    "movhpd (%[Ar0],%[LDA8],1), %%xmm0\n"
                    "prefetchnta (%[prefetch0],%[LDA8],1)\n"
                    "add %[LDA8x6], %[Ar0]\n"
                    "movhpd (%[Ar3]), %%xmm1\n"
                    "add %[LDA8x6], %[prefetch0]\n"
                    "prefetchnta (%[prefetch3])\n"
                    "movhpd (%[Ar3],%[LDA8],2), %%xmm2\n"
                    "prefetchnta (%[prefetch3],%[LDA8],2)\n"
                    "movsd (%[Ar0]), %%xmm3\n"
                    "prefetchnta (%[prefetch0])\n"
                    "movsd (%[Ar0],%[LDA8],4), %%xmm5\n"
                    "prefetchnta (%[prefetch0],%[LDA8],4)\n"
                    "movsd (%[Ar0],%[LDA8],8), %%xmm7\n"
                    "prefetchnta (%[prefetch0],%[LDA8],8)\n"
                    "add %[LDA8x6], %[Ar0]\n"
                    "movhpd (%[Ar3],%[LDA8],4), %%xmm3\n"
                    "add %[LDA8x6], %[prefetch0]\n"
                    "prefetchnta (%[prefetch3],%[LDA8],4)\n"
                    "movhpd (%[Ar3],%[LDA8],8), %%xmm5\n"
                    "prefetchnta (%[prefetch3],%[LDA8],8)\n"
                    "add %[LDA8x6], %[Ar3]\n"
                    "movhpd (%[Ar3]), %%xmm4\n"
                    "add %[LDA8x6], %[prefetch3]\n"
                    "prefetchnta (%[prefetch3])\n"
                    "add %[LDA8x6], %[Ar3]\n"
                    "movsd (%[Ar0]), %%xmm6\n"
                    "add %[LDA8x6], %[prefetch3]\n"
                    "prefetchnta (%[prefetch0])\n"
                    "movhpd (%[Ar0],%[LDA8],1), %%xmm6\n"
                    "prefetchnta (%[prefetch0],%[LDA8],1)\n"
                    "movhpd (%[Ar3]), %%xmm7\n"
                    "prefetchnta (%[prefetch3])\n"
                    "movntpd %%xmm0, (%[Uw])\n"
                    "movntpd %%xmm1, 0x10(%[Uw])\n"
                    "lea (%[Ar0],%[LDA8],4), %[Ar0]\n"
                    "movntpd %%xmm2, 0x20(%[Uw])\n"
                    "movntpd %%xmm3, 0x30(%[Uw])\n"
                    "lea (%[Ar3],%[LDA8],4), %[Ar3]\n"
                    "movntpd %%xmm4, 0x40(%[Uw])\n"
                    "movntpd %%xmm5, 0x50(%[Uw])\n"
                    "lea (%[prefetch0],%[LDA8],4), %[prefetch0]\n"
                    "movntpd %%xmm6, 0x60(%[Uw])\n"
                    "movntpd %%xmm7, 0x70(%[Uw])\n"
                    "lea (%[prefetch3],%[LDA8],4), %[prefetch3]\n"
                    : [Ar0]"+r"(Ar), [Ar3]"+r"(Ar3),
                      [prefetch0]"+r"(prefetch0),
                      [prefetch3]"+r"(prefetch3)
                    : [LDA8]"r"(static_cast<long>(LDA * 8)),
                      [LDA8x6]"r"(static_cast<long>(LDA * 8 * 6)),
                      [Uw]"r"(Uw)
                    : "xmm0", "xmm1", "xmm2", "xmm3",
                      "xmm4", "xmm5", "xmm6", "xmm7");
                Uw += 16;
#else
                Uw[ 0] = *Ar; Ar += LDA; Uw[ 1] = *Ar; Ar += LDA;
                Uw[ 2] = *Ar; Ar += LDA; Uw[ 3] = *Ar; Ar += LDA;
                Uw[ 4] = *Ar; Ar += LDA; Uw[ 5] = *Ar; Ar += LDA;
                Uw[ 6] = *Ar; Ar += LDA; Uw[ 7] = *Ar; Ar += LDA;
                Uw[ 8] = *Ar; Ar += LDA; Uw[ 9] = *Ar; Ar += LDA;
                Uw[10] = *Ar; Ar += LDA; Uw[11] = *Ar; Ar += LDA;
                Uw[12] = *Ar; Ar += LDA; Uw[13] = *Ar; Ar += LDA;
                Uw[14] = *Ar; Ar += LDA; Uw[15] = *Ar; Ar += LDA;
                Uw += 16;
#endif
            }
            for (int col = NN % 16; col; --col) {
                Uw[0] = Ar[0];
                ++Uw;
                Ar += LDA;
            }
        } else {
            // bad: copy one row from A to A
            const int rowAw = -LINDXAU[i];
            double *Aw = &A[rowAw];
            int col = N;
            do {
                //_mm_prefetch(Ar + prefetchStride1, _MM_HINT_NTA);
                //_m_prefetchw(Aw + prefetchStride1);
                *Aw = *Ar;
                Aw += LDA;
                Ar += LDA;
            } while (--col);
        }
    }
    //_mm_mfence();
#endif
}

#include "tsc.h"
#include <unistd.h>
#include <cstdlib>
#include <cstdio>

int main(int argc, char **argv)
{
    int nRows = 765;
    int nCols = 40449;
    int LDA = 41000;
    int LDU = 41000;
    int opt;
    while ((opt = getopt(argc, argv, "m:n:a:u:")) != -1) {
        switch (opt) {
        case 'm':
            nRows = atoi(optarg);
            break;
        case 'n':
            nCols = atoi(optarg);
            break;
        case 'a':
            LDA = atoi(optarg);
            break;
        case 'u':
            LDU = atoi(optarg);
            break;
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-m M] [-n N] [-a LDA] [-u LDU]\n",
                    argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    Matrix m1(LDA, LDA);
    Matrix m2 = m1;
    Matrix u1(LDU, LDU);
    Matrix u2 = u1;
    COMPARE(m1, m2);
    COMPARE(u1, u2);
    int *readRowsInA = new int[nRows];
    int *writeRowsInUA = new int[nRows];
    for (unsigned int i = 0; i < static_cast<unsigned int>(nRows); ++i) {
        readRowsInA[i] = (i * 363u + 102893u) % LDA;
        writeRowsInUA[i] = ((i * 7533059u + 8937563u) % (LDU + (LDA - 1))) - (LDA - 1);
    }
    TimeStampCounter tsc;
    tsc.start();
        dlaswp01T(nRows, nCols, m2.at(0, 0), m2.LDA(), u2.at(0, 0), u2.LDA(), readRowsInA, writeRowsInUA);
    tsc.stop();
    printf("new: %10lld cycles\n", tsc.cycles());
    tsc.start();
    HPL_dlaswp01T(nRows, nCols, m1.at(0, 0), m1.LDA(), u1.at(0, 0), u1.LDA(), readRowsInA, writeRowsInUA);
    tsc.stop();
    printf("old: %10lld cycles\n", tsc.cycles());
    COMPARE(m1, m2);
    COMPARE(u1, u2);
    return 0;
}
