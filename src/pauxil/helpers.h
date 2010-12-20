/*
 * Copyright 2010:
 *  - David Rohr (drohr@jwdt.org)
 *  - Matthias Bach (bach@compeng.uni-frankfurt.de)
 *  - Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 *
 * This file is part of HPL-GPU.
 *
 * HPL-GPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HPL-GPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HPL-GPU.  If not, see <http://www.gnu.org/licenses/>.
 *
 * In addition to the rules layed out by the GNU General Public License
 * the following exception is granted:
 *
 * Use with the Original BSD License.
 *
 * Notwithstanding any other provision of the GNU General Public License
 * Version 3, you have permission to link or combine any covered work with
 * a work licensed under the 4-clause BSD license into a single combined
 * work, and to convey the resulting work.  The terms of this License will
 * continue to apply to the part which is the covered work, but the special
 * requirements of the 4-clause BSD license, clause 3, concerning the
 * requirement of acknowledgement in advertising materials will apply to
 * the combination as such.
 */

#ifndef HELPERS_H
#define HELPERS_H

#ifdef HPL_HAVE_PREFETCHW
#include <mm3dnow.h>
#else
#define _m_prefetchw(addr) _mm_prefetch(addr, _MM_HINT_NTA)
#endif
#include <mmintrin.h>
#include <emmintrin.h>
#include <tbb/tbb_stddef.h>

namespace
{
    static inline void copy(double *__restrict__ dst, const double *__restrict__ src)
    {
        // use general purpose registers:
        //*dst = *src;

        // use MMX registers:
        //__m64 tmp;
        //asm("movq (%[src]), %[tmp]" : [tmp]"=y"(tmp) : [src]"r"(src));
        //asm("movq %[tmp], (%[dst])" :: [tmp]"y"(tmp), [dst]"r"(dst));

        // use SSE registers:
        __m128 tmp;
        asm("movq %[src], %[tmp]" : [tmp]"=x"(tmp) : [src]"m"(*src));
        asm("movq %[tmp], %[dst]" :: [tmp]"x"(tmp), [dst]"m"(*dst));
    }

    static inline void streamingCopy(double *__restrict__ dst, const double *__restrict__ src)
    {
        // use general purpose registers + streaming stores:
        //asm("movnti %0, (%1)"::"r"(*src), "r"(dst));

        // use MMX registers + streaming stores:
        const __m64 tmp = *reinterpret_cast<const __m64 *>(src);
        //asm("movq %[src], %[tmp]" : [tmp]"=y"(tmp) : [src]"m"(*src));
        //asm("movntq %[tmp], %[dst]" :: [tmp]"y"(tmp), [dst]"m"(*dst));
        asm("movntq %[src], %[dst]" :: [src]"y"(tmp), [dst]"m"(*dst));
        //_mm_stream_pi(reinterpret_cast<__m64 *>(dst), *reinterpret_cast<const __m64 *>(src));
    }

    enum {
        CachelineSize = 64, // Bytes
        DoublesInCacheline = CachelineSize / sizeof(double)
    };

    template<typename T> static inline T max(T a, T b) { return a > b ? a : b; }

    template<typename T> static inline void swap(T &__restrict__ a, T &__restrict__ b)
    {
        register T tmp = a;
        a = b;
        b = tmp;
    }

    static inline void swapSSE(double &__restrict__ a, double &__restrict__ b)
    {
        const __m128d va = _mm_load_pd(&a);
        const __m128d vb = _mm_load_pd(&b);
        _mm_store_pd(&a, vb);
        _mm_store_pd(&b, va);
    }

    template<size_t MultipleOf, size_t Blocksize>
    class MyRange
    {
        public:
            MyRange(size_t b, size_t n, size_t blocksize = Blocksize)
                : m_begin(b), m_n(n), m_blocksize(blocksize)
            {}

            MyRange(MyRange<MultipleOf, Blocksize> &r, tbb::split)
                : m_begin(r.m_begin),
                m_n((r.m_n / 2) & ~(MultipleOf - 1)),
                m_blocksize(r.m_blocksize)
            {
                r.m_begin += m_n;
                r.m_n -= m_n;
            }

            bool empty() const { return m_n == 0; }
            bool is_divisible() const { return m_n >= m_blocksize; }

            size_t N() const { return m_n; }
            size_t begin() const { return m_begin; }

        private:
            size_t m_begin;
            size_t m_n;
            size_t m_blocksize;
    };

} // anonymous namespace

#endif // HELPERS_H
