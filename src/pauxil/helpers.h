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

#ifndef HELPERS_H
#define HELPERS_H

#include <mm3dnow.h>
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
        asm("movq (%[src]), %[tmp]" : [tmp]"=x"(tmp) : [src]"r"(src));
        asm("movq %[tmp], (%[dst])" :: [tmp]"x"(tmp), [dst]"r"(dst));
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
