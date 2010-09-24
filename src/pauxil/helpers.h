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

namespace
{
    void streamingCopy(double *__restrict__ dst, const double *__restrict__ src)
    {
        // use general purpose registers:
        //*dst = *src;

        // use general purpose registers + streaming stores:
        //asm("movnti %0, (%1)"::"r"(*src), "r"(dst));

        // use MMX registers:
        //__m64 tmp;
        //asm("movq (%[src]), %[tmp]" : [tmp]"=y"(tmp) : [src]"r"(src));
        //asm("movq %[tmp], (%[dst])" :: [tmp]"y"(tmp), [dst]"r"(dst));

        // use MMX registers + streaming stores:
        //_mm_stream_pi(reinterpret_cast<__m64 *>(dst), _mm_cvtsi64_m64(*reinterpret_cast<const long long *>(src)));

        // use SSE registers:
        __m128 tmp;
        asm("movq (%[src]), %[tmp]" : [tmp]"=x"(tmp) : [src]"r"(src));
        asm("movq %[tmp], (%[dst])" :: [tmp]"x"(tmp), [dst]"r"(dst));
    }

    enum {
        CachelineSize = 64, // Bytes
        DoublesInCacheline = CachelineSize / sizeof(double)
    };

} // anonymous namespace

#endif // HELPERS_H
