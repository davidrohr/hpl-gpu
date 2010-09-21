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

#include <cstddef>
#include "util_timer.h"
#include "util_trace.h"
#include <tbb/parallel_for.h>
#include <emmintrin.h>
#include <mm3dnow.h>

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
    /* Trace:
                  M      N    LDA    LDU
HPL_dlaswp01T(  765, 40449, 40960, 40449) 3613527484 cycles
HPL_dlaswp01T(  765, 40960, 40960, 40960) 3788513393 cycles
HPL_dlaswp01T(  758, 40449, 40960, 40449) 3483699888 cycles
HPL_dlaswp01T(  758, 40448, 40960, 40448) 3723775397 cycles
HPL_dlaswp01T(  760, 40448, 40960, 40448) 3670094564 cycles
HPL_dlaswp01T(  760, 39937, 40960, 39937) 3736038645 cycles
HPL_dlaswp01T(  768, 39936, 40960, 39936) 3741281802 cycles
HPL_dlaswp01T(  768, 39937, 40960, 39937) 3783777078 cycles
HPL_dlaswp01T(  769, 39936, 40960, 39936) 3515353406 cycles
HPL_dlaswp01T(  769, 39425, 40960, 39425) 3734103373 cycles
HPL_dlaswp01T(  764, 39424, 40960, 39424) 3395104659 cycles
HPL_dlaswp01T(  764, 39425, 40960, 39425) 3708630845 cycles
HPL_dlaswp01T(  748, 38913, 40960, 38913) 3491349392 cycles
HPL_dlaswp01T(  748, 39424, 40960, 39424) 3495775492 cycles
HPL_dlaswp01T(  765, 38913, 40960, 38913) 3636834147 cycles
HPL_dlaswp01T(  765, 38912, 40960, 38912) 3674558332 cycles
HPL_dlaswp01T(  767, 38401, 40960, 38401) 3435449941 cycles
HPL_dlaswp01T(  767, 38912, 40960, 38912) 3660481737 cycles
HPL_dlaswp01T(  784, 38401, 40960, 38401) 3520025543 cycles
HPL_dlaswp01T(  784, 38400, 40960, 38400) 3746445120 cycles
HPL_dlaswp01T(  774, 37889, 40960, 37889) 3668879890 cycles
HPL_dlaswp01T(  774, 38400, 40960, 38400) 3690425691 cycles
HPL_dlaswp01T(  760, 37889, 40960, 37889) 3496570253 cycles
HPL_dlaswp01T(  760, 37888, 40960, 37888) 3508029130 cycles
HPL_dlaswp01T(  766, 37888, 40960, 37888) 3272301725 cycles
HPL_dlaswp01T(  766, 37377, 40960, 37377) 3561612846 cycles
HPL_dlaswp01T(  794, 37376, 40960, 37376) 3494222479 cycles
HPL_dlaswp01T(  794, 37377, 40960, 37377) 3760732852 cycles
HPL_dlaswp01T(  765, 36865, 40960, 36865) 3428482763 cycles
HPL_dlaswp01T(  765, 37376, 40960, 37376) 3452253006 cycles
HPL_dlaswp01T(  754, 36864, 40960, 36864) 3282475141 cycles
HPL_dlaswp01T(  754, 36865, 40960, 36865) 3385212882 cycles
HPL_dlaswp01T(  753, 36353, 40960, 36353) 3076836988 cycles
HPL_dlaswp01T(  753, 36864, 40960, 36864) 3397396362 cycles
HPL_dlaswp01T(  784, 36353, 40960, 36353) 3338748291 cycles
HPL_dlaswp01T(  784, 36352, 40960, 36352) 3561338115 cycles
HPL_dlaswp01T(  753, 35841, 40960, 35841) 3217706129 cycles
HPL_dlaswp01T(  753, 36352, 40960, 36352) 3272941242 cycles
HPL_dlaswp01T(  780, 35840, 40960, 35840) 3432933482 cycles
HPL_dlaswp01T(  780, 35841, 40960, 35841) 3492859841 cycles
HPL_dlaswp01T(  769, 35840, 40960, 35840) 3163217819 cycles
HPL_dlaswp01T(  769, 35329, 40960, 35329) 3333799489 cycles
HPL_dlaswp01T(  759, 35328, 40960, 35328) 3015082975 cycles
HPL_dlaswp01T(  759, 35329, 40960, 35329) 3237167838 cycles
HPL_dlaswp01T(  772, 35328, 40960, 35328) 3341786802 cycles
HPL_dlaswp01T(  772, 34817, 40960, 34817) 3378614438 cycles
HPL_dlaswp01T(  747, 34816, 40960, 34816) 3084773950 cycles
HPL_dlaswp01T(  747, 34817, 40960, 34817) 3130035082 cycles
HPL_dlaswp01T(  763, 34305, 40960, 34305) 3030669442 cycles
HPL_dlaswp01T(  763, 34816, 40960, 34816) 3205026486 cycles
HPL_dlaswp01T(  747, 34305, 40960, 34305) 2879958805 cycles
HPL_dlaswp01T(  747, 34304, 40960, 34304) 3036998021 cycles
HPL_dlaswp01T(  763, 33793, 40960, 33793) 3139629629 cycles
HPL_dlaswp01T(  763, 34304, 40960, 34304) 3164347712 cycles
HPL_dlaswp01T(  761, 33793, 40960, 33793) 3129222181 cycles
HPL_dlaswp01T(  761, 33792, 40960, 33792) 3142567679 cycles
HPL_dlaswp01T(  759, 33792, 40960, 33792) 2882523894 cycles
HPL_dlaswp01T(  759, 33281, 40960, 33281) 3086839486 cycles
HPL_dlaswp01T(  767, 33280, 40960, 33280) 2884463999 cycles
HPL_dlaswp01T(  767, 33281, 40960, 33281) 3174012214 cycles
HPL_dlaswp01T(  762, 33280, 40960, 33280) 3073419741 cycles
HPL_dlaswp01T(  762, 32769, 40960, 32769) 3079061465 cycles
HPL_dlaswp01T(  780, 32768, 40960, 32768) 3116570366 cycles
HPL_dlaswp01T(  780, 32769, 40960, 32769) 3173698594 cycles
HPL_dlaswp01T(  781, 32257, 40960, 32257) 2905498967 cycles
HPL_dlaswp01T(  781, 32768, 40960, 32768) 3191945317 cycles
HPL_dlaswp01T(  747, 32257, 40960, 32257) 2678274452 cycles
HPL_dlaswp01T(  747, 32256, 40960, 32256) 2851538132 cycles
HPL_dlaswp01T(  767, 31745, 40960, 31745) 2975797205 cycles
HPL_dlaswp01T(  767, 32256, 40960, 32256) 3059376134 cycles
HPL_dlaswp01T(  765, 31745, 40960, 31745) 2886882448 cycles
HPL_dlaswp01T(  765, 31744, 40960, 31744) 2909639837 cycles
HPL_dlaswp01T(  737, 31744, 40960, 31744) 2493187031 cycles
HPL_dlaswp01T(  737, 31233, 40960, 31233) 2761998009 cycles
HPL_dlaswp01T(  797, 31232, 40960, 31232) 2940255312 cycles
HPL_dlaswp01T(  797, 31233, 40960, 31233) 3167732533 cycles
HPL_dlaswp01T(  771, 30721, 40960, 30721) 2892225503 cycles
HPL_dlaswp01T(  771, 31232, 40960, 31232) 2901928637 cycles
HPL_dlaswp01T(  775, 30721, 40960, 30721) 2880926592 cycles
HPL_dlaswp01T(  775, 30720, 40960, 30720) 2907228517 cycles
HPL_dlaswp01T(  737, 30209, 40960, 30209) 2435517227 cycles
HPL_dlaswp01T(  737, 30720, 40960, 30720) 2638928752 cycles
HPL_dlaswp01T(  764, 30209, 40960, 30209) 2617745526 cycles
HPL_dlaswp01T(  764, 30208, 40960, 30208) 2782966349 cycles
HPL_dlaswp01T(  769, 29697, 40960, 29697) 2811612965 cycles
HPL_dlaswp01T(  769, 30208, 40960, 30208) 2830453262 cycles
HPL_dlaswp01T(  771, 29696, 40960, 29696) 2753515334 cycles
HPL_dlaswp01T(  771, 29697, 40960, 29697) 2789223816 cycles
HPL_dlaswp01T(  755, 29696, 40960, 29696) 2488264506 cycles
HPL_dlaswp01T(  755, 29185, 40960, 29185) 2678550722 cycles
HPL_dlaswp01T(  777, 29184, 40960, 29184) 2610575017 cycles
HPL_dlaswp01T(  777, 29185, 40960, 29185) 2819114795 cycles
HPL_dlaswp01T(  759, 29184, 40960, 29184) 2626931923 cycles
HPL_dlaswp01T(  759, 28673, 40960, 28673) 2653408629 cycles
HPL_dlaswp01T(  775, 28672, 40960, 28672) 2718252883 cycles
HPL_dlaswp01T(  775, 28673, 40960, 28673) 2732140618 cycles
HPL_dlaswp01T(  774, 28161, 40960, 28161) 2535678714 cycles
HPL_dlaswp01T(  774, 28672, 40960, 28672) 2720695981 cycles
HPL_dlaswp01T(  752, 28161, 40960, 28161) 2377171220 cycles
HPL_dlaswp01T(  752, 28160, 40960, 28160) 2510246896 cycles
HPL_dlaswp01T(  782, 27649, 40960, 27649) 2693178806 cycles
HPL_dlaswp01T(  782, 28160, 40960, 28160) 2704613312 cycles
HPL_dlaswp01T(  778, 27649, 40960, 27649) 2661270091 cycles
HPL_dlaswp01T(  778, 27648, 40960, 27648) 2661742594 cycles
HPL_dlaswp01T(  763, 27648, 40960, 27648) 2351451733 cycles
HPL_dlaswp01T(  763, 27137, 40960, 27137) 2566960123 cycles
HPL_dlaswp01T(  742, 27136, 40960, 27136) 2188262368 cycles
HPL_dlaswp01T(  742, 27137, 40960, 27137) 2382753140 cycles
HPL_dlaswp01T(  753, 26625, 40960, 26625) 2384274741 cycles
HPL_dlaswp01T(  753, 27136, 40960, 27136) 2421644086 cycles
HPL_dlaswp01T(  775, 26624, 40960, 26624) 2486619900 cycles
HPL_dlaswp01T(  775, 26625, 40960, 26625) 2528492734 cycles
HPL_dlaswp01T(  760, 26113, 40960, 26113) 2244751492 cycles
HPL_dlaswp01T(  760, 26624, 40960, 26624) 2412790220 cycles
HPL_dlaswp01T(  757, 26113, 40960, 26113) 2253519270 cycles
HPL_dlaswp01T(  757, 26112, 40960, 26112) 2353044930 cycles
HPL_dlaswp01T(  753, 26112, 40960, 26112) 2312370254 cycles
HPL_dlaswp01T(  753, 25601, 40960, 25601) 2323958862 cycles
HPL_dlaswp01T(  746, 25600, 40960, 25600) 2215853814 cycles
HPL_dlaswp01T(  746, 25601, 40960, 25601) 2283464807 cycles
HPL_dlaswp01T(  762, 25600, 40960, 25600) 2175167878 cycles
HPL_dlaswp01T(  762, 25089, 40960, 25089) 2332773671 cycles
HPL_dlaswp01T(  777, 25088, 40960, 25088) 2200506406 cycles
HPL_dlaswp01T(  777, 25089, 40960, 25089) 2380087092 cycles
HPL_dlaswp01T(  776, 24577, 40960, 24577) 2284126482 cycles
HPL_dlaswp01T(  776, 25088, 40960, 25088) 2371294860 cycles
HPL_dlaswp01T(  783, 24576, 40960, 24576) 2382433639 cycles
HPL_dlaswp01T(  783, 24577, 40960, 24577) 2392611972 cycles
HPL_dlaswp01T(  761, 24065, 40960, 24065) 2091258571 cycles
HPL_dlaswp01T(  761, 24576, 40960, 24576) 2271558238 cycles
HPL_dlaswp01T(  770, 24065, 40960, 24065) 2111046931 cycles
HPL_dlaswp01T(  770, 24064, 40960, 24064) 2249944641 cycles
HPL_dlaswp01T(  763, 23553, 40960, 23553) 2148248749 cycles
HPL_dlaswp01T(  763, 24064, 40960, 24064) 2193514187 cycles
HPL_dlaswp01T(  777, 23552, 40960, 23552) 2221226730 cycles
HPL_dlaswp01T(  777, 23553, 40960, 23553) 2228368042 cycles
HPL_dlaswp01T(  773, 23552, 40960, 23552) 2072811519 cycles
HPL_dlaswp01T(  773, 23041, 40960, 23041) 2188325824 cycles
HPL_dlaswp01T(  774, 23040, 40960, 23040) 1991401983 cycles
HPL_dlaswp01T(  774, 23041, 40960, 23041) 2170372470 cycles
HPL_dlaswp01T(  760, 23040, 40960, 23040) 2059405561 cycles
HPL_dlaswp01T(  760, 22529, 40960, 22529) 2078161790 cycles
HPL_dlaswp01T(  759, 22528, 40960, 22528) 2017269839 cycles
HPL_dlaswp01T(  759, 22529, 40960, 22529) 2055162550 cycles
HPL_dlaswp01T(  760, 22017, 40960, 22017) 1906763757 cycles
HPL_dlaswp01T(  760, 22528, 40960, 22528) 2057458829 cycles
HPL_dlaswp01T(  784, 22017, 40960, 22017) 2029461085 cycles
HPL_dlaswp01T(  784, 22016, 40960, 22016) 2133889811 cycles
HPL_dlaswp01T(  753, 22016, 40960, 22016) 1961580635 cycles
HPL_dlaswp01T(  753, 21505, 40960, 21505) 1964430986 cycles
HPL_dlaswp01T(  776, 21504, 40960, 21504) 2025183601 cycles
HPL_dlaswp01T(  776, 21505, 40960, 21505) 2034550948 cycles
HPL_dlaswp01T(  773, 21504, 40960, 21504) 1855063492 cycles
HPL_dlaswp01T(  773, 20993, 40960, 20993) 1978454698 cycles
HPL_dlaswp01T(  772, 20992, 40960, 20992) 1821649668 cycles
HPL_dlaswp01T(  772, 20993, 40960, 20993) 1960669226 cycles
HPL_dlaswp01T(  759, 20481, 40960, 20481) 1872966416 cycles
HPL_dlaswp01T(  759, 20992, 40960, 20992) 1882306151 cycles
HPL_dlaswp01T(  783, 20480, 40960, 20480) 1916396907 cycles
HPL_dlaswp01T(  783, 20481, 40960, 20481) 1964314831 cycles
HPL_dlaswp01T(  776, 19969, 40960, 19969) 1824699230 cycles
HPL_dlaswp01T(  776, 20480, 40960, 20480) 1941673496 cycles
HPL_dlaswp01T(  770, 19969, 40960, 19969) 1761423255 cycles
HPL_dlaswp01T(  770, 19968, 40960, 19968) 1826808519 cycles
HPL_dlaswp01T(  747, 19457, 40960, 19457) 1720117767 cycles
HPL_dlaswp01T(  747, 19968, 40960, 19968) 1739413150 cycles
HPL_dlaswp01T(  778, 19456, 40960, 19456) 1825174067 cycles
HPL_dlaswp01T(  778, 19457, 40960, 19457) 1862782090 cycles
HPL_dlaswp01T(  756, 19456, 40960, 19456) 1618824737 cycles
HPL_dlaswp01T(  756, 18945, 40960, 18945) 1711113287 cycles
HPL_dlaswp01T(  779, 18944, 40960, 18944) 1668720103 cycles
HPL_dlaswp01T(  779, 18945, 40960, 18945) 1798084346 cycles
HPL_dlaswp01T(  754, 18433, 40960, 18433) 1656142729 cycles
HPL_dlaswp01T(  754, 18944, 40960, 18944) 1688500064 cycles
HPL_dlaswp01T(  747, 18432, 40960, 18432) 1586131505 cycles
HPL_dlaswp01T(  747, 18433, 40960, 18433) 1603785700 cycles
HPL_dlaswp01T(  763, 17921, 40960, 17921) 1551953414 cycles
HPL_dlaswp01T(  763, 18432, 40960, 18432) 1656174540 cycles
HPL_dlaswp01T(  756, 17921, 40960, 17921) 1495307354 cycles
HPL_dlaswp01T(  756, 17920, 40960, 17920) 1578059147 cycles
HPL_dlaswp01T(  759, 17409, 40960, 17409) 1567160514 cycles
HPL_dlaswp01T(  759, 17920, 40960, 17920) 1592174696 cycles
HPL_dlaswp01T(  763, 17409, 40960, 17409) 1565360591 cycles
HPL_dlaswp01T(  763, 17408, 40960, 17408) 1584094046 cycles
HPL_dlaswp01T(  763, 17408, 40960, 17408) 1455445181 cycles
HPL_dlaswp01T(  763, 16897, 40960, 16897) 1575790444 cycles
HPL_dlaswp01T(  759, 16896, 40960, 16896) 1388923274 cycles
HPL_dlaswp01T(  759, 16897, 40960, 16897) 1532281402 cycles
HPL_dlaswp01T(  766, 16385, 40960, 16385) 1520793690 cycles
HPL_dlaswp01T(  766, 16896, 40960, 16896) 1526286528 cycles
HPL_dlaswp01T(  757, 16384, 40960, 16384) 1427995428 cycles
HPL_dlaswp01T(  757, 16385, 40960, 16385) 1471948120 cycles
HPL_dlaswp01T(  764, 15873, 40960, 15873) 1374614019 cycles
HPL_dlaswp01T(  764, 16384, 40960, 16384) 1461125574 cycles
HPL_dlaswp01T(  777, 15873, 40960, 15873) 1405605652 cycles
HPL_dlaswp01T(  777, 15872, 40960, 15872) 1480488922 cycles
HPL_dlaswp01T(  753, 15361, 40960, 15361) 1356465750 cycles
HPL_dlaswp01T(  753, 15872, 40960, 15872) 1402102161 cycles
HPL_dlaswp01T(  754, 15360, 40960, 15360) 1334708219 cycles
HPL_dlaswp01T(  754, 15361, 40960, 15361) 1372287491 cycles
HPL_dlaswp01T(  765, 15360, 40960, 15360) 1288583871 cycles
HPL_dlaswp01T(  765, 14849, 40960, 14849) 1374419989 cycles
HPL_dlaswp01T(  757, 14848, 40960, 14848) 1200954038 cycles
HPL_dlaswp01T(  757, 14849, 40960, 14849) 1355153053 cycles
HPL_dlaswp01T(  763, 14337, 40960, 14337) 1308527963 cycles
HPL_dlaswp01T(  763, 14848, 40960, 14848) 1337180239 cycles
HPL_dlaswp01T(  764, 14336, 40960, 14336) 1280731968 cycles
HPL_dlaswp01T(  764, 14337, 40960, 14337) 1314193904 cycles
HPL_dlaswp01T(  759, 13825, 40960, 13825) 1165271054 cycles
HPL_dlaswp01T(  759, 14336, 40960, 14336) 1269465309 cycles
HPL_dlaswp01T(  742, 13825, 40960, 13825) 1119928021 cycles
HPL_dlaswp01T(  742, 13824, 40960, 13824) 1181011361 cycles
HPL_dlaswp01T(  773, 13824, 40960, 13824) 1255736195 cycles
HPL_dlaswp01T(  773, 13313, 40960, 13313) 1257371440 cycles
HPL_dlaswp01T(  784, 13312, 40960, 13312) 1227786396 cycles
HPL_dlaswp01T(  784, 13313, 40960, 13313) 1273256432 cycles
HPL_dlaswp01T(  747, 13312, 40960, 13312) 1052823038 cycles
HPL_dlaswp01T(  747, 12801, 40960, 12801) 1136075920 cycles
HPL_dlaswp01T(  773, 12800, 40960, 12800) 1082266713 cycles
HPL_dlaswp01T(  773, 12801, 40960, 12801) 1197326998 cycles
HPL_dlaswp01T(  768, 12289, 40960, 12289) 1132137215 cycles
HPL_dlaswp01T(  768, 12800, 40960, 12800) 1151053667 cycles
HPL_dlaswp01T(  770, 12288, 40960, 12288) 1080514933 cycles
HPL_dlaswp01T(  770, 12289, 40960, 12289) 1140601180 cycles
HPL_dlaswp01T(  750, 11777, 40960, 11777) 986898727 cycles
HPL_dlaswp01T(  750, 12288, 40960, 12288) 1050705701 cycles
HPL_dlaswp01T(  784, 11777, 40960, 11777) 1053704262 cycles
HPL_dlaswp01T(  784, 11776, 40960, 11776) 1099839038 cycles
HPL_dlaswp01T(  761, 11265, 40960, 11265) 1020035573 cycles
HPL_dlaswp01T(  761, 11776, 40960, 11776) 1028461527 cycles
HPL_dlaswp01T(  771, 11264, 40960, 11264) 998336766 cycles
HPL_dlaswp01T(  771, 11265, 40960, 11265) 1036648529 cycles
HPL_dlaswp01T(  738, 11264, 40960, 11264) 856437103 cycles
HPL_dlaswp01T(  738, 10753, 40960, 10753) 920820884 cycles
HPL_dlaswp01T(  765, 10752, 40960, 10752) 874665921 cycles
HPL_dlaswp01T(  765, 10753, 40960, 10753) 980734929 cycles
HPL_dlaswp01T(  761, 10241, 40960, 10241) 912708262 cycles
HPL_dlaswp01T(  761, 10752, 40960, 10752) 920313466 cycles
HPL_dlaswp01T(  772, 10240, 40960, 10240) 916154057 cycles
HPL_dlaswp01T(  772, 10241, 40960, 10241) 939702498 cycles
HPL_dlaswp01T(  754,  9729, 40960,  9729) 809645121 cycles
HPL_dlaswp01T(  754, 10240, 40960, 10240) 884895339 cycles
HPL_dlaswp01T(  764,  9729, 40960,  9729) 824901049 cycles
HPL_dlaswp01T(  764,  9728, 40960,  9728) 856925635 cycles
HPL_dlaswp01T(  762,  9217, 40960,  9217) 831319789 cycles
HPL_dlaswp01T(  762,  9728, 40960,  9728) 854869110 cycles
HPL_dlaswp01T(  744,  9216, 40960,  9216) 758627096 cycles
HPL_dlaswp01T(  744,  9217, 40960,  9217) 777725900 cycles
HPL_dlaswp01T(  752,  9216, 40960,  9216) 714487301 cycles
HPL_dlaswp01T(  752,  8705, 40960,  8705) 759587616 cycles
HPL_dlaswp01T(  782,  8704, 40960,  8704) 721485696 cycles
HPL_dlaswp01T(  782,  8705, 40960,  8705) 809365769 cycles
HPL_dlaswp01T(  754,  8193, 40960,  8193) 719672611 cycles
HPL_dlaswp01T(  754,  8704, 40960,  8704) 721866679 cycles
HPL_dlaswp01T(  763,  8192, 40960,  8192) 690240859 cycles
HPL_dlaswp01T(  763,  8193, 40960,  8193) 721778109 cycles
HPL_dlaswp01T(  734,  7681, 40960,  7681) 608998007 cycles
HPL_dlaswp01T(  734,  8192, 40960,  8192) 659988671 cycles
HPL_dlaswp01T(  765,  7681, 40960,  7681) 642937765 cycles
HPL_dlaswp01T(  765,  7680, 40960,  7680) 656084207 cycles
HPL_dlaswp01T(  759,  7169, 40960,  7169) 622520102 cycles
HPL_dlaswp01T(  759,  7680, 40960,  7680) 629487864 cycles
HPL_dlaswp01T(  761,  7168, 40960,  7168) 591591779 cycles
HPL_dlaswp01T(  761,  7169, 40960,  7169) 626303075 cycles
HPL_dlaswp01T(  759,  7168, 40960,  7168) 557135907 cycles
HPL_dlaswp01T(  759,  6657, 40960,  6657) 580647073 cycles
HPL_dlaswp01T(  762,  6656, 40960,  6656) 516213652 cycles
HPL_dlaswp01T(  762,  6657, 40960,  6657) 578966591 cycles
HPL_dlaswp01T(  735,  6145, 40960,  6145) 470687553 cycles
HPL_dlaswp01T(  735,  6656, 40960,  6656) 521652782 cycles
HPL_dlaswp01T(  753,  6144, 40960,  6144) 463844377 cycles
HPL_dlaswp01T(  753,  6145, 40960,  6145) 515954815 cycles
HPL_dlaswp01T(  735,  5633, 40960,  5633) 405944140 cycles
HPL_dlaswp01T(  735,  6144, 40960,  6144) 476865984 cycles
HPL_dlaswp01T(  748,  5632, 40960,  5632) 418754769 cycles
HPL_dlaswp01T(  748,  5633, 40960,  5633) 458028859 cycles
HPL_dlaswp01T(  739,  5121, 40960,  5121) 388023757 cycles
HPL_dlaswp01T(  739,  5632, 40960,  5632) 438462696 cycles
HPL_dlaswp01T(  743,  5120, 40960,  5120) 379563882 cycles
HPL_dlaswp01T(  743,  5121, 40960,  5121) 422532880 cycles
HPL_dlaswp01T(  741,  4609, 40960,  4609) 351675524 cycles
HPL_dlaswp01T(  741,  5120, 40960,  5120) 374415591 cycles
HPL_dlaswp01T(  768,  4608, 40960,  4608) 331486675 cycles
HPL_dlaswp01T(  768,  4609, 40960,  4609) 404277485 cycles
HPL_dlaswp01T(  746,  4097, 40960,  4097) 316113717 cycles
HPL_dlaswp01T(  746,  4608, 40960,  4608) 363761049 cycles
HPL_dlaswp01T(  747,  4096, 40960,  4096) 296498757 cycles
HPL_dlaswp01T(  747,  4097, 40960,  4097) 331734493 cycles
HPL_dlaswp01T(  741,  3585, 40960,  3585) 266617398 cycles
HPL_dlaswp01T(  741,  4096, 40960,  4096) 312450598 cycles
HPL_dlaswp01T(  740,  3584, 40960,  3584) 249954355 cycles
HPL_dlaswp01T(  740,  3585, 40960,  3585) 261640477 cycles
HPL_dlaswp01T(  738,  3073, 40960,  3073) 227703587 cycles
HPL_dlaswp01T(  738,  3584, 40960,  3584) 244570257 cycles
HPL_dlaswp01T(  757,  3072, 40960,  3072) 218105271 cycles
HPL_dlaswp01T(  757,  3073, 40960,  3073) 233194396 cycles
HPL_dlaswp01T(  733,  2561, 40960,  2561) 181272547 cycles
HPL_dlaswp01T(  733,  3072, 40960,  3072) 185257044 cycles
HPL_dlaswp01T(  753,  2560, 40960,  2560) 163056045 cycles
HPL_dlaswp01T(  753,  2561, 40960,  2561) 185039304 cycles
HPL_dlaswp01T(  724,  2049, 40960,  2049) 134802567 cycles
HPL_dlaswp01T(  724,  2560, 40960,  2560) 159511115 cycles
HPL_dlaswp01T(  732,  2048, 40960,  2048) 131207195 cycles
HPL_dlaswp01T(  732,  2049, 40960,  2049) 138586960 cycles
HPL_dlaswp01T(  709,  1537, 40960,  1537) 96035884 cycles
HPL_dlaswp01T(  709,  2048, 40960,  2048) 119532712 cycles
HPL_dlaswp01T(  731,  1536, 40960,  1536) 91356051 cycles
HPL_dlaswp01T(  731,  1537, 40960,  1537) 95322886 cycles
HPL_dlaswp01T(  676,  1025, 40960,  1025) 57317937 cycles
HPL_dlaswp01T(  676,  1536, 40960,  1536) 75066301 cycles
HPL_dlaswp01T(  709,  1024, 40960,  1024) 53200438 cycles
HPL_dlaswp01T(  709,  1025, 40960,  1025) 57350342 cycles
HPL_dlaswp01T(  648,   513, 40960,   513) 23203072 cycles
HPL_dlaswp01T(  648,  1024, 40960,  1024) 38457833 cycles
HPL_dlaswp01T(  692,   512, 40960,   512) 21671135 cycles
HPL_dlaswp01T(  692,   513, 40960,   513) 24432011 cycles
HPL_dlaswp01T(  512,     1, 40960,     1) 14608 cycles
HPL_dlaswp01T(  512,   512, 40960,   512) 12082933 cycles
HPL_dlaswp01T(  512,     1, 40960,     1) 10444 cycles
     */
#ifdef TRACE_CALLS
   uint64_t tr_start, tr_end, tr_diff;
   tr_start = util_getTimestamp();
#endif /* TRACE_CALLS */

#ifdef USE_ORIGINAL_LASWP
#include "HPL_dlaswp01T.c"
#else
    if( N < 16 ) {
        if (N <= 0) {
            return;
        }

        for (int i = 0; i < M; ++i) {
            const double *Ar = &A[LINDXA[i]];
            if (LINDXAU[i] >= 0) {
                double *Uw = &U[LINDXAU[i] * static_cast<size_t>(LDU)];
                for (int col = 0; col < N; ++col) {
                    Uw[col] = *Ar;
                    Ar += LDA;
                }
            } else {
                double *Aw = &A[-LINDXAU[i]];
                for (int col = 0; col < N; ++col) {
                    *Aw = *Ar;
                    Ar += LDA;
                    Aw += LDA;
                }
            }
        }
    } else {
        tbb::parallel_for (tbb::blocked_range<size_t>(0, N, 48),
                dlaswp01T_impl(M, A, LDA, U, LDU, LINDXA, LINDXAU),
                tbb::simple_partitioner());
    }
#endif

#ifdef TRACE_CALLS
   tr_end = util_getTimestamp();
   tr_diff = util_getTimeDifference( tr_start, tr_end );

   fprintf( trace_dgemm, "DLASWP01T,M=%i,N=%i,LDA=%i,TIME=%lu,THRPT=%.2fGB/s\n", M, N, LDA, tr_diff,
           0.002 * sizeof(double) * M * N / tr_diff );
#endif /* TRACE_CALLS */
} 
