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
#include "helpers.h"

/*
 * Trace:
HPL_dlaswp06T(  258, 40449, 40960, 40449) 2796577975 cycles
HPL_dlaswp06T(  258, 40960, 40960, 40960) 2859118178 cycles
HPL_dlaswp06T(  266, 40449, 40960, 40449) 2910295893 cycles
HPL_dlaswp06T(  266, 40448, 40960, 40448) 2863313752 cycles
HPL_dlaswp06T(  261, 40448, 40960, 40448) 2849317322 cycles
HPL_dlaswp06T(  261, 39937, 40960, 39937) 2882081421 cycles
HPL_dlaswp06T(  256, 39936, 40960, 39936) 2753931278 cycles
HPL_dlaswp06T(  256, 39937, 40960, 39937) 2785773433 cycles
HPL_dlaswp06T(  251, 39936, 40960, 39936) 2709539751 cycles
HPL_dlaswp06T(  251, 39425, 40960, 39425) 2689442789 cycles
HPL_dlaswp06T(  258, 39424, 40960, 39424) 2744551886 cycles
HPL_dlaswp06T(  258, 39425, 40960, 39425) 2771001586 cycles
HPL_dlaswp06T(  275, 38913, 40960, 38913) 2905956311 cycles
HPL_dlaswp06T(  275, 39424, 40960, 39424) 2909001902 cycles
HPL_dlaswp06T(  256, 38912, 40960, 38912) 2630473412 cycles
HPL_dlaswp06T(  256, 38913, 40960, 38913) 2719029966 cycles
HPL_dlaswp06T(  257, 38401, 40960, 38401) 2633980907 cycles
HPL_dlaswp06T(  257, 38912, 40960, 38912) 2675164755 cycles
HPL_dlaswp06T(  236, 38401, 40960, 38401) 2424690643 cycles
HPL_dlaswp06T(  236, 38400, 40960, 38400) 2380554742 cycles
HPL_dlaswp06T(  244, 37889, 40960, 37889) 2483956349 cycles
HPL_dlaswp06T(  244, 38400, 40960, 38400) 2498005647 cycles
HPL_dlaswp06T(  259, 37888, 40960, 37888) 2591094842 cycles
HPL_dlaswp06T(  259, 37889, 40960, 37889) 2704659229 cycles
HPL_dlaswp06T(  252, 37888, 40960, 37888) 2523869886 cycles
HPL_dlaswp06T(  252, 37377, 40960, 37377) 2508528112 cycles
HPL_dlaswp06T(  223, 37376, 40960, 37376) 2243484345 cycles
HPL_dlaswp06T(  223, 37377, 40960, 37377) 2269336945 cycles
HPL_dlaswp06T(  254, 36865, 40960, 36865) 2533602673 cycles
HPL_dlaswp06T(  254, 37376, 40960, 37376) 2513782756 cycles
HPL_dlaswp06T(  266, 36864, 40960, 36864) 2624903473 cycles
HPL_dlaswp06T(  266, 36865, 40960, 36865) 2654764593 cycles
HPL_dlaswp06T(  264, 36353, 40960, 36353) 2596242642 cycles
HPL_dlaswp06T(  264, 36864, 40960, 36864) 2560015283 cycles
HPL_dlaswp06T(  238, 36353, 40960, 36353) 2363537032 cycles
HPL_dlaswp06T(  238, 36352, 40960, 36352) 2330636416 cycles
HPL_dlaswp06T(  270, 35841, 40960, 35841) 2559815333 cycles
HPL_dlaswp06T(  270, 36352, 40960, 36352) 2560957840 cycles
HPL_dlaswp06T(  240, 35840, 40960, 35840) 2327989538 cycles
HPL_dlaswp06T(  240, 35841, 40960, 35841) 2346368304 cycles
HPL_dlaswp06T(  254, 35840, 40960, 35840) 2427206482 cycles
HPL_dlaswp06T(  254, 35329, 40960, 35329) 2478372338 cycles
HPL_dlaswp06T(  259, 35328, 40960, 35328) 2399689304 cycles
HPL_dlaswp06T(  259, 35329, 40960, 35329) 2486145223 cycles
HPL_dlaswp06T(  249, 34817, 40960, 34817) 2329386597 cycles
HPL_dlaswp06T(  249, 35328, 40960, 35328) 2380552919 cycles
HPL_dlaswp06T(  272, 34816, 40960, 34816) 2508404817 cycles
HPL_dlaswp06T(  272, 34817, 40960, 34817) 2529980952 cycles
HPL_dlaswp06T(  259, 34305, 40960, 34305) 2408361179 cycles
HPL_dlaswp06T(  259, 34816, 40960, 34816) 2410083089 cycles
HPL_dlaswp06T(  271, 34305, 40960, 34305) 2476415374 cycles
HPL_dlaswp06T(  271, 34304, 40960, 34304) 2487194270 cycles
HPL_dlaswp06T(  257, 33793, 40960, 33793) 2312771978 cycles
HPL_dlaswp06T(  257, 34304, 40960, 34304) 2342802441 cycles
HPL_dlaswp06T(  259, 33793, 40960, 33793) 2350370954 cycles
HPL_dlaswp06T(  259, 33792, 40960, 33792) 2370924033 cycles
HPL_dlaswp06T(  260, 33792, 40960, 33792) 2338393927 cycles
HPL_dlaswp06T(  260, 33281, 40960, 33281) 2312151673 cycles
HPL_dlaswp06T(  255, 33280, 40960, 33280) 2268052440 cycles
HPL_dlaswp06T(  255, 33281, 40960, 33281) 2350170442 cycles
HPL_dlaswp06T(  259, 33280, 40960, 33280) 2263075458 cycles
HPL_dlaswp06T(  259, 32769, 40960, 32769) 2304381545 cycles
HPL_dlaswp06T(  240, 32768, 40960, 32768) 2158038728 cycles
HPL_dlaswp06T(  240, 32769, 40960, 32769) 2148758748 cycles
HPL_dlaswp06T(  241, 32257, 40960, 32257) 2093996825 cycles
HPL_dlaswp06T(  241, 32768, 40960, 32768) 2134951239 cycles
HPL_dlaswp06T(  269, 32257, 40960, 32257) 2348128029 cycles
HPL_dlaswp06T(  269, 32256, 40960, 32256) 2276840573 cycles
HPL_dlaswp06T(  251, 31745, 40960, 31745) 2140978740 cycles
HPL_dlaswp06T(  251, 32256, 40960, 32256) 2131800692 cycles
HPL_dlaswp06T(  256, 31744, 40960, 31744) 2155401220 cycles
HPL_dlaswp06T(  256, 31745, 40960, 31745) 2185165077 cycles
HPL_dlaswp06T(  284, 31744, 40960, 31744) 2365394179 cycles
HPL_dlaswp06T(  284, 31233, 40960, 31233) 2353471992 cycles
HPL_dlaswp06T(  223, 31232, 40960, 31232) 1832168457 cycles
HPL_dlaswp06T(  223, 31233, 40960, 31233) 1899819606 cycles
HPL_dlaswp06T(  251, 31232, 40960, 31232) 2059136575 cycles
HPL_dlaswp06T(  251, 30721, 40960, 30721) 2077440313 cycles
HPL_dlaswp06T(  245, 30721, 40960, 30721) 2016296947 cycles
HPL_dlaswp06T(  245, 30720, 40960, 30720) 2056158406 cycles
HPL_dlaswp06T(  282, 30209, 40960, 30209) 2266188038 cycles
HPL_dlaswp06T(  282, 30720, 40960, 30720) 2300682538 cycles
HPL_dlaswp06T(  255, 30209, 40960, 30209) 2090613938 cycles
HPL_dlaswp06T(  255, 30208, 40960, 30208) 2037542248 cycles
HPL_dlaswp06T(  252, 29697, 40960, 29697) 2011876049 cycles
HPL_dlaswp06T(  252, 30208, 40960, 30208) 2036208476 cycles
HPL_dlaswp06T(  250, 29696, 40960, 29696) 1955184062 cycles
HPL_dlaswp06T(  250, 29697, 40960, 29697) 1954492267 cycles
HPL_dlaswp06T(  267, 29696, 40960, 29696) 2118923296 cycles
HPL_dlaswp06T(  267, 29185, 40960, 29185) 2061874317 cycles
HPL_dlaswp06T(  240, 29184, 40960, 29184) 1829498884 cycles
HPL_dlaswp06T(  240, 29185, 40960, 29185) 1875474993 cycles
HPL_dlaswp06T(  256, 28673, 40960, 28673) 1920647440 cycles
HPL_dlaswp06T(  256, 29184, 40960, 29184) 1964268491 cycles
HPL_dlaswp06T(  245, 28672, 40960, 28672) 1844028938 cycles
HPL_dlaswp06T(  245, 28673, 40960, 28673) 1924091283 cycles
HPL_dlaswp06T(  245, 28161, 40960, 28161) 1863050046 cycles
HPL_dlaswp06T(  245, 28672, 40960, 28672) 1839046273 cycles
HPL_dlaswp06T(  268, 28161, 40960, 28161) 2029407379 cycles
HPL_dlaswp06T(  268, 28160, 40960, 28160) 1983505402 cycles
HPL_dlaswp06T(  241, 27649, 40960, 27649) 1772080109 cycles
HPL_dlaswp06T(  241, 28160, 40960, 28160) 1785123208 cycles
HPL_dlaswp06T(  236, 27648, 40960, 27648) 1744651419 cycles
HPL_dlaswp06T(  236, 27649, 40960, 27649) 1777667766 cycles
HPL_dlaswp06T(  260, 27648, 40960, 27648) 1859646747 cycles
HPL_dlaswp06T(  260, 27137, 40960, 27137) 1925851232 cycles
HPL_dlaswp06T(  279, 27136, 40960, 27136) 1957542502 cycles
HPL_dlaswp06T(  279, 27137, 40960, 27137) 2075788501 cycles
HPL_dlaswp06T(  264, 26625, 40960, 26625) 1865734960 cycles
HPL_dlaswp06T(  264, 27136, 40960, 27136) 1858968402 cycles
HPL_dlaswp06T(  244, 26624, 40960, 26624) 1729815810 cycles
HPL_dlaswp06T(  244, 26625, 40960, 26625) 1757727338 cycles
HPL_dlaswp06T(  259, 26113, 40960, 26113) 1826259144 cycles
HPL_dlaswp06T(  259, 26624, 40960, 26624) 1814888771 cycles
HPL_dlaswp06T(  261, 26113, 40960, 26113) 1829277775 cycles
HPL_dlaswp06T(  261, 26112, 40960, 26112) 1786274448 cycles
HPL_dlaswp06T(  265, 26112, 40960, 26112) 1806803910 cycles
HPL_dlaswp06T(  265, 25601, 40960, 25601) 1801705465 cycles
HPL_dlaswp06T(  271, 25600, 40960, 25600) 1813945047 cycles
HPL_dlaswp06T(  271, 25601, 40960, 25601) 1834152769 cycles
HPL_dlaswp06T(  255, 25600, 40960, 25600) 1671822312 cycles
HPL_dlaswp06T(  255, 25089, 40960, 25089) 1684753565 cycles
HPL_dlaswp06T(  243, 25088, 40960, 25088) 1636484133 cycles
HPL_dlaswp06T(  243, 25089, 40960, 25089) 1636162093 cycles
HPL_dlaswp06T(  246, 24577, 40960, 24577) 1640308896 cycles
HPL_dlaswp06T(  246, 25088, 40960, 25088) 1595639389 cycles
HPL_dlaswp06T(  234, 24576, 40960, 24576) 1547132841 cycles
HPL_dlaswp06T(  234, 24577, 40960, 24577) 1558355128 cycles
HPL_dlaswp06T(  257, 24065, 40960, 24065) 1646333427 cycles
HPL_dlaswp06T(  257, 24576, 40960, 24576) 1670410826 cycles
HPL_dlaswp06T(  240, 24065, 40960, 24065) 1557591382 cycles
HPL_dlaswp06T(  240, 24064, 40960, 24064) 1529590103 cycles
HPL_dlaswp06T(  252, 23553, 40960, 23553) 1601718621 cycles
HPL_dlaswp06T(  252, 24064, 40960, 24064) 1608682188 cycles
HPL_dlaswp06T(  242, 23552, 40960, 23552) 1462748537 cycles
HPL_dlaswp06T(  242, 23553, 40960, 23553) 1521616068 cycles
HPL_dlaswp06T(  249, 23552, 40960, 23552) 1544155399 cycles
HPL_dlaswp06T(  249, 23041, 40960, 23041) 1523571799 cycles
HPL_dlaswp06T(  241, 23040, 40960, 23040) 1470376679 cycles
HPL_dlaswp06T(  241, 23041, 40960, 23041) 1528967824 cycles
HPL_dlaswp06T(  261, 23040, 40960, 23040) 1568343847 cycles
HPL_dlaswp06T(  261, 22529, 40960, 22529) 1599311616 cycles
HPL_dlaswp06T(  260, 22528, 40960, 22528) 1528072153 cycles
HPL_dlaswp06T(  260, 22529, 40960, 22529) 1593479777 cycles
HPL_dlaswp06T(  262, 22017, 40960, 22017) 1557564769 cycles
HPL_dlaswp06T(  262, 22528, 40960, 22528) 1543232294 cycles
HPL_dlaswp06T(  234, 22017, 40960, 22017) 1406166161 cycles
HPL_dlaswp06T(  234, 22016, 40960, 22016) 1355352676 cycles
HPL_dlaswp06T(  264, 21505, 40960, 21505) 1488000421 cycles
HPL_dlaswp06T(  264, 22016, 40960, 22016) 1531447608 cycles
HPL_dlaswp06T(  240, 21504, 40960, 21504) 1346225694 cycles
HPL_dlaswp06T(  240, 21505, 40960, 21505) 1374707144 cycles
HPL_dlaswp06T(  247, 21504, 40960, 21504) 1345250875 cycles
HPL_dlaswp06T(  247, 20993, 40960, 20993) 1391557129 cycles
HPL_dlaswp06T(  248, 20992, 40960, 20992) 1351653985 cycles
HPL_dlaswp06T(  248, 20993, 40960, 20993) 1397348929 cycles
HPL_dlaswp06T(  259, 20481, 40960, 20481) 1400914946 cycles
HPL_dlaswp06T(  259, 20992, 40960, 20992) 1441746565 cycles
HPL_dlaswp06T(  237, 20480, 40960, 20480) 1276594121 cycles
HPL_dlaswp06T(  237, 20481, 40960, 20481) 1262844266 cycles
HPL_dlaswp06T(  241, 19969, 40960, 19969) 1308016720 cycles
HPL_dlaswp06T(  241, 20480, 40960, 20480) 1316594259 cycles
HPL_dlaswp06T(  246, 19969, 40960, 19969) 1334231650 cycles
HPL_dlaswp06T(  246, 19968, 40960, 19968) 1279944053 cycles
HPL_dlaswp06T(  267, 19457, 40960, 19457) 1396867148 cycles
HPL_dlaswp06T(  267, 19968, 40960, 19968) 1388014927 cycles
HPL_dlaswp06T(  242, 19456, 40960, 19456) 1240135740 cycles
HPL_dlaswp06T(  242, 19457, 40960, 19457) 1272369482 cycles
HPL_dlaswp06T(  259, 19456, 40960, 19456) 1304787799 cycles
HPL_dlaswp06T(  259, 18945, 40960, 18945) 1307408339 cycles
HPL_dlaswp06T(  239, 18944, 40960, 18944) 1184525332 cycles
HPL_dlaswp06T(  239, 18945, 40960, 18945) 1181418176 cycles
HPL_dlaswp06T(  265, 18433, 40960, 18433) 1295873516 cycles
HPL_dlaswp06T(  265, 18944, 40960, 18944) 1283296173 cycles
HPL_dlaswp06T(  272, 18432, 40960, 18432) 1291489199 cycles
HPL_dlaswp06T(  272, 18433, 40960, 18433) 1333920010 cycles
HPL_dlaswp06T(  257, 17921, 40960, 17921) 1203168882 cycles
HPL_dlaswp06T(  257, 18432, 40960, 18432) 1247936256 cycles
HPL_dlaswp06T(  255, 17921, 40960, 17921) 1237005861 cycles
HPL_dlaswp06T(  255, 17920, 40960, 17920) 1199990567 cycles
HPL_dlaswp06T(  259, 17409, 40960, 17409) 1215510591 cycles
HPL_dlaswp06T(  259, 17920, 40960, 17920) 1220798190 cycles
HPL_dlaswp06T(  253, 17408, 40960, 17408) 1163338535 cycles
HPL_dlaswp06T(  253, 17409, 40960, 17409) 1196007979 cycles
HPL_dlaswp06T(  248, 17408, 40960, 17408) 1136810401 cycles
HPL_dlaswp06T(  248, 16897, 40960, 16897) 1112016934 cycles
HPL_dlaswp06T(  259, 16896, 40960, 16896) 1121042796 cycles
HPL_dlaswp06T(  259, 16897, 40960, 16897) 1165110870 cycles
HPL_dlaswp06T(  254, 16896, 40960, 16896) 1111955358 cycles
HPL_dlaswp06T(  254, 16385, 40960, 16385) 1127952564 cycles
HPL_dlaswp06T(  254, 16384, 40960, 16384) 1068270168 cycles
HPL_dlaswp06T(  254, 16385, 40960, 16385) 1086961006 cycles
HPL_dlaswp06T(  252, 15873, 40960, 15873) 1037195586 cycles
HPL_dlaswp06T(  252, 16384, 40960, 16384) 1084464540 cycles
HPL_dlaswp06T(  236, 15873, 40960, 15873) 981968413 cycles
HPL_dlaswp06T(  236, 15872, 40960, 15872) 954383383 cycles
HPL_dlaswp06T(  261, 15361, 40960, 15361) 1053867237 cycles
HPL_dlaswp06T(  261, 15872, 40960, 15872) 1092233101 cycles
HPL_dlaswp06T(  255, 15360, 40960, 15360) 1001901012 cycles
HPL_dlaswp06T(  255, 15361, 40960, 15361) 1055029154 cycles
HPL_dlaswp06T(  251, 15360, 40960, 15360) 1009972688 cycles
HPL_dlaswp06T(  251, 14849, 40960, 14849) 998567462 cycles
HPL_dlaswp06T(  255, 14848, 40960, 14848) 977007190 cycles
HPL_dlaswp06T(  255, 14849, 40960, 14849) 984726358 cycles
HPL_dlaswp06T(  253, 14337, 40960, 14337) 968942242 cycles
HPL_dlaswp06T(  253, 14848, 40960, 14848) 967498874 cycles
HPL_dlaswp06T(  251, 14336, 40960, 14336) 918910560 cycles
HPL_dlaswp06T(  251, 14337, 40960, 14337) 927193680 cycles
HPL_dlaswp06T(  252, 13825, 40960, 13825) 915280814 cycles
HPL_dlaswp06T(  252, 14336, 40960, 14336) 934300839 cycles
HPL_dlaswp06T(  276, 13825, 40960, 13825) 1010834572 cycles
HPL_dlaswp06T(  276, 13824, 40960, 13824) 974132141 cycles
HPL_dlaswp06T(  242, 13824, 40960, 13824) 855271754 cycles
HPL_dlaswp06T(  242, 13313, 40960, 13313) 859185496 cycles
HPL_dlaswp06T(  234, 13312, 40960, 13312) 797439564 cycles
HPL_dlaswp06T(  234, 13313, 40960, 13313) 841443310 cycles
HPL_dlaswp06T(  260, 13312, 40960, 13312) 895240519 cycles
HPL_dlaswp06T(  260, 12801, 40960, 12801) 884505125 cycles
HPL_dlaswp06T(  244, 12800, 40960, 12800) 796722682 cycles
HPL_dlaswp06T(  244, 12801, 40960, 12801) 842990516 cycles
HPL_dlaswp06T(  246, 12289, 40960, 12289) 787598212 cycles
HPL_dlaswp06T(  246, 12800, 40960, 12800) 810920909 cycles
HPL_dlaswp06T(  243, 12288, 40960, 12288) 771942383 cycles
HPL_dlaswp06T(  243, 12289, 40960, 12289) 799651747 cycles
HPL_dlaswp06T(  264, 11777, 40960, 11777) 810863354 cycles
HPL_dlaswp06T(  264, 12288, 40960, 12288) 833828671 cycles
HPL_dlaswp06T(  225, 11777, 40960, 11777) 695947928 cycles
HPL_dlaswp06T(  225, 11776, 40960, 11776) 662105097 cycles
HPL_dlaswp06T(  251, 11265, 40960, 11265) 736946251 cycles
HPL_dlaswp06T(  251, 11776, 40960, 11776) 745377482 cycles
HPL_dlaswp06T(  247, 11264, 40960, 11264) 699807609 cycles
HPL_dlaswp06T(  247, 11265, 40960, 11265) 726902346 cycles
HPL_dlaswp06T(  268, 11264, 40960, 11264) 745104406 cycles
HPL_dlaswp06T(  268, 10753, 40960, 10753) 755032057 cycles
HPL_dlaswp06T(  251, 10752, 40960, 10752) 661959913 cycles
HPL_dlaswp06T(  251, 10753, 40960, 10753) 702616642 cycles
HPL_dlaswp06T(  247, 10241, 40960, 10241) 647134123 cycles
HPL_dlaswp06T(  247, 10752, 40960, 10752) 663555790 cycles
HPL_dlaswp06T(  242, 10240, 40960, 10240) 621968424 cycles
HPL_dlaswp06T(  242, 10241, 40960, 10241) 647912310 cycles
HPL_dlaswp06T(  252,  9729, 40960,  9729) 643357460 cycles
HPL_dlaswp06T(  252, 10240, 40960, 10240) 668343431 cycles
HPL_dlaswp06T(  247,  9728, 40960,  9728) 595556118 cycles
HPL_dlaswp06T(  247,  9729, 40960,  9729) 636797460 cycles
HPL_dlaswp06T(  249,  9217, 40960,  9217) 586812770 cycles
HPL_dlaswp06T(  249,  9728, 40960,  9728) 586773529 cycles
HPL_dlaswp06T(  268,  9216, 40960,  9216) 593971347 cycles
HPL_dlaswp06T(  268,  9217, 40960,  9217) 631977495 cycles
HPL_dlaswp06T(  264,  9216, 40960,  9216) 598423475 cycles
HPL_dlaswp06T(  264,  8705, 40960,  8705) 598175520 cycles
HPL_dlaswp06T(  222,  8704, 40960,  8704) 502648411 cycles
HPL_dlaswp06T(  222,  8705, 40960,  8705) 509036654 cycles
HPL_dlaswp06T(  254,  8193, 40960,  8193) 533074591 cycles
HPL_dlaswp06T(  254,  8704, 40960,  8704) 550433682 cycles
HPL_dlaswp06T(  253,  8192, 40960,  8192) 507036793 cycles
HPL_dlaswp06T(  253,  8193, 40960,  8193) 532903445 cycles
HPL_dlaswp06T(  274,  7681, 40960,  7681) 539133785 cycles
HPL_dlaswp06T(  274,  8192, 40960,  8192) 554918756 cycles
HPL_dlaswp06T(  241,  7680, 40960,  7680) 455987330 cycles
HPL_dlaswp06T(  241,  7681, 40960,  7681) 478484175 cycles
HPL_dlaswp06T(  245,  7169, 40960,  7169) 446460856 cycles
HPL_dlaswp06T(  245,  7680, 40960,  7680) 457986584 cycles
HPL_dlaswp06T(  249,  7168, 40960,  7168) 433447545 cycles
HPL_dlaswp06T(  249,  7169, 40960,  7169) 467703684 cycles
HPL_dlaswp06T(  244,  7168, 40960,  7168) 421554238 cycles
HPL_dlaswp06T(  244,  6657, 40960,  6657) 409749204 cycles
HPL_dlaswp06T(  244,  6656, 40960,  6656) 394589292 cycles
HPL_dlaswp06T(  244,  6657, 40960,  6657) 415194818 cycles
HPL_dlaswp06T(  263,  6145, 40960,  6145) 388730466 cycles
HPL_dlaswp06T(  263,  6656, 40960,  6656) 415073790 cycles
HPL_dlaswp06T(  245,  6144, 40960,  6144) 354897745 cycles
HPL_dlaswp06T(  245,  6145, 40960,  6145) 378024708 cycles
HPL_dlaswp06T(  264,  5633, 40960,  5633) 373753544 cycles
HPL_dlaswp06T(  264,  6144, 40960,  6144) 393102910 cycles
HPL_dlaswp06T(  255,  5632, 40960,  5632) 319982029 cycles
HPL_dlaswp06T(  255,  5633, 40960,  5633) 360426758 cycles
HPL_dlaswp06T(  260,  5121, 40960,  5121) 338536317 cycles
HPL_dlaswp06T(  260,  5632, 40960,  5632) 359189339 cycles
HPL_dlaswp06T(  256,  5120, 40960,  5120) 291622678 cycles
HPL_dlaswp06T(  256,  5121, 40960,  5121) 324827090 cycles
HPL_dlaswp06T(  254,  4609, 40960,  4609) 292388966 cycles
HPL_dlaswp06T(  254,  5120, 40960,  5120) 295056838 cycles
HPL_dlaswp06T(  224,  4608, 40960,  4608) 247297363 cycles
HPL_dlaswp06T(  224,  4609, 40960,  4609) 272497155 cycles
HPL_dlaswp06T(  251,  4097, 40960,  4097) 256996894 cycles
HPL_dlaswp06T(  251,  4608, 40960,  4608) 275428640 cycles
HPL_dlaswp06T(  243,  4096, 40960,  4096) 233434947 cycles
HPL_dlaswp06T(  243,  4097, 40960,  4097) 254398015 cycles
HPL_dlaswp06T(  253,  3585, 40960,  3585) 217285781 cycles
HPL_dlaswp06T(  253,  4096, 40960,  4096) 237558692 cycles
HPL_dlaswp06T(  253,  3584, 40960,  3584) 201815483 cycles
HPL_dlaswp06T(  253,  3585, 40960,  3585) 223223990 cycles
HPL_dlaswp06T(  248,  3073, 40960,  3073) 183299424 cycles
HPL_dlaswp06T(  248,  3584, 40960,  3584) 196250313 cycles
HPL_dlaswp06T(  237,  3072, 40960,  3072) 155436817 cycles
HPL_dlaswp06T(  237,  3073, 40960,  3073) 174887676 cycles
HPL_dlaswp06T(  262,  2561, 40960,  2561) 149872433 cycles
HPL_dlaswp06T(  262,  3072, 40960,  3072) 159653558 cycles
HPL_dlaswp06T(  217,  2560, 40960,  2560) 120927448 cycles
HPL_dlaswp06T(  217,  2561, 40960,  2561) 131182043 cycles
HPL_dlaswp06T(  252,  2049, 40960,  2049) 117143688 cycles
HPL_dlaswp06T(  252,  2560, 40960,  2560) 141379131 cycles
HPL_dlaswp06T(  237,  2048, 40960,  2048) 98245882 cycles
HPL_dlaswp06T(  237,  2049, 40960,  2049) 107072959 cycles
HPL_dlaswp06T(  257,  1537, 40960,  1537) 83647461 cycles
HPL_dlaswp06T(  257,  2048, 40960,  2048) 108189319 cycles
HPL_dlaswp06T(  230,  1536, 40960,  1536) 61084856 cycles
HPL_dlaswp06T(  230,  1537, 40960,  1537) 68536528 cycles
HPL_dlaswp06T(  255,  1025, 40960,  1025) 50153576 cycles
HPL_dlaswp06T(  255,  1536, 40960,  1536) 70727002 cycles
HPL_dlaswp06T(  211,  1024, 40960,  1024) 32148076 cycles
HPL_dlaswp06T(  211,  1025, 40960,  1025) 39014981 cycles
HPL_dlaswp06T(  252,   513, 40960,   513) 21220935 cycles
HPL_dlaswp06T(  252,  1024, 40960,  1024) 37582438 cycles
HPL_dlaswp06T(  171,   512, 40960,   512) 11876023 cycles
HPL_dlaswp06T(  171,   513, 40960,   513) 12655857 cycles
HPL_dlaswp06T(  272,   512, 40960,   512) 16224702 cycles
HPL_dlaswp06T: 110902526005
HPL_dlaswp06T: 108769729258
HPL_dlaswp06T: 106189048145
HPL_dlaswp06T: 111609275039
 */

class dlaswp06T_impl
{
    private:
        const size_t M;
        const size_t LDA;
        const size_t LDU;
        double *__restrict__ const A;
        double *__restrict__ const U;
        const int *__restrict__ const LINDXA;
    public:
        dlaswp06T_impl(size_t _M, double *_A, size_t _LDA, double *_U, size_t _LDU,
                const int *_LINDXA)
            : M(_M), LDA(_LDA), LDU(_LDU), A(_A), U(_U), LINDXA(_LINDXA) {}
        void operator()(const tbb::blocked_range<size_t> &range) const
        {
            const size_t begin = range.begin();
            const size_t columns = range.end() - begin;
            double *u = &U[begin];
            double *uNext = u + LDU;
            double *A2 = &A[begin * LDA];
            for (size_t i = 0; i < M; ++i) {
                double *a = &A2[LINDXA[i]];
                ptrdiff_t aNext = &A2[LINDXA[i + 1]] - a;
                size_t j = 7;
                for (; j < columns; j += 8) {
                    _m_prefetchw(&uNext[j - 7]);
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 7]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 6]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 5]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 4]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 3]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 2]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j - 1]); a += LDA;
                    _m_prefetchw(&a[aNext]); swap(*a, u[j    ]); a += LDA;
                }
                _m_prefetchw(&uNext[j - 7]);
                double *u2 = &u[columns];
                j = (j - columns) * 0x15;
                long tmp1;
                asm volatile(
                    "lea 0x6(%%rip), %[tmp1]\n"
                    "lea (%[tmp1], %[tmp0], 1), %[tmp1]\n"
                    "jmpq *%[tmp1]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x38(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x38(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x30(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x30(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x28(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x28(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x20(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x20(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x18(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x18(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x10(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x10(%[u])\n"
                    "add %[LDA], %[a]\n"

                    "prefetchw (%[a],%[aNext],1)\n"
                    "mov (%[a]), %[tmp0]\n"
                    "mov -0x8(%[u]), %[tmp1]\n"
                    "mov %[tmp1], (%[a])\n"
                    "mov %[tmp0], -0x8(%[u])\n"
                    "add %[LDA], %[a]\n"
                    : [a]"+d"(a),
                      [tmp0]"+S"(j),
                      [tmp1]"=a"(tmp1)
                    : [u]"D"(u2),
                      [aNext]"c"(aNext * 8),
                      [LDA]"b"(LDA * 8)
                   );
                u = uNext;
                uNext += LDU;
            }
        }
};

extern "C" void HPL_dlaswp06T(const int M, const int N, double *A,
        const int LDA, double *U, const int LDU, const int *LINDXA)
{
#ifdef TRACE_CALLS
   uint64_t tr_start, tr_end, tr_diff;
   tr_start = util_getTimestamp();
#endif /* TRACE_CALLS */

#ifdef USE_ORIGINAL_LASWP
#include "HPL_dlaswp06T.c"
#else
    if (M <= 0 || N <= 0) {
        return;
    }

    tbb::parallel_for (tbb::blocked_range<size_t>(0, N, 32),
            dlaswp06T_impl(M, A, LDA, U, LDU, LINDXA));
#endif

#ifdef TRACE_CALLS
   tr_end = util_getTimestamp();
   tr_diff = util_getTimeDifference( tr_start, tr_end );

   fprintf( trace_dgemm, "DLASWP06T,M=%i,N=%i,LDA=%i,LDU=%i,TIME=%lu,THRPT=%.2fGB/s\n", M, N, LDA, LDU, tr_diff,
           0.004 * sizeof(double) * M * N / tr_diff );
#ifdef TRACE_PERMDATA
   char filename[256];
   snprintf(filename, 256, "dlaswp06T.%04d.%05d.%05d.%05d.dat", M, N, LDA, LDU);
   FILE *permdata = fopen(filename, "w");
   fwrite(LINDXA, sizeof(LINDXA[0]), M, permdata);
   fclose(permdata);
#endif
#endif /* TRACE_CALLS */
}
