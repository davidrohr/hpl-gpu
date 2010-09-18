/**
 * A utility for timing stuff
 *
 * The source code is property of the Frankfurt Institute for Advanced Studies
 * (FIAS). None of the material may be copied, reproduced, distributed,
 * republished, downloaded, displayed, posted or transmitted in any form or by
 * any means, including, but not limited to, electronic, mechanical,
 * photocopying, recording, or otherwise, without the prior written permission
 * of FIAS.
 *
 * Authors:
 * David Rohr (drohr@jwdt.org)
 * Matthias Bach (bach@compeng.uni-frankfurt.de)
 * Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 */

#include "util_timer.h"

#ifdef __APPLE__
#include <mach/mach_time.h>

// method to get monotonic mac time, inspired by 
// http://www.wand.net.nz/~smr26/wordpress/2009/01/19/monotonic-time-in-mac-os-x/

#else

#include <time.h>

#endif

/**
 * Retrieves a timestamt that never walks backwards.
 *
 * Use util_getTimeDifference to compare Timestamps
 */
uint64_t util_getTimestamp( )
{
#ifdef __APPLE__
        return mach_absolute_time();
#else
        struct timespec now;
        clock_gettime( CLOCK_MONOTONIC, &now );

        long nanos = now.tv_nsec;
        time_t seconds = now.tv_sec;

        uint64_t mus = ((uint64_t) seconds ) * 1000 * 1000 + nanos / 1000;

        return mus;
#endif
}

/**
 * Retrieves the difference in mus between two timestamps
 */
uint64_t util_getTimeDifference( uint64_t start, uint64_t end )
{
        uint64_t mus;
#ifdef __APPLE__
        uint64_t difference = end - start;
        mach_timebase_info_data_t info = {0,0};

        if (info.denom == 0)
                mach_timebase_info(&info);

        uint64_t nanos = difference * (info.numer / info.denom);
        mus = nanos / 1000;
#else
        mus = end - start;
#endif
        return mus;
}

