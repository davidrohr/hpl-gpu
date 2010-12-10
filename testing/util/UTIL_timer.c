/**
 * A utility for timing stuff
 *
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
 * Use util_getTimeDifference to compare Walltimes
 */
uint64_t util_getWalltime( )
{
        struct timespec now;
        clock_gettime( CLOCK_REALTIME, &now );

        long nanos = now.tv_nsec;
        time_t seconds = now.tv_sec;

        uint64_t mus = ((uint64_t) seconds ) * 1000 * 1000 + nanos / 1000;

        return mus;
}

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

