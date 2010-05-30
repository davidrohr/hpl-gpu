/**
 * A utility for timing stuff
 *
 * (c) Matthias Bach <bach@compeng.uni-frankfurt.de> 2010
 */

#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

#include <stdint.h>

/**
 * Retrieves a timestamt that never walks backwards.
 *
 * Use util_getTimeDifference to compare Timestamps
 */
uint64_t util_getTimestamp( );

/**
 * Retrieves the difference in mus between two timestamps
 */
uint64_t util_getTimeDifference( uint64_t start, uint64_t end );

#endif

