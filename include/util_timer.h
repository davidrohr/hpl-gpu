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

#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Retrieves Timestamp in museconds since the Epoch
 * ( can walk backwards )
 */
uint64_t util_getWalltime( );

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

#ifdef __cplusplus
}
#endif

#endif

