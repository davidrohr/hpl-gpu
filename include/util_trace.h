/**
 * Utility header for tracing some stuff
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

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * File pointers to write traces to
 */

#ifdef TRACE_CALLS
extern FILE* trace_dgemm;
#endif

FILE* openTraceFile( const char *basename, const int run, const int rank );

#ifdef __cplusplus
}
#endif

