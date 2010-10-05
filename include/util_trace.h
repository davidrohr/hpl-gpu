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

#include <util_timer.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Structure representing the trace values of one function
 */
typedef struct trace_counters
{
	uint64_t walltime;
	uint64_t cputime;
	char * func_name;
} trace_counters_t;

/**
 * Resets all trace counters.
 */
void resetTraceCounters();

/**
 * Writes all trace counters to a file whose name is derived from the arguments.
 */
void writeTraceCounters( const char *basename, const int run, const int rank );

/**
 * Cleans up memory if trace counters are no longer needed.
 *
 * Be aware that this invalidates all pointers that have been handed out via
 * createTraceCounter() and there is pretty much no way to notice this
 * invalidation for anybody that called createTraceCounter();
 * Therefore only call it if you KNOW that no trace counters will be used
 * anymore, e.g. right before exit(..).
 */
void releaseTraceCounters();

/**
 * Aquire a new trace counter. This will be globally cleared by the
 * releaseTraceCounters function.
 */
trace_counters_t * aquireTraceCounter();

/**
 * Utility macro to be used at the start of the function to start the trace
 */
#ifdef TRACE_CALLS

#define START_TRACE( FUNC_NAME ) \
uint64_t start_wall, end_wall; \
struct timespec start_cpu, end_cpu; \
static trace_counters_t * counter = 0; \
if( counter == 0 ) \
{ \
	counter = aquireTraceCounter(); \
	counter->func_name = #FUNC_NAME; \
} \
start_wall = util_getTimestamp(); \
clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_cpu);

#define END_TRACE \
end_wall = util_getTimestamp(); \
clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_cpu); \
counter->walltime += util_getTimeDifference( start_wall, end_wall ); \
counter->cputime += end_cpu.tv_sec * 1000000ull + end_cpu.tv_nsec / 1000ull - start_cpu.tv_sec * 1000000ull - start_cpu.tv_nsec / 1000ull;

#else

#define START_TRACE( FUNC_NAME )
#define END_TRACE

#endif /* TRACE_CALLS */

#ifdef __cplusplus
}
#endif

