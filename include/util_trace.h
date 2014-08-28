/**
 * Utility header for tracing some stuff
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
	uint32_t invocations;
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
#ifdef TRACE_LASWP
#include "../caldgemm/cmodules/timer.h"
extern HighResTimer TimerLASWP;
#define LASWP_TIMER_START HighResTimer.ResetStart();
#define LASWP_TIMER_STOP  double laswp_time = HighResTimer.GetCurrentElapsedTime();
#else
#define LASWP_TIMER_START
#define LASWP_TIMER_SOP
#endif

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
clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_cpu); \
	LASWP_TIMER_START

#define END_TRACE \
end_wall = util_getTimestamp(); \
clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_cpu); \
++(counter->invocations); \
counter->walltime += util_getTimeDifference( start_wall, end_wall ); \
counter->cputime += end_cpu.tv_sec * 1000000ull + end_cpu.tv_nsec / 1000ull - start_cpu.tv_sec * 1000000ull - start_cpu.tv_nsec / 1000ull; \
	LASWP_TIMER_STOP

#else

#define START_TRACE( FUNC_NAME ) LASWP_TIMER_START
#define END_TRACE LASWP_TIMER_SOP

#endif /* TRACE_CALLS */

#ifdef __cplusplus
}
#endif

