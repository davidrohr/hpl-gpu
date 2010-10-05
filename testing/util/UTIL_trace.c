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

#include "util_trace.h"

#include <string.h>

#include "hpl.h"

FILE* openTraceFile( const char *basename, const int run, const int rank );

/**
 * Definition of the counters.
 *
 * Having this dynamically created would be nice, but they are small and I am
 * a lazy person ...
 */
#define MAX_COUNTERS 100
size_t used_counters = 0;
trace_counters_t counters[ MAX_COUNTERS ];

/**
 * Resets all trace counters.
 */
void resetTraceCounters() 
{
	size_t i;
	for( i = 0; i < used_counters; ++i )
	{
		trace_counters_t * counter = &(counters[i]);
		counter->walltime = 0;
		counter->cputime = 0;
		counter->invocations = 0;
	}
}

void writeTraceCounters( const char *basename, const int run, const int rank )
{
	float ratio;
	FILE* tracefile = openTraceFile( basename, run, rank );

	fprintf( tracefile, "#FUNCTION NAME,INVOCATIONS,WALLTIME,CPUTIME,RATIO CPU vs. WALL\n" );

	size_t i;
	for( i = 0; i < used_counters; ++i )
	{
		trace_counters_t * counter = &(counters[i]);
		ratio = (float) counter->cputime / (float) counter->walltime;
		fprintf( tracefile, "%s,%u,%lu,%lu,%.2f\n", counter->func_name, counter->invocations, counter->walltime, counter->cputime, ratio );
	}

	fclose( tracefile );
}

void releaseTraceCounters()
{
	// we are using statically allocated memory, so we don't really have to do something
	// at least we reset the counter, this way maybe nobody notices
	used_counters = 0;
}

trace_counters_t * aquireTraceCounter()
{
	if( used_counters >= MAX_COUNTERS )
		HPL_pabort( __LINE__, "aquireTraceCounter", "Ran out of counters for tracing. Increase the variable MAX_COUNTERS in %s and make sure you are properly resusing your trace counters.", __FILE__ );

	trace_counters_t * new_counter = &counters[used_counters];
	++used_counters;

	new_counter->func_name = "";
	new_counter->walltime = 0;
	new_counter->cputime = 0;
	new_counter->invocations = 0;

	return new_counter;
}

FILE* openTraceFile( const char *basename, const int run, const int rank )
{
    FILE* res;

    int baseSize = strlen( basename );
    const int rankSize = 5;
    const int addSize = 6;

    char* filename = malloc( baseSize + 2 * rankSize + addSize );
    if( ! filename )
        HPL_pabort( __LINE__, openTraceFile, "Failed to allocate mem for filename generation" );

    int printed = sprintf( filename, "%s.%.5d.%.5d.csv", basename, run, rank );
    if( ! printed )
        HPL_pabort( __LINE__, openTraceFile, "Failed to generate filename" );

    res = fopen( filename, "w" );
    if( ! res )
        HPL_pabort( __LINE__, openTraceFile, "Failed to open trace file" );

    free( filename );

    return res;
}

