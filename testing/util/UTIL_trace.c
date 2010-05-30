/**
 * Utility header for tracing some stuff
 */

#include "util_trace.h"

#include <string.h>

#include "hpl.h"

/*
 * File pointers to write traces to
 */

#ifdef TRACE_DGEMM
FILE* trace_dgemm = 0;
#endif

FILE* openTraceFile( const char *basename, const int run, const int rank )
{
	FILE* res;

	int baseSize = strlen( basename );
	const int rankSize = 5;
	const int addSize = 6;

	char* filename = malloc( baseSize + 2 * rankSize + addSize );
	if( ! filename )
		HPL_abort( __LINE__, __FILE__, "Failed to allocate mem for filename generation" );

	int printed = sprintf( filename, "%s.%.5d.%.5d.csv", basename, run, rank );
	if( ! printed )
		HPL_abort( __LINE__, __FILE__, "Failed to generate filename" );

	res = fopen( filename, "w" );
	if( ! res )
		HPL_abort( __LINE__, __FILE__, "Failed to open trace file" );

	free( filename );

	return res;
}

