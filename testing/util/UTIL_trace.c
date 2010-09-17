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

/*
 * File pointers to write traces to
 */

#ifdef TRACE_CALLS
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

