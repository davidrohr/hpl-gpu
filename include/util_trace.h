/**
 * Utility header for tracing some stuff
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

