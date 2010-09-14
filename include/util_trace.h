/**
 * Utility header for tracing some stuff
 */

#include <stdio.h>

/*
 * File pointers to write traces to
 */

#ifdef TRACE_BLAS
extern FILE* trace_dgemm;
#endif

FILE* openTraceFile( const char *basename, const int run, const int rank );

