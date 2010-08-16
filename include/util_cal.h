/**
 * A wrapper for the C++ CALDGEMM
 *
 * (c) Matthias Bach <bach@compeng.uni-frankfurt.de> 2010
 */

#ifndef UTIL_CAL_H
#define UTIL_CAL_H

#ifdef __cplusplus
extern "C"
{
#endif

void CALDGEMM_Init();
void CALDGEMM_Shutdown();
void CALDGEMM_dgemm( const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE,
                     const enum CBLAS_TRANSPOSE, const int, const int,
   const int,       const double,    const double *,  const int,
   const double *,  const int,       const double,    double *,
   const int );

#ifdef __cplusplus
}
#endif

#endif


