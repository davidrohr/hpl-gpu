/**
 * A wrapper for the C++ CALDGEMM
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

#ifndef UTIL_CAL_H
#define UTIL_CAL_H

#ifdef __cplusplus
extern "C"
{
#endif

int CALDGEMM_Init();
void CALDGEMM_Shutdown();
void CALDGEMM_dgemm( const enum CBLAS_ORDER, const enum CBLAS_TRANSPOSE,
                     const enum CBLAS_TRANSPOSE, const int, const int,
   const int,       const double,    const double *,  const int,
   const double *,  const int,       const double,    double *,
   const int, int LinpackCallbacks );
void* CALDGEMM_alloc(size_t size);
void CALDGEMM_free(void* ptr);
void CALDGEMM_set_num_nodes(int num);

#ifdef __cplusplus
}
#endif

#endif


