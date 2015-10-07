/*
 *  -- High Performance Computing Linpack Benchmark (HPL-GPU)
 *     HPL-GPU - 2.0 - 2015
 *
 *     David Rohr
 *     Matthias Kretz
 *     Matthias Bach
 *     Goethe Universität, Frankfurt am Main
 *     Frankfurt Institute for Advanced Studies
 *     (C) Copyright 2010 All Rights Reserved
 *
 *     Antoine P. Petitet
 *     University of Tennessee, Knoxville
 *     Innovative Computing Laboratory
 *     (C) Copyright 2000-2008 All Rights Reserved
 *
 *  -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 *  1. Redistributions  of  source  code  must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce  the above copyright
 *  notice, this list of conditions,  and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. All  advertising  materials  mentioning  features  or  use of this
 *  software must display the following acknowledgements:
 *  This  product  includes  software  developed  at  the  University  of
 *  Tennessee, Knoxville, Innovative Computing Laboratory.
 *  This product  includes software  developed at the Frankfurt Institute
 *  for Advanced Studies.
 *
 *  4. The name of the  University,  the name of the  Laboratory,  or the
 *  names  of  its  contributors  may  not  be used to endorse or promote
 *  products  derived   from   this  software  without  specific  written
 *  permission.
 *
 *  -- Disclaimer:
 *
 *  THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 *  OR  CONTRIBUTORS  BE  LIABLE FOR ANY  DIRECT,  INDIRECT,  INCIDENTAL,
 *  SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES  (INCLUDING,  BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,  OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * ======================================================================
 */

#ifndef HPL_PMISC_H
#define HPL_PMISC_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_misc.h"
#ifndef HPL_NO_MPI_LIB
#include "mpi.h"
#else
#include "mpi_disabled.h"
#endif

#ifdef VTRACE
#include <vt_user.h>
#include <pthread.h>
extern pthread_mutex_t global_vt_mutex;
#define VT_USER_START_A(a) {pthread_mutex_lock(&global_vt_mutex);VT_USER_START(a);pthread_mutex_unlock(&global_vt_mutex);}
#define VT_USER_END_A(a) {pthread_mutex_lock(&global_vt_mutex);VT_USER_END(a);pthread_mutex_unlock(&global_vt_mutex);}
#else
#define VT_USER_START_A(a)
#define VT_USER_END_A(a)
#endif

#endif
/*
 * End of hpl_pmisc.h
 */
