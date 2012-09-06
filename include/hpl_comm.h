/*
 *  -- High Performance Computing Linpack Benchmark (HPL-GPU)
 *     HPL-GPU - 1.1 - 2011
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

#ifndef HPL_COMM_H
#define HPL_COMM_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_pmisc.h"
#include "hpl_panel.h"
/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */
typedef enum
{
 HPL_1RING = 401, /* Increasing ring */
 HPL_1RING_M = 402, /* Increasing ring (modified) */
 HPL_2RING = 403, /* Increasing 2-ring */
 HPL_2RING_M = 404, /* Increasing 2-ring (modified) */
 HPL_BLONG = 405, /* long broadcast */
 HPL_BLONG_M = 406, /* long broadcast (modified) */
 HPL_MPI_BCAST = 407 /* Use MPI_BCAST with subcommunicators */
} HPL_T_TOP;
/*
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define HPL_FAILURE 0
#define HPL_SUCCESS 1
#define HPL_KEEP_TESTING 2
/*
 * ---------------------------------------------------------------------
 * comm function prototypes
 * ---------------------------------------------------------------------
 */
int HPL_send( double *, int, int, int, MPI_Comm );
int HPL_recv( double *, int, int, int, MPI_Comm );
int HPL_sdrv( double *, int, int, double *, int, int, int, MPI_Comm );
int HPL_binit( HPL_T_panel *);
int HPL_bcast( HPL_T_panel *);
int HPL_packL( HPL_T_panel *, const int, const int, const int );
void HPL_copyL( HPL_T_panel * );
 
int HPL_binit_1ring ( HPL_T_panel * );
int HPL_bcast_1ring ( HPL_T_panel *);
 
int HPL_binit_1rinM ( HPL_T_panel * );
int HPL_bcast_1rinM ( HPL_T_panel *);
 
int HPL_binit_2ring ( HPL_T_panel * );
int HPL_bcast_2ring ( HPL_T_panel *);
 
int HPL_binit_2rinM ( HPL_T_panel * );
int HPL_bcast_2rinM ( HPL_T_panel *);
 
int HPL_binit_blong ( HPL_T_panel * );
int HPL_bcast_blong ( HPL_T_panel *);
 
int HPL_binit_blonM ( HPL_T_panel * );
int HPL_bcast_blonM ( HPL_T_panel *);

int HPL_binit_mpi ( HPL_T_panel * );
int HPL_bcast_mpi ( HPL_T_panel *);

#endif
/*
 * End of hpl_comm.h
 */
