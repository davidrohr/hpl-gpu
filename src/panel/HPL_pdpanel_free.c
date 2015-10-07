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

#include "hpl.h"

#include "util_timer.h"
#include "util_trace.h"

#include "util_cal.h"

int HPL_pdpanel_free
(
   HPL_T_panel *                    PANEL
)
{
/* 
 * Purpose
 * =======
 *
 * HPL_pdpanel_free deallocates  the panel resources  and  stores the error
 * code returned by the panel factorization.
 *
 * Arguments
 * =========
 *
 * PANEL   (local input/output)          HPL_T_panel *
 *         On entry,  PANEL  points  to  the  panel data  structure from
 *         which the resources should be deallocated.
 *
 * ---------------------------------------------------------------------
 */ 
START_TRACE( PDPANEL_FREE )

/* ..
 * .. Executable Statements ..
 */
   if( PANEL->pmat->info == 0 ) PANEL->pmat->info = *(PANEL->DINFO);

   //if( PANEL->WORK  ) free( PANEL->WORK  );
   //if( PANEL->IWORK ) free( PANEL->IWORK );

END_TRACE

   return( MPI_SUCCESS );
/*
 * End of HPL_pdpanel_free
 */
}
