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

/*
 * Include files
 */
#include "hpl.h"

#include "util_timer.h"
#include "util_trace.h"

#define   I_SEND    0
#define   I_RECV    1

void HPL_rollT
(
   HPL_T_panel *                    PANEL,
   const int                        N,
   double *                         U,
   const int                        LDU,
   const int *                      IPLEN,
   const int *                      IPMAP,
   const int *                      IPMAPM1
)
{
/* 
 * Purpose
 * =======
 *
 * HPL_rollT rolls the local arrays containing the local pieces of U, so
 * that on exit to this function  U  is replicated in every process row.
 * Arguments
 * =========
 *
 * PANEL   (local input/output)          HPL_T_panel *
 *         On entry,  PANEL  points to the data structure containing the
 *         panel (to be rolled) information.
 *
 * N       (local input)                 const int
 *         On entry, N specifies the local number of rows of  U.  N must
 *         be at least zero.
 *
 * U       (local input/output)          double *
 *         On entry,  U  is an array of dimension (LDU,*) containing the
 *         local pieces of U in each process row.
 *
 * LDU     (local input)                 const int
 *         On entry, LDU specifies the local leading dimension of U. LDU
 *         should be at least  MAX(1,N).
 *
 * IPLEN   (global input)                const int *
 *         On entry, IPLEN is an array of dimension NPROW+1.  This array
 *         is such that IPLEN[i+1] - IPLEN[i] is the number of rows of U
 *         in each process row.
 *
 * IPMAP   (global input)                const int *
 *         On entry, IMAP  is an array of dimension  NPROW.  This  array
 *         contains  the  logarithmic mapping of the processes. In other
 *         words,  IMAP[myrow]  is the absolute coordinate of the sorted
 *         process.
 *
 * IPMAPM1 (global input)                const int *
 *         On entry,  IMAPM1  is an array of dimension NPROW. This array
 *         contains  the inverse of the logarithmic mapping contained in
 *         IMAP: For i in [0.. NPROW) IMAPM1[IMAP[i]] = i.
 *
 * ---------------------------------------------------------------------
 */ 
START_TRACE( ROLLT )

/*
 * .. Local Variables ..
 */
#ifndef HPL_SEND_U_PADDING
   MPI_Datatype               type[2];
#endif
   MPI_Comm                   comm;
   int                        Cmsgid=MSGID_BEGIN_PFACT, ibufR, ibufS,
                              ierr=MPI_SUCCESS, il, k, l, lengthR, 
                              lengthS, mydist, myrow, next, npm1, nprow,
                              partner, prev;
/* ..
 * .. Executable Statements ..
 */
   if( N <= 0 ) return;

   npm1 = ( nprow = PANEL->grid->nprow ) - 1; myrow = PANEL->grid->myrow;
   comm = PANEL->grid->col_comm;
/*
 * Rolling phase
 */
   mydist = IPMAPM1[myrow];
   prev   = IPMAP[MModSub1( mydist, nprow )];
   next   = IPMAP[MModAdd1( mydist, nprow )];
 
   for( k = 0; k < npm1; k++ )
   {
      l = (int)( (unsigned int)(k) >> 1 );
 
      if( ( ( mydist + k ) & 1 ) != 0 )
      {
         il      = MModAdd( mydist, l,   nprow );
         lengthS = IPLEN[il+1] - ( ibufS = IPLEN[il] );
         il    = MModSub( mydist, l+1, nprow );
         lengthR = IPLEN[il+1] - ( ibufR = IPLEN[il] ); partner = prev;
      }
      else
      {
         il    = MModSub( mydist, l,   nprow );
         lengthS = IPLEN[il+1] - ( ibufS = IPLEN[il] );
         il    = MModAdd( mydist, l+1, nprow );
         lengthR = IPLEN[il+1] - ( ibufR = IPLEN[il] ); partner = next;
      }

	  //Create Datatypes if necesary
#ifndef HPL_SEND_U_PADDING
      if( lengthR > 0 )
      {
         if( LDU == N ) {
            MPI_Type_contiguous( lengthR * LDU, MPI_DOUBLE, &type[I_RECV] );
         } else {
            MPI_Type_vector( lengthR, N, LDU, MPI_DOUBLE, &type[I_RECV] );
         }
         MPI_Type_commit( &type[I_RECV] );
	  }
      if( lengthS > 0 )
      {
         if( LDU == N ) {
            ierr =   MPI_Type_contiguous( lengthS*LDU, MPI_DOUBLE, &type[I_SEND] );
         } else {
            ierr =   MPI_Type_vector( lengthS, N, LDU, MPI_DOUBLE, &type[I_SEND] );
         }
         MPI_Type_commit( &type[I_SEND] );
	  }
#endif
 
	  //Communication
	  if (lengthR > 0 && lengthS > 0)
	  {
#ifndef HPL_SEND_U_PADDING
         MPI_Sendrecv_Mod( Mptr( U, 0, ibufS, LDU ), 1, type[I_SEND], partner, Cmsgid, Mptr( U, 0, ibufR, LDU ), 1, type[I_RECV], partner, Cmsgid, comm, MPI_STATUS_IGNORE );
#else
         MPI_Sendrecv_Mod( Mptr( U, 0, ibufS, LDU ), lengthS*LDU, MPI_DOUBLE, partner, Cmsgid, Mptr( U, 0, ibufR, LDU ), lengthR*LDU, MPI_DOUBLE, partner, Cmsgid, comm, MPI_STATUS_IGNORE );
#endif
	  }
	  else if (lengthR > 0)
	  {
#ifndef HPL_SEND_U_PADDING
         MPI_Recv_Mod( Mptr( U, 0, ibufR, LDU ), 1, type[I_RECV], partner, Cmsgid, comm, MPI_STATUS_IGNORE );
#else
         MPI_Recv_Mod( Mptr( U, 0, ibufR, LDU ), lengthR*LDU, MPI_DOUBLE, partner, Cmsgid, comm, MPI_STATUS_IGNORE );
#endif
	  }
	  else if (lengthS > 0)
	  {
#ifndef HPL_SEND_U_PADDING
         MPI_Send_Mod( Mptr( U, 0, ibufS, LDU ), 1, type[I_SEND], partner, Cmsgid, comm );
#else
         MPI_Send_Mod( Mptr( U, 0, ibufS, LDU ), lengthS*LDU, MPI_DOUBLE, partner, Cmsgid, comm );
#endif
	  }

	  //Destroy Data Types
#ifndef HPL_SEND_U_PADDING
      if( lengthR > 0 )
      {
         MPI_Type_free( &type[I_RECV] );
	  }
      if( lengthS > 0 )
      {
         MPI_Type_free( &type[I_SEND] );
	  }
#endif
   }

   if( ierr != MPI_SUCCESS )
   { HPL_pabort( __LINE__, "HPL_rollT", "MPI call failed" ); }

END_TRACE
/*
 * End of HPL_rollT
 */
}
