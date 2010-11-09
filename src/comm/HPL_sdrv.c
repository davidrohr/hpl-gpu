/* 
 * This is a modified version of the High Performance Computing Linpack
 * Benchmark (HPL). All code not contained in the original HPL version
 * 2.0 is property of the Frankfurt Institute for Advanced Studies
 * (FIAS). None of the material may be copied, reproduced, distributed,
 * republished, downloaded, displayed, posted or transmitted in any form
 * or by any means, including, but not limited to, electronic,
 * mechanical, photocopying, recording, or otherwise, without the prior
 * written permission of FIAS. For those parts contained in the
 * unmodified version of the HPL the below copyright notice applies.
 * 
 * Authors:
 * David Rohr (drohr@jwdt.org)
 * Matthias Bach (bach@compeng.uni-frankfurt.de)
 * Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 * 
 * -- High Performance Computing Linpack Benchmark (HPL)                
 *    HPL - 2.0 - September 10, 2008                          
 *    Antoine P. Petitet                                                
 *    University of Tennessee, Knoxville                                
 *    Innovative Computing Laboratory                                 
 *    (C) Copyright 2000-2008 All Rights Reserved                       
 *                                                                      
 * -- Copyright notice and Licensing terms:                             
 *                                                                      
 * Redistribution  and  use in  source and binary forms, with or without
 * modification, are  permitted provided  that the following  conditions
 * are met:                                                             
 *                                                                      
 * 1. Redistributions  of  source  code  must retain the above copyright
 * notice, this list of conditions and the following disclaimer.        
 *                                                                      
 * 2. Redistributions in binary form must reproduce  the above copyright
 * notice, this list of conditions,  and the following disclaimer in the
 * documentation and/or other materials provided with the distribution. 
 *                                                                      
 * 3. All  advertising  materials  mentioning  features  or  use of this
 * software must display the following acknowledgement:                 
 * This  product  includes  software  developed  at  the  University  of
 * Tennessee, Knoxville, Innovative Computing Laboratory.             
 *                                                                      
 * 4. The name of the  University,  the name of the  Laboratory,  or the
 * names  of  its  contributors  may  not  be used to endorse or promote
 * products  derived   from   this  software  without  specific  written
 * permission.                                                          
 *                                                                      
 * -- Disclaimer:                                                       
 *                                                                      
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR  CONTRIBUTORS  BE  LIABLE FOR ANY  DIRECT,  INDIRECT,  INCIDENTAL,
 * SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES  (INCLUDING,  BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA OR PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT,  STRICT LIABILITY,  OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 * ---------------------------------------------------------------------
 */ 
/*
 * Include files
 */
#include "hpl.h"

#include "util_timer.h"
#include "util_trace.h"
/*
 * Do not use  MPI user-defined data types no matter what.  This routine
 * is used for small contiguous messages.
 */
#ifdef HPL_USE_MPI_DATATYPE
#undef HPL_USE_MPI_DATATYPE
#endif

int HPL_sdrv
(
   double *                         SBUF,
   int                              SCOUNT,
   int                              STAG,
   double *                         RBUF,
   int                              RCOUNT,
   int                              RTAG,
   int                              PARTNER,
   MPI_Comm                         COMM
)
{
/* 
 * Purpose
 * =======
 *
 * HPL_sdrv is a simple wrapper around MPI_Sendrecv. Its main purpose is
 * to allow for some experimentation and tuning of this simple function.
 * Messages  of  length  less than  or  equal to zero  are not sent  nor
 * received.  Successful completion  is  indicated by the returned error
 * code HPL_SUCCESS.
 *
 * Arguments
 * =========
 *
 * SBUF    (local input)                 double *
 *         On entry, SBUF specifies the starting address of buffer to be
 *         sent.
 *
 * SCOUNT  (local input)                 int
 *         On entry,  SCOUNT  specifies  the number  of double precision
 *         entries in SBUF. SCOUNT must be at least zero.
 *
 * STAG    (local input)                 int
 *         On entry,  STAG  specifies the message tag to be used for the
 *         sending communication operation.
 *
 * RBUF    (local output)                double *
 *         On entry, RBUF specifies the starting address of buffer to be
 *         received.
 *
 * RCOUNT  (local input)                 int
 *         On entry,  RCOUNT  specifies  the number  of double precision
 *         entries in RBUF. RCOUNT must be at least zero.
 *
 * RTAG    (local input)                 int
 *         On entry,  RTAG  specifies the message tag to be used for the
 *         receiving communication operation.
 *
 * PARTNER (local input)                 int
 *         On entry,  PARTNER  specifies  the rank of the  collaborative
 *         process in the communication space defined by COMM.
 *
 * COMM    (local input)                 MPI_Comm
 *         The MPI communicator identifying the communication space.
 *
 * ---------------------------------------------------------------------
 */ 
START_TRACE( SDRV )

/*
 * .. Local Variables ..
 */
#ifdef HPL_USE_MPI_DATATYPE
   MPI_Datatype               type[2];
#endif
   MPI_Request                request;
   MPI_Status                 status;
   int                        ierr;
/* ..
 * .. Executable Statements ..
 */
   if( RCOUNT > 0 )
   {
      if( SCOUNT > 0 )
      {
#ifdef HPL_USE_MPI_DATATYPE
/*
 * Post asynchronous receive
 */
checkMpiThread          ierr =      MPI_Type_contiguous( RCOUNT, MPI_DOUBLE, &type[0] );
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Type_commit( &type[0] );
         }
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Irecv( (void *)(RBUF), 1, type[0], PARTNER,
                                RTAG, COMM, &request );
         }
/*
 * Blocking send
 */
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Type_contiguous( SCOUNT, MPI_DOUBLE, &type[1] );
         }
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Type_commit( &type[1] );
         }
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Send( (void *)(SBUF), 1, type[1], PARTNER,
                               STAG, COMM );
         }
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Type_free( &type[1] );
         }
/*
 * Wait for the receive to complete
 */
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Wait( &request, &status );
         }
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Type_free( &type[0] );
         }
#else
/*
 * Post asynchronous receive
 */
checkMpiThread          ierr =      MPI_Irecv( (void *)(RBUF), RCOUNT, MPI_DOUBLE,
                                PARTNER, RTAG, COMM, &request );
/*
 * Blocking send
 */
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Send( (void *)(SBUF), SCOUNT, MPI_DOUBLE,
                               PARTNER, STAG, COMM );
         }
/*
 * Wait for the receive to complete
 */
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Wait( &request, &status );
         }
#endif
      }
      else
      {
/*
 * Blocking receive
 */
#ifdef HPL_USE_MPI_DATATYPE
checkMpiThread          ierr =      MPI_Type_contiguous( RCOUNT, MPI_DOUBLE, &type[0] );
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Type_commit( &type[0] );
         }
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Recv( (void *)(RBUF), 1, type[0], PARTNER, RTAG,
                               COMM, &status );
         }
         if( ierr == MPI_SUCCESS ) {
checkMpiThread             ierr =   MPI_Type_free( &type[0] );
         }
#else
checkMpiThread          ierr =      MPI_Recv( (void *)(RBUF), RCOUNT, MPI_DOUBLE,
                               PARTNER, RTAG, COMM, &status );
#endif
      }
   }
   else if( SCOUNT > 0 )
   {
/*
 * Blocking send
 */
#ifdef HPL_USE_MPI_DATATYPE
checkMpiThread       ierr =      MPI_Type_contiguous( SCOUNT, MPI_DOUBLE, &type[1] );
      if( ierr == MPI_SUCCESS ) {
checkMpiThread          ierr =   MPI_Type_commit( &type[1] );
      }
      if( ierr == MPI_SUCCESS ) {
checkMpiThread          ierr =   MPI_Send( (void *)(SBUF), 1, type[1], PARTNER, STAG,
                          COMM );
      }
      if( ierr == MPI_SUCCESS ) {
checkMpiThread          ierr =   MPI_Type_free( &type[1] ) );
      }
#else
checkMpiThread       ierr =      MPI_Send( (void *)(SBUF), SCOUNT, MPI_DOUBLE, PARTNER,
                            STAG, COMM );
#endif
   }
   else { ierr = MPI_SUCCESS; }

END_TRACE
   return( ( ierr == MPI_SUCCESS ? HPL_SUCCESS : HPL_FAILURE ) );
/*
 * End of HPL_sdrv
 */
}
