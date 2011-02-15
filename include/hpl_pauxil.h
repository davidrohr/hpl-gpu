/*
 *  -- High Performance Computing Linpack Benchmark (HPL-GPU)
 *     HPL-GPU - 1.0 - 2010
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

#ifndef HPL_PAUXIL_H
#define HPL_PAUXIL_H
/*
 * ---------------------------------------------------------------------
 * Include files
 * ---------------------------------------------------------------------
 */
#include "hpl_misc.h"
#include "hpl_blas.h"
#include "hpl_auxil.h"

#include "hpl_pmisc.h"
#include "hpl_grid.h"
/*
 * ---------------------------------------------------------------------
 * #define macros definitions
 * ---------------------------------------------------------------------
 */
/*
 * Mindxg2p returns the process coodinate owning the entry globally in-
 * dexed by ig_.
 */
#define Mindxg2p( ig_, inb_, nb_, proc_, nprocs_ ) \
		 { \
		 if( ( (ig_) >= (inb_) ) && ( (nprocs_) > 1 ) ) \
			 { \
			 proc_ = 1 + ( (ig_)-(inb_) ) / (nb_); \
			 proc_ -= ( proc_ / (nprocs_) ) * (nprocs_); \
			 } \
			 else \
			 { \
			 proc_ = 0; \
			 } \
		 }

#define Mindxg2l( il_, ig_, inb_, nb_, proc_, nprocs_ ) \
		 { \
		 if( ( (ig_) < (inb_) ) || ( (nprocs_) == 1 ) ) { il_ = (ig_); } \
			 else \
			 { \
			 int i__, j__; \
			 j__ = ( i__ = ( (ig_)-(inb_) ) / (nb_) ) / (nprocs_); \
			 il_ = (nb_)*( j__ - i__ ) + \
			 ( (i__ + 1 - ( j__ + 1 ) * (nprocs_) ) ? \
			 (ig_) - (inb_) : (ig_) ); \
			 } \
		 }

/*
 * Mindxl2g computes the global index ig_ corresponding to the local
 * index il_ in process proc_.
 */
#define Mindxl2g( ig_, il_, inb_, nb_, proc_, nprocs_ ) \
		 { \
		 if( ( (nprocs_) > 1 ) ) \
			 { \
			 if( (proc_) == 0 ) \
				 { \
				 if( (il_) < (inb_) ) ig_ = (il_); \
					else ig_ = (il_) + \
					(nb_)*((nprocs_)-1)*(((il_)-(inb_))/(nb_) + 1); \
				 } \
				 else \
				 { \
				 ig_ = (il_) + (inb_) + \
				 (nb_)*( ((nprocs_)-1)*((il_)/(nb_)) + \
				 (proc_)--1 ); \
				 } \
			 } \
			 else \
			 { \
			 ig_ = (il_); \
			 } \
		 }
/*
 * MnumrocI computes the # of local indexes np_ residing in the process
 * of coordinate proc_ corresponding to the interval of global indexes
 * i_:i_+n_-1 assuming that the global index 0 resides in the process
 * src_, and that the indexes are distributed from src_ using the para-
 * meters inb_, nb_ and nprocs_.
 */
#define MnumrocI( np_, n_, i_, inb_, nb_, proc_, nprocs_ ) \
		 { \
		 if( ( (nprocs_) > 1 ) ) \
			 { \
			 int inb__, mydist__, n__, nblk__, quot__, src__; \
			 if( ( inb__ = (inb_) - (i_) ) <= 0 ) \
				 { \
				 nblk__ = (-inb__) / (nb_) + 1; \
				 src__ = nblk__; \
				 src__ -= ( src__ / (nprocs_) ) * (nprocs_); \
				 inb__ += nblk__*(nb_); \
				 if( ( n__ = (n_) - inb__ ) <= 0 ) \
					{ \
					if( (proc_) == src__ ) np_ = (n_); \
					 else np_ = 0; \
					} \
					else \
					{ \
					if( ( mydist__ = (proc_) - src__ ) < 0 ) \
					mydist__ += (nprocs_); \
					nblk__ = n__ / (nb_) + 1; \
					mydist__ -= nblk__ - \
					(quot__ = (nblk__ / (nprocs_))) * (nprocs_); \
					if( mydist__ < 0 ) \
					 { \
					 if( (proc_) != src__ ) \
					 np_ = (nb_) + (nb_) * quot__; \
						 else \
						 np_ = inb__ + (nb_) * quot__; \
					 } \
					 else if( mydist__ > 0 ) \
					 { \
					 np_ = (nb_) * quot__; \
					 } \
					 else \
					 { \
					 if( (proc_) != src__ ) \
					 np_ = n__ +(nb_)+(nb_)*(quot__ - nblk__); \
						 else \
						 np_ = (n_)+ (nb_)*(quot__ - nblk__); \
					 } \
					} \
				 } \
				 else \
				 { \
				 if( ( n__ = (n_) - inb__ ) <= 0 ) \
					{ \
					if( (proc_) == 0 ) np_ = (n_); \
					 else np_ = 0; \
					} \
					else \
					{ \
					if( ( mydist__ = (proc_) ) < 0 ) \
					mydist__ += (nprocs_); \
					nblk__ = n__ / (nb_) + 1; \
					mydist__ -= nblk__ - \
					( quot__ = (nblk__ / (nprocs_)) )*(nprocs_); \
					if( mydist__ < 0 ) \
					 { \
					 if( (proc_) != 0 ) \
					 np_ = (nb_) + (nb_) * quot__; \
						 else \
						 np_ = inb__ + (nb_) * quot__; \
					 } \
					 else if( mydist__ > 0 ) \
					 { \
					 np_ = (nb_) * quot__; \
					 } \
					 else \
					 { \
					 if( (proc_) != 0 ) \
					 np_ = n__ +(nb_)+(nb_)*(quot__ - nblk__); \
						 else \
						 np_ = (n_)+ (nb_)*(quot__ - nblk__); \
					 } \
					} \
				 } \
			 } \
			 else \
			 { \
			 np_ = (n_); \
			 } \
		 }

#define Mnumroc( np_, n_, inb_, nb_, proc_, nprocs_ ) \
	MnumrocI( np_, n_, 0, inb_, nb_, proc_, nprocs_ )
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
int HPL_indxg2p( const int, const int, const int, const int );
void HPL_infog2l( int, int, const int, const int, const int, const int, const int, const int, const int, const int, const int, const int, int *, int *, int *, int * );
int HPL_numroc( const int, const int, const int, const int, const int );
int HPL_numrocI( const int, const int, const int, const int, const int, const int );

void HPL_dlaswp00N( const int, const int, double *, const int, const int * );
void HPL_dlaswp10N( const int, const int, double *, const int, const int * );
void HPL_dlaswp01N( const int, const int, double *, const int, double *, const int, const int *, const int * );
void HPL_dlaswp01T( const int, const int, double *, const int, double *, const int, const int *, const int * );
void HPL_dlaswp02N( const int, const int, const double *, const int, double *, double *, const int, const int *, const int * );
void HPL_dlaswp03N( const int, const int, double *, const int, const double *, const double *, const int );
void HPL_dlaswp03T( const int, const int, double *, const int, const double *, const double *, const int );
void HPL_dlaswp04N( const int, const int, const int, double *, const int, double *, const int, const double *, const double *, const int, const int *, const int * );
void HPL_dlaswp04T( const int, const int, const int, double *, const int, double *, const int, const double *, const double *, const int, const int *, const int * );
void HPL_dlaswp05N( const int, const int, double *, const int, const double *, const int, const int *, const int * );
void HPL_dlaswp05T( const int, const int, double *, const int, const double *, const int, const int *, const int * );
void HPL_dlaswp06N( const int, const int, double *, const int, double *, const int, const int * );
void HPL_dlaswp06T( const int, const int, double *, const int, double *, const int, const int * );

void HPL_pabort( int, const char *, const char *, ...);
void HPL_pwarn( FILE *, int, const char *, const char *, ... );
double HPL_pdlamch( MPI_Comm, const HPL_T_MACH );
double HPL_pdlange( const HPL_T_grid *, const HPL_T_NORM, const int, const int, const int, const double *, const int );

#endif
/*
 * End of hpl_pauxil.h
 */
