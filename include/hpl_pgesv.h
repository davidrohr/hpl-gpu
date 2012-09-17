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

#ifndef HPL_PGESV_H
#define HPL_PGESV_H
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
#include "hpl_comm.h"
#include "hpl_pauxil.h"
#include "hpl_panel.h"
#include "hpl_pfact.h"
/*
 * ---------------------------------------------------------------------
 * #typedefs and data structures
 * ---------------------------------------------------------------------
 */

typedef struct HPL_S_palg
{
	HPL_T_TOP btopo; /* row broadcast topology */
	int depth; /* look-ahead depth */
	int nbdiv; /* recursive division factor */
	int nbmin; /* recursion stopping criterium */
	HPL_T_FACT pfact; /* panel fact variant */
	HPL_T_FACT rfact; /* recursive fact variant */
	HPL_T_PFA_FUN pffun; /* panel fact function ptr */
	HPL_T_RFA_FUN rffun; /* recursive fact function ptr */
	int align; /* data alignment constant */
} HPL_T_palg;

typedef struct HPL_S_pmat
{
	double * A; /* pointer to local piece of A */
	double * X; /* pointer to solution vector */
	int n; /* global problem size */
	int nb; /* blocking factor */
	int ld; /* local leading dimension */
	int mp; /* local number of rows */
	int nq; /* local number of columns */
	int info; /* computational flag */
} HPL_T_pmat;
/*
 * ---------------------------------------------------------------------
 * #define macro constants
 * ---------------------------------------------------------------------
 */
#define MSGID_BEGIN_PFACT 1001 /* message id ranges */
#define MSGID_END_PFACT 2000
#define MSGID_BEGIN_FACT 2001
#define MSGID_END_FACT 3000
#define MSGID_BEGIN_PTRSV 3001
#define MSGID_END_PTRSV 4000
 
#define MSGID_BEGIN_COLL 9001
#define MSGID_END_COLL 10000
/*
 * ---------------------------------------------------------------------
 * #define macros definitions
 * ---------------------------------------------------------------------
 */
#define MNxtMgid( id_, beg_, end_ ) \
	(( (id_)+1 > (end_) ? (beg_) : (id_)+1 ))
/*
 * ---------------------------------------------------------------------
 * Function prototypes
 * ---------------------------------------------------------------------
 */
void HPL_pipid( HPL_T_panel *, int *, int * );
void HPL_pdlaswp00N( HPL_T_panel *, int *, HPL_T_panel *, const int );
void HPL_pdlaswp00T( HPL_T_panel *, int *, HPL_T_panel *, const int );

void HPL_perm( const int, int *, int *, int * );
void HPL_logsort( const int, const int, int *, int *, int * );
void HPL_plindx10( HPL_T_panel *, const int, const int *, int *, int *, int * );
void HPL_plindx1( HPL_T_panel *, const int, const int *, int *, int *, int *, int *, int *, int *, int *, int * );
void HPL_spreadT( HPL_T_panel *, const enum HPL_SIDE, const int, double *, const int, const int, const int *, const int *, const int * );
void HPL_equil( HPL_T_panel *, const int, double *, const int, int *, const int *, const int *, int * );
void HPL_rollT( HPL_T_panel *, const int, double *, const int, const int *, const int *, const int * );
int* HPL_pdlaswp01T( HPL_T_panel *, const int );

void HPL_pdgesv( HPL_T_grid *, HPL_T_palg *, HPL_T_pmat * );
 
void HPL_pdtrsv( HPL_T_grid *, HPL_T_pmat * );

#endif
/*
 * End of hpl_pgesv.h
 */ 
