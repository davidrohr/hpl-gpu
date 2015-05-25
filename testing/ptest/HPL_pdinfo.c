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
#include <unistd.h>
#include <stdlib.h>

struct runtime_config_options global_runtime_config;

extern int max_gpu_nb;
extern int max_gpu_nb_factor;

void HPL_pdinfo
(
   HPL_T_test *                     TEST,
   int *                            NS,
   int *                            N,
   int *                            NBS,
   int *                            NB,
   HPL_T_ORDER *                    PMAPPIN,
   int *                            NPQS,
   int *                            P,
   int *                            Q,
   int *                            NPFS,
   HPL_T_FACT *                     PF,
   int *                            NBMS,
   int *                            NBM,
   int *                            NDVS,
   int *                            NDV,
   int *                            NRFS,
   HPL_T_FACT *                     RF,
   int *                            NTPS,
   HPL_T_TOP *                      TP,
   int *                            NDHS,
   int *                            DH,
   int *                            ALIGN,
   int *                            SEED
)
{
/* 
 * Purpose
 * =======
 *
 * HPL_pdinfo reads  the  startup  information for the various tests and
 * transmits it to all processes.
 *
 * Arguments
 * =========
 *
 * TEST    (global output)               HPL_T_test *
 *         On entry, TEST  points to a testing data structure.  On exit,
 *         the fields of this data structure are initialized as follows:
 *         TEST->outfp  specifies the output file where the results will
 *         be printed.  It is only defined and used by  the process 0 of
 *         the grid.  TEST->thrsh specifies the threshhold value for the
 *         test ratio.  TEST->epsil is the relative machine precision of
 *         the distributed computer.  Finally  the test counters, kfail,
 *         kpass, kskip, ktest are initialized to zero.
 *
 * NS      (global output)               int *
 *         On exit,  NS  specifies the number of different problem sizes
 *         to be tested. NS is less than or equal to HPL_MAX_PARAM.
 *
 * N       (global output)               int *
 *         On entry, N is an array of dimension HPL_MAX_PARAM.  On exit,
 *         the first NS entries of this array contain the  problem sizes
 *         to run the code with.
 *
 * NBS     (global output)               int *
 *         On exit,  NBS  specifies the number of different distribution
 *         blocking factors to be tested. NBS must be less than or equal
 *         to HPL_MAX_PARAM.
 *
 * NB      (global output)               int *
 *         On exit,  PMAPPIN  specifies the process mapping onto the no-
 *         des of the  MPI machine configuration.  PMAPPIN  defaults  to
 *         row-major ordering.
 *
 * PMAPPIN (global output)               HPL_T_ORDER *
 *         On entry, NB is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first NBS entries of this array contain the values of the
 *         various distribution blocking factors, to run the code with.
 *
 * NPQS    (global output)               int *
 *         On exit, NPQS  specifies the  number of different values that
 *         can be used for P and Q, i.e., the number of process grids to
 *         run  the  code with.  NPQS must be  less  than  or  equal  to
 *         HPL_MAX_PARAM.
 *
 * P       (global output)               int *
 *         On entry, P  is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first NPQS entries of this array contain the values of P,
 *         the number of process rows of the  NPQS grids to run the code
 *         with.
 *
 * Q       (global output)               int *
 *         On entry, Q  is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first NPQS entries of this array contain the values of Q,
 *         the number of process columns of the  NPQS  grids to  run the
 *         code with.
 *
 * NPFS    (global output)               int *
 *         On exit, NPFS  specifies the  number of different values that
 *         can be used for PF : the panel factorization algorithm to run
 *         the code with. NPFS is less than or equal to HPL_MAX_PARAM.
 *
 * PF      (global output)               HPL_T_FACT *
 *         On entry, PF is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first  NPFS  entries  of this array  contain  the various
 *         panel factorization algorithms to run the code with.
 *
 * NBMS    (global output)               int *
 *         On exit,  NBMS  specifies  the  number  of  various recursive
 *         stopping criteria  to be tested.  NBMS  must be  less than or
 *         equal to HPL_MAX_PARAM.
 *
 * NBM     (global output)               int *
 *         On entry,  NBM  is an array of  dimension  HPL_MAX_PARAM.  On
 *         exit, the first NBMS entries of this array contain the values
 *         of the various recursive stopping criteria to be tested.
 *
 * NDVS    (global output)               int *
 *         On exit,  NDVS  specifies  the number  of various numbers  of
 *         panels in recursion to be tested.  NDVS is less than or equal
 *         to HPL_MAX_PARAM.
 *
 * NDV     (global output)               int *
 *         On entry,  NDV  is an array of  dimension  HPL_MAX_PARAM.  On
 *         exit, the first NDVS entries of this array contain the values
 *         of the various numbers of panels in recursion to be tested.
 *
 * NRFS    (global output)               int *
 *         On exit, NRFS  specifies the  number of different values that
 *         can be used for RF : the recursive factorization algorithm to
 *         be tested. NRFS is less than or equal to HPL_MAX_PARAM.
 *
 * RF      (global output)               HPL_T_FACT *
 *         On entry, RF is an array of dimension HPL_MAX_PARAM. On exit,
 *         the first  NRFS  entries  of  this array contain  the various
 *         recursive factorization algorithms to run the code with.
 *
 * NTPS    (global output)               int *
 *         On exit, NTPS  specifies the  number of different values that
 *         can be used for the  broadcast topologies  to be tested. NTPS
 *         is less than or equal to HPL_MAX_PARAM.
 *
 * TP      (global output)               HPL_T_TOP *
 *         On entry, TP is an array of dimension HPL_MAX_PARAM. On exit,
 *         the  first NTPS  entries of this  array  contain  the various
 *         broadcast (along rows) topologies to run the code with.
 *
 * NDHS    (global output)               int *
 *         On exit, NDHS  specifies the  number of different values that
 *         can be used for the  lookahead depths to be  tested.  NDHS is
 *         less than or equal to HPL_MAX_PARAM.
 *
 * DH      (global output)               int *
 *         On entry,  DH  is  an array of  dimension  HPL_MAX_PARAM.  On
 *         exit, the first NDHS entries of this array contain the values
 *         of lookahead depths to run the code with.  Such a value is at
 *         least 0 (no-lookahead) or greater than zero.
 *
 * ALIGN   (global output)               int *
 *         On exit,  ALIGN  specifies the alignment  of  the dynamically
 *         allocated buffers in double precision words. ALIGN is greater
 *         than zero.
 *
 * SEED    (global output)               int *
 *         On exit,  SEED specifies the seed to be used for matrix
 *         generation. SEED is greater than zero.
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   char                       file[HPL_LINE_MAX], line[HPL_LINE_MAX],
                              auth[HPL_LINE_MAX], num [HPL_LINE_MAX];
   FILE                       * infp;
   int                        * iwork = NULL;
   char                       * lineptr;
   int                        error=0, fid, i, j, lwork, maxp, nprocs,
                              rank, size;

   char                       output_buffer[16384];
/* ..
 * .. Executable Statements ..
 * 
 * 
 */
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   MPI_Comm_size( MPI_COMM_WORLD, &size );
/*
 * Initialize the TEST data structure with default values
 */
   TEST->outfp = stderr; TEST->epsil = 2.0e-16; TEST->thrsh = 16.0;
   TEST->kfail = TEST->kpass = TEST->kskip = TEST->ktest = 0;
   
   //Read node-performance files
   TEST->node_perf = (float*) malloc(size * sizeof(float));
   float node_perf = 1.;
   if ((infp = fopen("node-perf.dat", "r")) != 0)
   {
     char hostname[256];
     gethostname(hostname, 255);
     
     char buffer[256];
     char tmpname[256];
     float tmpperf;
     int tmprank;
     while (!feof(infp))
     {
#pragma GCC diagnostic ignored "-Wunused-result"
		fgets(buffer, 256, infp);
        if (buffer[0] == '#') continue;
        if (buffer[0] == '/')
        {
           sscanf(buffer, "/%d %f", &tmprank, &tmpperf);
           if (tmprank == rank)
           {
              node_perf = tmpperf;
              break;
           }
        }
        else
        {
           sscanf(buffer, "%s %f", tmpname, &tmpperf);
           if (strcmp(tmpname, hostname) == 0 && tmpperf > 0 && tmpperf <= 1.)
           {
              node_perf = tmpperf;
              break;
           }
        }
     }
     fclose(infp);
   }
   MPI_Allgather(&node_perf, 1, MPI_FLOAT, TEST->node_perf, 1, MPI_FLOAT, MPI_COMM_WORLD);
   
   if ( rank == 0)
   {
      for (i = 0;i < size;i++)
      {
         fprintfctd( TEST->outfp, "Performance of rank %d: %f\n", i, TEST->node_perf[i]);
      }
   }
   
/*
 * Process 0 reads the input data, broadcasts to other processes and
 * writes needed information to TEST->outfp.
 */
   if( rank == 0 )
   {
/*
 * Open file and skip data file header
 */
      if( ( infp = fopen( "HPL.dat", "r" ) ) == NULL )
      { 
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                    "cannot open file HPL.dat" );
         error = 1; goto label_error;
      }

      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) fgets( auth, HPL_LINE_MAX - 2, infp );
/*
 * Read name and unit number for summary output file
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", file );
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num  );
      fid = atoi( num );
      if     ( fid == 6 ) TEST->outfp = stdout;
      else if( fid == 7 ) TEST->outfp = stderr;
      else if( ( TEST->outfp = fopen( file, "w" ) ) == NULL )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "cannot open file %s.",
                    file );
         error = 1; goto label_error;
      }
/*
 * Read and check the parameter values for the tests.
 *
 * Problem size (>=0) (N)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); 
      (void) sscanf( line, "%s", num ); *NS = atoi( num );
      if( ( *NS < 1 ) || ( *NS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %d",
                    "Number of values of N is less than 1 or greater than",
                    HPL_MAX_PARAM );
         error = 1; goto label_error;
      }

      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( N[ i ] = atoi( num ) ) < 0 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of N less than 0" );
            error = 1; goto label_error;
         }
      }
/*
 * Block size (>=1) (NB)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NBS = atoi( num );
      if( ( *NBS < 1 ) || ( *NBS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of NB is less than 1 or",
                    "greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }

      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NBS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( NB[ i ] = atoi( num ) ) < 1 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", 
                       "Value of NB less than 1" );
            error = 1; goto label_error;
         }
#ifndef HPL_GPU_MAX_NB
		 if (NB[i] > max_gpu_nb) max_gpu_nb = NB[i];
#else
		 if (NB[i] > HPL_GPU_MAX_NB)
		 {
			HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "Value of NB exceeds HPL_GPU_MAX_NB" );
			error = 1;
			goto label_error;
		 }
#endif
      }
/*
 * Process grids, mapping, (>=1) (P, Q)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num );
      *PMAPPIN = ( atoi( num ) == 1 ? HPL_COLUMN_MAJOR : HPL_ROW_MAJOR );

      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NPQS = atoi( num );
      if( ( *NPQS < 1 ) || ( *NPQS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of grids is less",
                    "than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }

      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NPQS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( P[ i ] = atoi( num ) ) < 1 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of P less than 1" );
            error = 1; goto label_error;
         }
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NPQS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( Q[ i ] = atoi( num ) ) < 1 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of Q less than 1" );
            error = 1; goto label_error;
         }
      }
/*
 * Check for enough processes in machine configuration
 */
      maxp = 0;
      for( i = 0; i < *NPQS; i++ )
      { nprocs   = P[i] * Q[i]; maxp = Mmax( maxp, nprocs ); }
      if( maxp > size )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                    "Need at least %d processes for these tests", maxp );
         error = 1; goto label_error;
      }
/*
 * Checking threshold value (TEST->thrsh)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
#if defined(HPL_START_PERCENTAGE) | defined(HPL_END_N)
      TEST->thrsh = -1;
#else
	  if (global_runtime_config.fastrand == 1)
	  {
		  TEST->thrsh = -1;
	  }
	  else
	  {
         (void) sscanf( line, "%s", num ); TEST->thrsh = atof( num );
	  }
#endif

/*
 * Panel factorization algorithm (PF)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NPFS = atoi( num );
      if( ( *NPFS < 1 ) || ( *NPFS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "number of values of PFACT",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NPFS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         j = atoi( num );
         if(      j == 0 ) PF[ i ] = HPL_LEFT_LOOKING;
         else if( j == 1 ) PF[ i ] = HPL_CROUT;
         else if( j == 2 ) PF[ i ] = HPL_RIGHT_LOOKING;
         else              PF[ i ] = HPL_RIGHT_LOOKING;
      }
/*
 * Recursive stopping criterium (>=1) (NBM)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NBMS = atoi( num );
      if( ( *NBMS < 1 ) || ( *NBMS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of NBMIN",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NBMS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( NBM[ i ] = atoi( num ) ) < 1 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of NBMIN less than 1" );
            error = 1; goto label_error;
         }
      }
/*
 * Number of panels in recursion (>=2) (NDV)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NDVS = atoi( num );
      if( ( *NDVS < 1 ) || ( *NDVS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of NDIV",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NDVS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         if( ( NDV[ i ] = atoi( num ) ) < 2 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of NDIV less than 2" );
            error = 1; goto label_error;
         }
      }
/*
 * Recursive panel factorization (RF)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NRFS = atoi( num );
      if( ( *NRFS < 1 ) || ( *NRFS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of RFACT",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NRFS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         j = atoi( num );
         if(      j == 0 ) RF[ i ] = HPL_LEFT_LOOKING;
         else if( j == 1 ) RF[ i ] = HPL_CROUT;
         else if( j == 2 ) RF[ i ] = HPL_RIGHT_LOOKING;
         else              RF[ i ] = HPL_RIGHT_LOOKING;
      }
/*
 * Broadcast topology (TP) (0=rg, 1=2rg, 2=rgM, 3=2rgM, 4=L)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NTPS = atoi( num );
      if( ( *NTPS < 1 ) || ( *NTPS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of BCAST",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NTPS; i++ )
      {
         (void) sscanf( lineptr, "%s", num ); lineptr += strlen( num ) + 1;
         j = atoi( num );
         if(      j == 0 ) TP[ i ] = HPL_1RING;
         else if( j == 1 ) TP[ i ] = HPL_1RING_M;
         else if( j == 2 ) TP[ i ] = HPL_2RING;
         else if( j == 3 ) TP[ i ] = HPL_2RING_M;
         else if( j == 4 ) TP[ i ] = HPL_BLONG;
         else if( j == 5 ) TP[ i ] = HPL_BLONG_M;
         else if( j == 6 ) TP[ i ] = HPL_MPI_BCAST;
         else              TP[ i ] = HPL_1RING_M;
      }
/*
 * Lookahead depth (>=0) (NDH)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *NDHS = atoi( num );
      if( ( *NDHS < 1 ) || ( *NDHS > HPL_MAX_PARAM ) )
      {
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo", "%s %s %d",
                    "Number of values of DEPTH",
                    "is less than 1 or greater than", HPL_MAX_PARAM );
         error = 1; goto label_error;
      }
      (void) fgets( line, HPL_LINE_MAX - 2, infp ); lineptr = line;
      for( i = 0; i < *NDHS; i++ )
      {
         (void) sscanf( lineptr, "%s", num );
         lineptr += strlen( num ) + 1;
         if( ( DH[ i ] = atoi( num ) ) < 0  || DH[i] > 3 )
         {
            HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                       "Value of DEPTH less than 0 or greater than 3" );
            error = 1; goto label_error;
         }
      }
/*
 * Memory alignment in bytes (> 0) (ALIGN)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *ALIGN = atoi( num );
      if( *ALIGN <= 0 ) *ALIGN = 4;
/*
 * Matrix seed (> 0) (SEED)
 */
      (void) fgets( line, HPL_LINE_MAX - 2, infp );
      (void) sscanf( line, "%s", num ); *SEED = atoi( num );
      if( *SEED <= 0 ) *SEED = HPL_IDEFSEED;
/*
 * Close input file
 */
label_error:
      (void) fclose( infp );
   }
   else { TEST->outfp = NULL; }
/*
 * Check for error on reading input file
 */
   (void) HPL_all_reduce( (void *)(&error), 1, HPL_INT, HPL_max,
                          MPI_COMM_WORLD );
   if( error )
   {
      if( rank == 0 )
         HPL_pwarn( stderr, __LINE__, "HPL_pdinfo",
                    "Illegal input in file HPL.dat. Exiting ..." );
      MPI_Finalize();
      exit( 1 );
   }
/*
 * Compute and broadcast machine epsilon
 */
   TEST->epsil = HPL_pdlamch( MPI_COMM_WORLD, HPL_MACH_EPS );
/*
 * Pack information arrays and broadcast
 */
   (void) HPL_broadcast( (void *)(&(TEST->thrsh)), 1, HPL_DOUBLE, 0,
                         MPI_COMM_WORLD );
/*
 * Broadcast array sizes
 */
   iwork = (int *)malloc( (size_t)(12) * sizeof( int ) );
   if( rank == 0 )
   {
      iwork[ 0] = *NS;      iwork[ 1] = *NBS;
      iwork[ 2] = ( *PMAPPIN == HPL_ROW_MAJOR ? 0 : 1 );
      iwork[ 3] = *NPQS;    iwork[ 4] = *NPFS;     iwork[ 5] = *NBMS;
      iwork[ 6] = *NDVS;    iwork[ 7] = *NRFS;     iwork[ 8] = *NTPS;
      iwork[ 9] = *NDHS;    iwork[10] = *ALIGN;    iwork[11] = *SEED;
   }
   (void) HPL_broadcast( (void *)iwork, 12, HPL_INT, 0, MPI_COMM_WORLD );
   if( rank != 0 )
   {
      *NS       = iwork[ 0]; *NBS   = iwork[ 1];
      *PMAPPIN  = ( iwork[ 2] == 0 ?  HPL_ROW_MAJOR : HPL_COLUMN_MAJOR );
      *NPQS     = iwork[ 3]; *NPFS  = iwork[ 4]; *NBMS     = iwork[ 5];
      *NDVS     = iwork[ 6]; *NRFS  = iwork[ 7]; *NTPS     = iwork[ 8];
      *NDHS     = iwork[ 9]; *ALIGN = iwork[10]; *SEED     = iwork[11];
   }
   if( iwork ) free( iwork );
/*
 * Pack information arrays and broadcast
 */
   lwork = (*NS) + (*NBS) + 2 * (*NPQS) + (*NPFS) + (*NBMS) + 
           (*NDVS) + (*NRFS) + (*NTPS) + (*NDHS);
   iwork = (int *)malloc( (size_t)(lwork) * sizeof( int ) );
   if( rank == 0 )
   {
      j = 0;
      for( i = 0; i < *NS;   i++ ) { iwork[j] = N [i]; j++; }
      for( i = 0; i < *NBS;  i++ ) { iwork[j] = NB[i]; j++; }
      for( i = 0; i < *NPQS; i++ ) { iwork[j] = P [i]; j++; }
      for( i = 0; i < *NPQS; i++ ) { iwork[j] = Q [i]; j++; }
      for( i = 0; i < *NPFS; i++ )
      {
         if(      PF[i] == HPL_LEFT_LOOKING  ) iwork[j] = 0;
         else if( PF[i] == HPL_CROUT         ) iwork[j] = 1;
         else if( PF[i] == HPL_RIGHT_LOOKING ) iwork[j] = 2;
         j++;
      }
      for( i = 0; i < *NBMS; i++ ) { iwork[j] = NBM[i]; j++; }
      for( i = 0; i < *NDVS; i++ ) { iwork[j] = NDV[i]; j++; }
      for( i = 0; i < *NRFS; i++ )
      {
         if(      RF[i] == HPL_LEFT_LOOKING  ) iwork[j] = 0;
         else if( RF[i] == HPL_CROUT         ) iwork[j] = 1;
         else if( RF[i] == HPL_RIGHT_LOOKING ) iwork[j] = 2;
         j++;
      }
      for( i = 0; i < *NTPS; i++ )
      {
         if(      TP[i] == HPL_1RING     ) iwork[j] = 0;
         else if( TP[i] == HPL_1RING_M   ) iwork[j] = 1;
         else if( TP[i] == HPL_2RING     ) iwork[j] = 2;
         else if( TP[i] == HPL_2RING_M   ) iwork[j] = 3;
         else if( TP[i] == HPL_BLONG     ) iwork[j] = 4;
         else if( TP[i] == HPL_BLONG_M   ) iwork[j] = 5;
         else if( TP[i] == HPL_MPI_BCAST ) iwork[j] = 6;
         j++;
      }
      for( i = 0; i < *NDHS; i++ ) { iwork[j] = DH[i]; j++; }
   }
   (void) HPL_broadcast( (void*)iwork, lwork, HPL_INT, 0,
                         MPI_COMM_WORLD );
   if( rank != 0 )
   {
      j = 0;
      for( i = 0; i < *NS;   i++ ) { N [i] = iwork[j]; j++; }
      for( i = 0; i < *NBS;  i++ ) { NB[i] = iwork[j]; j++; }
      for( i = 0; i < *NPQS; i++ ) { P [i] = iwork[j]; j++; }
      for( i = 0; i < *NPQS; i++ ) { Q [i] = iwork[j]; j++; }

      for( i = 0; i < *NPFS; i++ )
      {
         if(      iwork[j] == 0 ) PF[i] = HPL_LEFT_LOOKING;
         else if( iwork[j] == 1 ) PF[i] = HPL_CROUT;
         else if( iwork[j] == 2 ) PF[i] = HPL_RIGHT_LOOKING;
         j++;
      }
      for( i = 0; i < *NBMS; i++ ) { NBM[i] = iwork[j]; j++; }
      for( i = 0; i < *NDVS; i++ ) { NDV[i] = iwork[j]; j++; }
      for( i = 0; i < *NRFS; i++ )
      {
         if(      iwork[j] == 0 ) RF[i] = HPL_LEFT_LOOKING;
         else if( iwork[j] == 1 ) RF[i] = HPL_CROUT;
         else if( iwork[j] == 2 ) RF[i] = HPL_RIGHT_LOOKING;
         j++;
      }
      for( i = 0; i < *NTPS; i++ )
      {
         if(      iwork[j] == 0 ) TP[i] = HPL_1RING;
         else if( iwork[j] == 1 ) TP[i] = HPL_1RING_M;
         else if( iwork[j] == 2 ) TP[i] = HPL_2RING;
         else if( iwork[j] == 3 ) TP[i] = HPL_2RING_M;
         else if( iwork[j] == 4 ) TP[i] = HPL_BLONG;
         else if( iwork[j] == 5 ) TP[i] = HPL_BLONG_M;
         else if( iwork[j] == 6 ) TP[i] = HPL_MPI_BCAST;
         j++;
      }
      for( i = 0; i < *NDHS; i++ ) { DH[i] = iwork[j]; j++; }
   }
   if( iwork ) free( iwork );

    //Broadcast maximum Nb value that will appear during HPL runs (if not defined via HPL_GPU_MAX_NB)
#ifndef HPL_GPU_MAX_NB
   HPL_broadcast( (void*) &max_gpu_nb, 1, HPL_INT, 0, MPI_COMM_WORLD );
#endif


/*
 * regurgitate input
 */
   if( rank == 0 )
   {
      HPL_fprintf( TEST->outfp, "%s%s\n",
                   "========================================",
                   "========================================" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "HPL-GPU " HPL_VERSION "  --  High-Performance Linpack benchmark for GPUs  --  ",
          " 2015" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "Written by D. Rohr, M. Kretz and M. Bach, ",
          "Frankfurt Institute for Advanced Studies" );
      HPL_fprintf( TEST->outfp, "%s\n",
          "This program is open source software licensed partly under GPL and partly\nunder BSD license." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "You can obtain the source at: http://code.compeng.uni-frankfurt.de/projects/hpl." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "Based on:" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "HPLinpack 2.1  --  High-Performance Linpack benchmark  --  ",
          " October 26, 2012" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "Written by A. Petitet and R. Clint Whaley,  ",
          "Innovative Computing Laboratory, UTK" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "Modified by Piotr Luszczek, ",
          "Innovative Computing Laboratory, UTK" );
      HPL_fprintf( TEST->outfp, "%s%s\n",
          "Modified by Julien Langou, ",
          "University of Colorado Denver");
      HPL_fprintf( TEST->outfp, "%s%s\n",
                   "========================================",
                   "========================================" );

      HPL_fprintf( TEST->outfp, "\n%s\n",
          "An explanation of the input/output parameters follows:" );
      HPL_fprintf( TEST->outfp, "%s\n",
          "T/V    : Wall time / encoded variant." );
      HPL_fprintf( TEST->outfp, "%s\n",
         "N      : The order of the coefficient matrix A." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "NB     : The partitioning blocking factor." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "P      : The number of process rows." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "Q      : The number of process columns." );
      HPL_fprintf( TEST->outfp, "%s\n",
         "Time   : Time in seconds to solve the linear system." );
      HPL_fprintf( TEST->outfp, "%s\n\n",
         "Gflops : Rate of execution for solving the linear system." );
      HPL_fprintf( TEST->outfp, "%s\n",
          "The following parameter values will be used:" );
/*
 * Problem size
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nN      :" );
      for( i = 0; i < Mmin( 8, *NS ); i++ )
         sprintf( output_buffer + strlen(output_buffer),       "%8d ", N[i]  );
      if( *NS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NS ); i++ )
            sprintf( output_buffer + strlen(output_buffer),    "%8d ", N[i]  );
         if( *NS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NS; i++ )
               sprintf( output_buffer + strlen(output_buffer), "%8d ", N[i]  );
         }
      }
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Distribution blocking factor
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nNB     :" );
      for( i = 0; i < Mmin( 8, *NBS ); i++ )
         sprintf( output_buffer + strlen(output_buffer),       "%8d ", NB[i] );
      if( *NBS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NBS ); i++ )
            sprintf( output_buffer + strlen(output_buffer),    "%8d ", NB[i] );
         if( *NBS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NBS; i++ )
               sprintf( output_buffer + strlen(output_buffer), "%8d ", NB[i] );
         }
      }
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Process mapping
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nPMAP   :" );
      if(      *PMAPPIN == HPL_ROW_MAJOR    )
         sprintf( output_buffer + strlen(output_buffer), " Row-major process mapping" );
      else if( *PMAPPIN == HPL_COLUMN_MAJOR )
         sprintf( output_buffer + strlen(output_buffer), " Column-major process mapping" );
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Process grid
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nP      :" );
      for( i = 0; i < Mmin( 8, *NPQS ); i++ )
         sprintf( output_buffer + strlen(output_buffer),       "%8d ", P[i]  );
      if( *NPQS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NPQS ); i++ )
            sprintf( output_buffer + strlen(output_buffer),    "%8d ", P[i]  );
         if( *NPQS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NPQS; i++ )
               sprintf( output_buffer + strlen(output_buffer), "%8d ", P[i]  );
         }
      }
      sprintf( output_buffer + strlen(output_buffer),       "\nQ      :" );
      for( i = 0; i < Mmin( 8, *NPQS ); i++ )
         sprintf( output_buffer + strlen(output_buffer),       "%8d ", Q[i]  );
      if( *NPQS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NPQS ); i++ )
            sprintf( output_buffer + strlen(output_buffer),    "%8d ", Q[i]  );
         if( *NPQS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NPQS; i++ )
               sprintf( output_buffer + strlen(output_buffer), "%8d ", Q[i]  );
         }
      }
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Panel Factorization
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nPFACT  :" );
      for( i = 0; i < Mmin( 8, *NPFS ); i++ )
      {
         if(      PF[i] == HPL_LEFT_LOOKING  )
            sprintf( output_buffer + strlen(output_buffer),       "    Left " );
         else if( PF[i] == HPL_CROUT         )
            sprintf( output_buffer + strlen(output_buffer),       "   Crout " );
         else if( PF[i] == HPL_RIGHT_LOOKING )
            sprintf( output_buffer + strlen(output_buffer),       "   Right " );
      }
      if( *NPFS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NPFS ); i++ )
         {
            if(      PF[i] == HPL_LEFT_LOOKING  )
               sprintf( output_buffer + strlen(output_buffer),       "    Left " );
            else if( PF[i] == HPL_CROUT         )
               sprintf( output_buffer + strlen(output_buffer),       "   Crout " );
            else if( PF[i] == HPL_RIGHT_LOOKING )
               sprintf( output_buffer + strlen(output_buffer),       "   Right " );
         }
         if( *NPFS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NPFS; i++ )
            {
               if(      PF[i] == HPL_LEFT_LOOKING  )
                  sprintf( output_buffer + strlen(output_buffer),       "    Left " );
               else if( PF[i] == HPL_CROUT         )
                  sprintf( output_buffer + strlen(output_buffer),       "   Crout " );
               else if( PF[i] == HPL_RIGHT_LOOKING )
                  sprintf( output_buffer + strlen(output_buffer),       "   Right " );
            }
         }
      }
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Recursive stopping criterium
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nNBMIN  :" );
      for( i = 0; i < Mmin( 8, *NBMS ); i++ )
         sprintf( output_buffer + strlen(output_buffer),       "%8d ", NBM[i]  );
      if( *NBMS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NBMS ); i++ )
            sprintf( output_buffer + strlen(output_buffer),    "%8d ", NBM[i]  );
         if( *NBMS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NBMS; i++ )
               sprintf( output_buffer + strlen(output_buffer), "%8d ", NBM[i]  );
         }
      }
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Number of panels in recursion
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nNDIV   :" );
      for( i = 0; i < Mmin( 8, *NDVS ); i++ )
         sprintf( output_buffer + strlen(output_buffer),       "%8d ", NDV[i]  );
      if( *NDVS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NDVS ); i++ )
            sprintf( output_buffer + strlen(output_buffer),    "%8d ", NDV[i]  );
         if( *NDVS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NDVS; i++ )
               sprintf( output_buffer + strlen(output_buffer), "%8d ", NDV[i]  );
         }
      }
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Recursive Factorization
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nRFACT  :" );
      for( i = 0; i < Mmin( 8, *NRFS ); i++ )
      {
         if(      RF[i] == HPL_LEFT_LOOKING  )
            sprintf( output_buffer + strlen(output_buffer),       "    Left " );
         else if( RF[i] == HPL_CROUT         )
            sprintf( output_buffer + strlen(output_buffer),       "   Crout " );
         else if( RF[i] == HPL_RIGHT_LOOKING )
            sprintf( output_buffer + strlen(output_buffer),       "   Right " );
      }
      if( *NRFS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NRFS ); i++ )
         {
            if(      RF[i] == HPL_LEFT_LOOKING  )
               sprintf( output_buffer + strlen(output_buffer),       "    Left " );
            else if( RF[i] == HPL_CROUT         )
               sprintf( output_buffer + strlen(output_buffer),       "   Crout " );
            else if( RF[i] == HPL_RIGHT_LOOKING )
               sprintf( output_buffer + strlen(output_buffer),       "   Right " );
         }
         if( *NRFS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NRFS; i++ )
            {
               if(      RF[i] == HPL_LEFT_LOOKING  )
                  sprintf( output_buffer + strlen(output_buffer),       "    Left " );
               else if( RF[i] == HPL_CROUT         )
                  sprintf( output_buffer + strlen(output_buffer),       "   Crout " );
               else if( RF[i] == HPL_RIGHT_LOOKING )
                  sprintf( output_buffer + strlen(output_buffer),       "   Right " );
            }
         }
      }
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Broadcast topology
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nBCAST  :" );
      for( i = 0; i < Mmin( 8, *NTPS ); i++ )
      {
         if(      TP[i] == HPL_1RING   )
            sprintf( output_buffer + strlen(output_buffer),       "   1ring " );
         else if( TP[i] == HPL_1RING_M )
            sprintf( output_buffer + strlen(output_buffer),       "  1ringM " );
         else if( TP[i] == HPL_2RING   )
            sprintf( output_buffer + strlen(output_buffer),       "   2ring " );
         else if( TP[i] == HPL_2RING_M )
            sprintf( output_buffer + strlen(output_buffer),       "  2ringM " );
         else if( TP[i] == HPL_BLONG   )
            sprintf( output_buffer + strlen(output_buffer),       "   Blong " );
         else if( TP[i] == HPL_BLONG_M )
            sprintf( output_buffer + strlen(output_buffer),       "  BlongM " );
         else if( TP[i] == HPL_MPI_BCAST )
            sprintf( output_buffer + strlen(output_buffer),       "     MPI " );
      }
      if( *NTPS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NTPS ); i++ )
         {
            if(      TP[i] == HPL_1RING   )
               sprintf( output_buffer + strlen(output_buffer),       "   1ring " );
            else if( TP[i] == HPL_1RING_M )
               sprintf( output_buffer + strlen(output_buffer),       "  1ringM " );
            else if( TP[i] == HPL_2RING   )
               sprintf( output_buffer + strlen(output_buffer),       "   2ring " );
            else if( TP[i] == HPL_2RING_M )
               sprintf( output_buffer + strlen(output_buffer),       "  2ringM " );
            else if( TP[i] == HPL_BLONG   )
               sprintf( output_buffer + strlen(output_buffer),       "   Blong " );
            else if( TP[i] == HPL_BLONG_M )
               sprintf( output_buffer + strlen(output_buffer),       "  BlongM " );
            else if( TP[i] == HPL_MPI_BCAST )
               sprintf( output_buffer + strlen(output_buffer),       "     MPI " );
         }
         if( *NTPS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NTPS; i++ )
            {
               if(      TP[i] == HPL_1RING   )
                  sprintf( output_buffer + strlen(output_buffer),       "   1ring " );
               else if( TP[i] == HPL_1RING_M )
                  sprintf( output_buffer + strlen(output_buffer),       "  1ringM " );
               else if( TP[i] == HPL_2RING   )
                  sprintf( output_buffer + strlen(output_buffer),       "   2ring " );
               else if( TP[i] == HPL_2RING_M )
                  sprintf( output_buffer + strlen(output_buffer),       "  2ringM " );
               else if( TP[i] == HPL_BLONG   )
                  sprintf( output_buffer + strlen(output_buffer),       "   Blong " );
               else if( TP[i] == HPL_BLONG_M )
                  sprintf( output_buffer + strlen(output_buffer),       "  BlongM " );
               else if( TP[i] == HPL_MPI_BCAST )
                  sprintf( output_buffer + strlen(output_buffer),       "     MPI " );
            }
         }
      }
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Lookahead depths
 */
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nDEPTH  :" );
      for( i = 0; i < Mmin( 8, *NDHS ); i++ )
         sprintf( output_buffer + strlen(output_buffer),       "%8d ", DH[i]  );
      if( *NDHS > 8 )
      {
         sprintf( output_buffer + strlen(output_buffer),    "\n        " );
         for( i = 8; i < Mmin( 16, *NDHS ); i++ )
            sprintf( output_buffer + strlen(output_buffer),    "%8d ", DH[i]  );
         if( *NDHS > 16 )
         {
            sprintf( output_buffer + strlen(output_buffer), "\n        " );
            for( i = 16; i < *NDHS; i++ )
               sprintf( output_buffer + strlen(output_buffer), "%8d ", DH[i]  );
         }
      }
      HPL_fprintf( TEST->outfp, "%s", output_buffer );
/*
 * Swapping algorithm
 */
      HPL_fprintf( TEST->outfp,       "\nSWAP   : Spread-roll (long)" );
/*
 * L1 storage form
 */
      HPL_fprintf( TEST->outfp,       "\nL1     : transposed form" );
/*
 * U  storage form
 */
      HPL_fprintf( TEST->outfp,       "\nU      : transposed form" );
/*
 * Equilibration
 */
#ifndef NO_EQUILIBRATION
      HPL_fprintf( TEST->outfp,       "\nEQUIL  : yes" );
#else
      HPL_fprintf( TEST->outfp,       "\nEQUIL  : no" );
#endif
/*
 * Alignment
 */
      HPL_fprintf( TEST->outfp,       "\nALIGN  : %d double precision words",
                   *ALIGN );
/*
 * Seed
 */
      HPL_fprintf( TEST->outfp,       "\nSEED   :%8d",
                   *SEED );

/*
 * Config Options
 */
#ifdef HPL_PRINT_CONFIG
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nConfig : ");
#include "hpl_config_option_list.h"
      fprintf( TEST->outfp, "%s", output_buffer );
#define MXSTR(x) #x
	if (strcmp(MXSTR(HPL_CALDGEMM_BACKEND), "cal") == 0)
	{
      output_buffer[0] = 0;
      sprintf( output_buffer + strlen(output_buffer),       "\nCALDGEM: ");
#include "caldgemm_config_option_list.h"
      fprintf( TEST->outfp, "%s", output_buffer );
        }
#endif

      HPL_fprintf( TEST->outfp, "\n\n" );
/*
 * For testing only
 */
      if( TEST->thrsh > HPL_rzero )
      {
         HPL_fprintf( TEST->outfp, "%s%s\n\n",
                      "----------------------------------------",
                      "----------------------------------------" );
         HPL_fprintf( TEST->outfp, "%s\n",
            "- The matrix A is randomly generated for each test." );
         HPL_fprintf( TEST->outfp, "%s\n",
            "- The following scaled residual check will be computed:" );
         HPL_fprintf( TEST->outfp, "%s\n",
            "      ||Ax-b||_oo / ( eps * ( || x ||_oo * || A ||_oo + || b ||_oo ) * N )" );
         HPL_fprintf( TEST->outfp, "%s %21.6e\n",
            "- The relative machine precision (eps) is taken to be     ",
            TEST->epsil );
         HPL_fprintf( TEST->outfp, "%s   %11.1f\n\n",
            "- Computational tests pass if scaled residuals are less than      ",
            TEST->thrsh );
      }
   }


    global_runtime_config.paramdefs = (char*) malloc(1);
    global_runtime_config.paramdefs[0] = 0;
#ifdef HPL_FASTINIT
#ifdef HPL_FASTVERIFY
	global_runtime_config.fastrand = 2;
#else
	global_runtime_config.fastrand = 1;
#endif
#else
	global_runtime_config.fastrand = 0;
#endif
#ifdef HPL_WARMUP
	global_runtime_config.warmup = 1;
#else
	global_runtime_config.warmup = 0;
#endif
#ifdef HPL_DISABLE_LOOKAHEAD
    global_runtime_config.disable_lookahead = HPL_DISABLE_LOOKAHEAD;
#else
    global_runtime_config.disable_lookahead = 0;
#endif
#ifdef HPL_LOOKAHEAD2_TURNOFF
    global_runtime_config.lookahead2_turnoff = HPL_LOOKAHEAD2_TURNOFF;
#else
    global_runtime_config.lookahead2_turnoff = 0;
#endif
#ifdef HPL_DURATION_FIND_HELPER 
    global_runtime_config.duration_find_helper = 1;
#else
    global_runtime_config.duration_find_helper = 0;
#endif
#ifdef HPL_CALDGEMM_ASYNC_FACT_DGEMM
    global_runtime_config.caldgemm_async_fact_dgemm = HPL_CALDGEMM_ASYNC_FACT_DGEMM;
#else
    global_runtime_config.caldgemm_async_fact_dgemm = 0;
#endif
#ifdef HPL_CALDGEMM_ASYNC_FACT_FIRST
    global_runtime_config.caldgemm_async_fact_first = HPL_CALDGEMM_ASYNC_FACT_FIRST ;
#else
    global_runtime_config.caldgemm_async_fact_first = 0;
#endif
#ifdef HPL_CALDGEMM_ASYNC_DTRSM
    global_runtime_config.caldgemm_async_dtrsm = HPL_CALDGEMM_ASYNC_DTRSM;
#else
    global_runtime_config.caldgemm_async_dtrsm = 0;
#endif
#ifdef HPL_CALDGEMM_ASYNC_FACT_DTRSM
    global_runtime_config.caldgemm_async_fact_dtrsm = HPL_CALDGEMM_ASYNC_FACT_DTRSM;
#else
    global_runtime_config.caldgemm_async_fact_dtrsm = 0;
#endif
    global_runtime_config.hpl_nb_multiplier_count = 0;
    for (i = 0;i < HPL_NB_MULTIPLIER_MAX;i++) global_runtime_config.hpl_nb_multiplier_factor[i] = 1;
                                        
    //Read HPL_GPU_CONFIG runtime config file
#ifdef HPL_GPU_RUNTIME_CONFIG
  char* Buffer;
  if ( rank == 0)
  {
    FILE* fRuntime = fopen("HPL-GPU.conf", "r");
    if (fRuntime == 0)
    {
      HPL_fprintf( TEST->outfp, "Error Opening Runtime Config File HPL-GPU.conf!\n");
      exit(1);
    }
    fseek(fRuntime, 0, SEEK_END);
    int filesize = ftell(fRuntime);
    fseek(fRuntime, 0, SEEK_SET);
    Buffer = (char*) malloc(filesize + 1);
    int nread = fread(Buffer, 1, filesize, fRuntime);
    fclose(fRuntime);
    if (nread != filesize)
    {
      HPL_fprintf( TEST->outfp, "Error readubg File HPL-GPU.conf!\n");
      exit(1);
    }
    Buffer[filesize] = 0;
    HPL_broadcast( (void*) &filesize, 1, HPL_INT, 0, MPI_COMM_WORLD );
    HPL_broadcast( (void*) Buffer, filesize + 1, MPI_BYTE, 0, MPI_COMM_WORLD );
  }
  else
  {
    int filesize;
    HPL_broadcast( (void*) &filesize, 1, HPL_INT, 0, MPI_COMM_WORLD );
    Buffer = (char*) malloc(filesize + 1);
    HPL_broadcast( (void*) Buffer, filesize + 1, MPI_BYTE, 0, MPI_COMM_WORLD );
  }
  
  char* ptr = Buffer;
  while (*ptr)
  {
    char *cmd, *option = (char*) "";
    while (*ptr == ' ' || *ptr == '	') ptr++;
    if (*ptr == '#')
    {
	while (*ptr != 10 && *ptr != 13 && *ptr) ptr++;
	while (*ptr == 10 || *ptr == 13) ptr++;
	continue;
    }
    cmd = ptr;
    while (*ptr != ' ' && *ptr != '	' && *ptr != ':' && *ptr != 10 && *ptr != 13 && *ptr) ptr++;
    char* ptr2 = ptr;
    while (*ptr2 == ' ' || *ptr2 == '	') ptr2++;
    if (*ptr2 == ':')
    {
	*ptr = 0;
	ptr = ptr2 + 1;
	while (*ptr == ' ' || *ptr == '	') ptr++;
	option = ptr2 = ptr;
	while (*ptr2 != 10 && *ptr2 != 13 && *ptr2)
	{
	    if (*ptr2 != ' ' && *ptr2 != '	') ptr = ptr2;
	    ptr2++;
	}
	if (*ptr && *ptr != 10 && *ptr != 13) ptr++;
    }
    else if (*ptr != 10 && *ptr != 13)
    {
	HPL_fprintf( TEST->outfp, "Error parsing runtime config file\n");
	exit(1);
    }
    while (*ptr2 == 10 || *ptr2 == 13) ptr2++;
    *ptr = 0;
    ptr = ptr2;

	if (rank == 0 && strcmp(cmd, "HPL_PARAMDEFS") != 0) HPL_fprintf( TEST->outfp, "Runtime Option \"%s\", Parameter \"%s\"\n", cmd, option);
	if (strcmp(cmd, "HPL_WARMUP") == 0)
	{
		global_runtime_config.warmup = option[0] ? atoi(option) : 1;
	}
	else if (strcmp(cmd, "HPL_FASTRAND") == 0)
	{
		global_runtime_config.fastrand = option[0] ? atoi(option) : 2;
	}
	else if (strcmp(cmd, "HPL_DISABLE_LOOKAHEAD") == 0)
	{
		global_runtime_config.disable_lookahead = atoi(option);
	}
	else if (strcmp(cmd, "HPL_LOOKAHEAD2_TURNOFF") == 0)
	{
		global_runtime_config.lookahead2_turnoff = atoi(option);
	}
	else if (strcmp(cmd, "HPL_DURATION_FIND_HELPER") == 0)
	{
		global_runtime_config.duration_find_helper = option[0] ? atoi(option) : 1;
	}
	else if (strcmp(cmd, "HPL_CALDGEMM_ASYNC_FACT_DGEMM") == 0)
	{
		global_runtime_config.caldgemm_async_fact_dgemm = atoi(option);
	}
	else if (strcmp(cmd, "HPL_CALDGEMM_ASYNC_FACT_FIRST") == 0)
	{
		global_runtime_config.caldgemm_async_fact_first = option[0] ? atoi(option) : 1;
	}
	else if (strcmp(cmd, "HPL_CALDGEMM_ASYNC_DTRSM") == 0)
	{
		global_runtime_config.caldgemm_async_dtrsm = atoi(option);
	}
	else if (strcmp(cmd, "HPL_CALDGEMM_ASYNC_FACT_DTRSM") == 0)
	{
		global_runtime_config.caldgemm_async_fact_dtrsm = atoi(option);
	}
	else if (strcmp(cmd, "HPL_NB_MULTIPLIER") == 0)
	{
		char* nbptr = option;
		j = 0;
		int a = strlen(nbptr);
		int nbnum = 0;
		for (i = 0;i <= a;i++)
		{
			if (nbptr[i] == ',' || nbptr[i] == ';' || nbptr[i] == 0)
			{
				if (i > j)
				{
					int tmpval;
					nbptr[i] = 0;
					sscanf(&nbptr[j], "%d", &tmpval);
					j = i + 1;
					if (nbnum >= (signed) HPL_NB_MULTIPLIER_MAX)
					{
						fprintf(STD_OUT, "Please increase HPL_NB_MULTIPLIER_MAX\n");
						break;
					}
					global_runtime_config.hpl_nb_multiplier_factor[nbnum] = tmpval;
					nbnum++;
				}
			}
		}
		global_runtime_config.hpl_nb_multiplier_count = nbnum;
	}
	else if (strcmp(cmd, "HPL_NB_MULTIPLIER_THRESHOLD") == 0)
	{
		char* nbptr = option;
		j = 0;
		int a = strlen(nbptr);
		int nbnum = 0;
		for (i = 0;i <= a;i++)
		{
			if (nbptr[i] == ',' || nbptr[i] == ';' || nbptr[i] == 0)
			{
				if (i > j)
				{
					int tmpval;
					nbptr[i] = 0;
					sscanf(&nbptr[j], "%d", &tmpval);
					j = i + 1;
					if (nbnum >= (signed) HPL_NB_MULTIPLIER_MAX)
					{
						fprintf(STD_OUT, "Please increase HPL_NB_MULTIPLIER_MAX\n");
						break;
					}
					global_runtime_config.hpl_nb_multiplier_threshold[nbnum] = tmpval;
					nbnum++;
				}
			}
		}
	}
	else if (strcmp(cmd, "HPL_PARAMDEFS") == 0)
	{
		int len = strlen(option);
		if (len)
		{
			if (strlen(global_runtime_config.paramdefs)) len++;
			len += strlen(global_runtime_config.paramdefs);
			global_runtime_config.paramdefs = (char*) realloc(global_runtime_config.paramdefs, len + 1);
			if (strlen(global_runtime_config.paramdefs)) strcat(global_runtime_config.paramdefs, " ");
			strcat(global_runtime_config.paramdefs, option);
		}
	}
	else
	{
		HPL_fprintf(TEST->outfp, "Unknown HPL Runtime option: %s\n", cmd);
	}
  }
  free(Buffer);

	char* envPtr;
	if ((envPtr = getenv("HPL_WARMUP")))
	{
		global_runtime_config.warmup = atoi(envPtr);
	}
	if ((envPtr = getenv("HPL_FASTRAND")))
	{
		global_runtime_config.fastrand = atoi(envPtr);
	}
	if ((envPtr = getenv("HPL_DISABLE_LOOKAHEAD")))
	{
		global_runtime_config.disable_lookahead = atoi(envPtr);
	}
	if ((envPtr = getenv("HPL_LOOKAHEAD2_TURNOFF")))
	{
		global_runtime_config.lookahead2_turnoff = atoi(envPtr);
	}
	if ((envPtr = getenv("HPL_DURATION_FIND_HELPER")))
	{
		global_runtime_config.duration_find_helper = atoi(envPtr);
	}
	if ((envPtr = getenv("HPL_CALDGEMM_ASYNC_FACT_DGEMM")))
	{
		global_runtime_config.caldgemm_async_fact_dgemm = atoi(envPtr);
	}
	if ((envPtr = getenv("HPL_CALDGEMM_ASYNC_FACT_FIRST")))
	{
		global_runtime_config.caldgemm_async_fact_first = atoi(envPtr);
	}
	if ((envPtr = getenv("HPL_CALDGEMM_ASYNC_DTRSM")))
	{
		global_runtime_config.caldgemm_async_dtrsm = atoi(envPtr);
	}
	if ((envPtr = getenv("HPL_CALDGEMM_ASYNC_FACT_DTRSM")))
	{
		global_runtime_config.caldgemm_async_fact_dtrsm = atoi(envPtr);
	}
	if ((envPtr = getenv("HPL_PARAMDEFS")))
	{
		int len = strlen(envPtr);
		if (len)
		{
			if (strlen(global_runtime_config.paramdefs)) len++;
			len += strlen(global_runtime_config.paramdefs);
			global_runtime_config.paramdefs = (char*) realloc(global_runtime_config.paramdefs, len + 1);
			if (strlen(global_runtime_config.paramdefs)) strcat(global_runtime_config.paramdefs, " ");
			strcat(global_runtime_config.paramdefs, envPtr);
		}
	}
#endif

    if (global_runtime_config.hpl_nb_multiplier_count)
    {
        max_gpu_nb_factor = 1;
	printf("NB Multipliers:");
	for (i = 0;i < global_runtime_config.hpl_nb_multiplier_count;i++)
	{
		printf(" %d/%d", global_runtime_config.hpl_nb_multiplier_threshold[i], global_runtime_config.hpl_nb_multiplier_factor[i]);
		if (global_runtime_config.hpl_nb_multiplier_factor[i] > max_gpu_nb_factor) max_gpu_nb_factor = global_runtime_config.hpl_nb_multiplier_factor[i];
	}
	max_gpu_nb *= max_gpu_nb_factor;
	
	printf(" MaxNB: %d\n", max_gpu_nb);
    }
/*
 * End of HPL_pdinfo
 */
}
