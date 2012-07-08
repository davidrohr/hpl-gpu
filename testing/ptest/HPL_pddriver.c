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
#include "util_trace.h"
#include "util_cal.h"

int HPL_init_laswp(void* ptr);

#if defined(HPL_INTERLEAVE_MEMORY) & !defined(HPL_CALL_CALDGEMM)
#define MPOL_DEFAULT 0
#define MPOL_PREFERRED 1
#define MPOL_BIND 2
#define MPOL_INTERLEAVE 3
#include <syscall.h>
#endif

int main
(
   int                        ARGC,
   char                       * * ARGV
)
{
/* 
 * Purpose
 * =======
 *
 * main is the main driver program for testing the HPL routines.
 * This  program is  driven  by  a short data file named  "HPL.dat".
 *
 * ---------------------------------------------------------------------
 */ 
/*
 * .. Local Variables ..
 */
   int                        nval  [HPL_MAX_PARAM],
                              nbval [HPL_MAX_PARAM],
                              pval  [HPL_MAX_PARAM],
                              qval  [HPL_MAX_PARAM],
                              nbmval[HPL_MAX_PARAM],
                              ndvval[HPL_MAX_PARAM],
                              ndhval[HPL_MAX_PARAM];

   HPL_T_FACT                 pfaval[HPL_MAX_PARAM],
                              rfaval[HPL_MAX_PARAM];

   HPL_T_TOP                  topval[HPL_MAX_PARAM];

   HPL_T_grid                 grid;
   HPL_T_palg                 algo;
   HPL_T_test                 test;
   int                        align, in, inb,
                              inbm, indh, indv, ipfa, ipq, irfa, itop,
                              mycol, myrow, ns, nbs, nbms, ndhs, ndvs,
                              npcol, npfs, npqs, nprow, nrfs, ntps, 
                              rank, size, seed;
   HPL_T_ORDER                pmapping;
   HPL_T_FACT                 rpfa;

   int run;
/* ..
 * .. Executable Statements ..
 */
#ifdef HPL_NO_MPI_THREAD_CHECK
   MPI_Init( &ARGC, &ARGV );
#else
   int mpiavail = 0;
   
#ifdef HPL_MPI_FUNNELED_THREADING
#define MPI_REQUIRE_THREAD_SAFETY MPI_THREAD_FUNNELED
#else
#define MPI_REQUIRE_THREAD_SAFETY MPI_THREAD_SERIALIZED
#endif

#if defined(HPL_INTERLEAVE_MEMORY) & !defined(HPL_CALL_CALDGEMM)
   unsigned long nodemask = 0xffffff;
   syscall(SYS_set_mempolicy, MPOL_INTERLEAVE, &nodemask, sizeof(nodemask) * 8);
   #error a
#endif

    if (MPI_Init_thread( &ARGC, &ARGV, MPI_REQUIRE_THREAD_SAFETY, &mpiavail ) != MPI_SUCCESS)
   {
	printf("Error initializing MPI\n");
	return(1);
   }
   if (mpiavail != MPI_REQUIRE_THREAD_SAFETY)
   {
	printf("MPI does not provide the required thread safety\n");
	return(1);
   }
#endif

#ifdef HPL_GPU_FACTORIZE
    double* A = (double*) malloc(1024);
    double* B = (double*) malloc(1024);
    dtrsm('L', 'U', 'T', 'U', 1, 1, 1.0, A, 1, B, 1);
    free(A);
    free(B);
#endif

#ifdef HPL_CALL_CALDGEMM
   if (CALDGEMM_Init())
   {
	printf("Error initializing CALDGEMM, abborting run\n");
	return(1);
   }
#endif

#ifndef USE_ORIGINAL_LASWP
#ifdef HPL_CALL_CALDGEMM
	HPL_init_laswp(CALDGEMM_GetObject());
#else
	HPL_init_laswp(NULL);
#endif
#endif

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   MPI_Comm_size( MPI_COMM_WORLD, &size );
/*
HPLinpack benchmark input file
6            device out (6=stdout,7=stderr,file)
1            # of problems sizes (N)
89088        Ns
1            # of NBs
1024         NBs
0            PMAP process mapping (0=Row-,1=Column-major)
1            # of process grids (P x Q)
1            Ps
1            Qs
0.1          threshold
1            # of panel fact
1            PFACTs (0=left, 1=Crout, 2=Right)
1            # of recursive stopping criterium
64           NBMINs (>= 1)
1            # of panels in recursion
2            NDIVs
1            # of recursive panel fact.
0            RFACTs (0=left, 1=Crout, 2=Right)
1            # of broadcast
5            BCASTs (0=1rg,1=1rM,2=2rg,3=2rM,4=Lng,5=LnM,6=MPI)
1            # of lookahead options
1            LOOKAHEADs (enable = 1)
8            memory alignment in double (> 0)
100          Seed for the matrix generation 
 */
   HPL_pdinfo( &test, &ns, nval, &nbs, nbval, &pmapping, &npqs, pval, qval,
               &npfs, pfaval, &nbms, nbmval, &ndvs, ndvval, &nrfs, rfaval,
               &ntps, topval, &ndhs, ndhval, &align, &seed );
/*
 * Loop over different process grids - Define process grid. Go to bottom
 * of process grid loop if this case does not use my process.
 */
   run = 0;
   for( ipq = 0; ipq < npqs; ipq++ )
   {
      (void) HPL_grid_init( MPI_COMM_WORLD, pmapping, pval[ipq], qval[ipq],
                            &grid );
      (void) HPL_grid_info( &grid, &nprow, &npcol, &myrow, &mycol );

      if( ( myrow < 0 ) || ( myrow >= nprow ) ||
          ( mycol < 0 ) || ( mycol >= npcol ) ) goto label_end_of_npqs;

#ifdef HPL_CALL_CALDGEMM
      CALDGEMM_set_num_nodes(pval[ipq] * qval[ipq], grid.iam);
#endif

      for( in = 0; in < ns; in++ )
      {                            /* Loop over various problem sizes */
       for( inb = 0; inb < nbs; inb++ )
       {                        /* Loop over various blocking factors */
        for( indh = 0; indh < ndhs; indh++ )
        {                       /* Loop over various lookahead depths */
         for( itop = 0; itop < ntps; itop++ )
         {                  /* Loop over various broadcast topologies */
          for( irfa = 0; irfa < nrfs; irfa++ )
          {             /* Loop over various recursive factorizations */
           for( ipfa = 0; ipfa < npfs; ipfa++ )
           {                /* Loop over various panel factorizations */
            for( inbm = 0; inbm < nbms; inbm++ )
            {        /* Loop over various recursive stopping criteria */
             for( indv = 0; indv < ndvs; indv++ )
             {          /* Loop over various # of panels in recursion */
/*
 * Set up the algorithm parameters
 */
              algo.btopo = topval[itop]; algo.depth = ndhval[indh];
              algo.nbmin = nbmval[inbm]; algo.nbdiv = ndvval[indv];

              algo.pfact = rpfa = pfaval[ipfa];

              if( rpfa == HPL_LEFT_LOOKING ) algo.pffun = HPL_pdpanllT;
              else if( rpfa == HPL_CROUT   ) algo.pffun = HPL_pdpancrT;
              else                           algo.pffun = HPL_pdpanrlT;

              algo.rfact = rpfa = rfaval[irfa];
              if( rpfa == HPL_LEFT_LOOKING ) algo.rffun = HPL_pdrpanllT;
              else if( rpfa == HPL_CROUT   ) algo.rffun = HPL_pdrpancrT;
              else                           algo.rffun = HPL_pdrpanrlT;

              algo.align = align;

              ++run;
#ifdef TRACE_CALLS
              resetTraceCounters();
#endif

			  fprintfct(stderr, "(Problem: N %d NB %d)(Network: BCAST %d LOOKAHEAD %d) (Factorization: NBMIN %d NBDIV %d PFACT %d RFACT %d)\n", nval[in], nbval[inb], algo.btopo, algo.depth, algo.nbmin, algo.nbdiv, algo.pfact, algo.rfact);
			  
   int mcols = (nval[in] + nbval[inb]) / nbval[inb];
   grid.col_mapping = (int*) malloc(mcols * sizeof(int));
   grid.mcols_per_pcol = (int*) malloc(npcol * sizeof(int));
   
   if (rank == 0)
   {
      float* cols = malloc(npcol * sizeof(float));
      int nprocs = npcol * nprow;
      for (int i = 0;i < npcol;i++) cols[i] = 1.;
      for (int i = 0;i < nprocs;i++)
      {
         if (pmapping == HPL_ROW_MAJOR)
         {
            fprintfctd(test.outfp, "Node %d col %d perf %f/%f\n", i, i % npcol, test.node_perf[i], cols[i % npcol]);
            if (test.node_perf[i] < cols[i % npcol]) cols[i % npcol] = test.node_perf[i];
         }
         else
         {
            fprintfctd(test.outfp, "Node %d col %d perf %f/%f\n", i, i / nprow, test.node_perf[i], cols[i / nprow]);
            if (test.node_perf[i] < cols[i / nprow]) cols[i / nprow] = test.node_perf[i];
         }
      }
      float max_perf = 0;
      for (int i = 0;i < npcol;i++) if (cols[i] > max_perf) max_perf = cols[i];
      for (int i = 0;i < npcol;i++) cols[i] /= max_perf;
      max_perf = 0;
      for (int i = 0;i < npcol;i++) max_perf += cols[i];
      for (int i = 0;i < npcol;i++)
      {
         fprintfctd(test.outfp, "Process Col %d Performance %f (of %f total)\n", i, cols[i], max_perf);
      }
      
      for (int i = 0;i < npcol;i++) grid.mcols_per_pcol[i] = 0;
      int j = 0;
      int lastcol = -1;
      float relax = 0;
      for (int i = 0;i < mcols;i++)
      {
         if (npcol == 1)
         {
            grid.col_mapping[i] = i % npcol;
            grid.mcols_per_pcol[i % npcol]++;
	    continue;
         }
         int jstart = j;
         int round1 = 1;
         while (i && (j == lastcol || cols[j] / max_perf * (float) (i + 1) < (float) grid.mcols_per_pcol[j] + 0.5 * (float) round1 - relax))
         {
            fprintfctd(test.outfp, "Skipping process col %d (desired mcols %f, present mcols %d)\n", j, cols[j] / max_perf * (float) (i + 1), grid.mcols_per_pcol[j]);
            j++;
            j = j % npcol;
            if (j == jstart)
            {
               if (round1 > 0) round1 = 0;
               else relax += 0.1;
            }
         }
         grid.col_mapping[i] = j;
         grid.mcols_per_pcol[j]++;
         lastcol = j;
         relax = 0;
         fprintfctd(test.outfp, "Matrix col %d processed by process col %d (%d total matrix cols)\n", i, j, grid.mcols_per_pcol[j]);
         j++;
         j = j % npcol;
      }
      
      for (int i = 0;i < npcol;i++)
      {
         fprintfct(test.outfp, "Process col %d process %d matrix cols\n", i, grid.mcols_per_pcol[i]);
      }
      
      free(cols);
   }
   
   MPI_Bcast(grid.col_mapping, mcols, MPI_INT, 0, grid.all_comm);
   MPI_Bcast(grid.mcols_per_pcol, npcol, MPI_INT, 0, grid.all_comm);
			  

              HPL_pdtest( &test, &grid, &algo, nval[in], nbval[inb], seed );
              free(grid.col_mapping);
              free(grid.mcols_per_pcol);

#ifdef TRACE_CALLS
              writeTraceCounters( "trace_counters", run, rank );
#endif
             }
            }
           }
          }
         }
        }
       }
      }
      (void) HPL_grid_exit( &grid );
label_end_of_npqs: ;
   }
/*
 * Print ending messages, close output file, exit.
 */
   free(test.node_perf);
   if( rank == 0 )
   {
      test.ktest = test.kpass + test.kfail + test.kskip;
#ifndef HPL_DETAILED_TIMING
      HPL_fprintf( test.outfp, "%s%s\n",
                   "========================================",
                   "========================================" );
#else
      if( test.thrsh > HPL_rzero )
         HPL_fprintf( test.outfp, "%s%s\n",
                      "========================================",
                      "========================================" );
#endif

      HPL_fprintf( test.outfp, "\n%s %6d %s\n", "Finished", test.ktest,
                   "tests with the following results:" );
      if( test.thrsh > HPL_rzero )
      {
         HPL_fprintf( test.outfp, "         %6d %s\n", test.kpass,
                      "tests completed and passed residual checks," );
         HPL_fprintf( test.outfp, "         %6d %s\n", test.kfail,
                      "tests completed and failed residual checks," );
         HPL_fprintf( test.outfp, "         %6d %s\n", test.kskip,
                      "tests skipped because of illegal input values." );
      }
      else
      {
         HPL_fprintf( test.outfp, "         %6d %s\n", test.kpass,
                      "tests completed without checking," );
         HPL_fprintf( test.outfp, "         %6d %s\n", test.kskip,
                      "tests skipped because of illegal input values." );
      }

      HPL_fprintf( test.outfp, "%s%s\n",
                   "----------------------------------------",
                   "----------------------------------------" );
      HPL_fprintf( test.outfp, "\nEnd of tests.\n" );
      HPL_fprintf( test.outfp, "%s%s\n",
                   "========================================",
                   "========================================" );

      if( ( test.outfp != stdout ) && ( test.outfp != stderr ) )
         (void) fclose( test.outfp );
   }
#ifdef HPL_CALL_CALDGEMM
   CALDGEMM_Shutdown();
#endif
#ifdef TRACE_CALLS
   releaseTraceCounters();
#endif
   MPI_Finalize();
   exit( 0 );

   return( 0 );
/*
 * End of main
 */
}
