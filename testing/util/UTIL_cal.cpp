/**
 * C wrapper for CALDGEMM
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

#include <caldgemm.h>
#include <pthread.h>
#include <errno.h>
#include "util_cal.h"

#define fprintfdvv	//Disable verbose verbose debug output
#define fprintfdvv fprintf


extern "C"
{
    extern void HPL_CALLDGEMM_wrapper_factorize();
    extern void HPL_CALLDGEMM_wrapper_broadcast();
}

static caldgemm::SampleInfo cal_info;
static caldgemm cal_dgemm;

#ifdef HPL_MPI_FUNNELED_THREADING
pthread_mutex_t startgpudgemm, gpudgemmdone, startfactorize, factorizedone, startbroadcast, broadcastdone;
int exitgpudgemm = 0;

struct gpudgemmparamstruct
{
	double *A, *B, *C;
	double ALPHA, BETA;
	int M, K, N, LDA, LDB, LDC;
	CBLAS_ORDER ORDER;
	CBLAS_TRANSPOSE TRANSA, TRANSB;
	int LinpackCallbacks;
} gpudgemmparams;

void* gpudgemm_wrapper(void* arg)
{
	fprintfdvv(stderr, "GPU DGEMM Thread Started\n");
	pthread_mutex_lock(&startgpudgemm);
	fprintfdvv(stderr, "GPU DGEMM Thread waiting for commands\n");
	while (pthread_mutex_lock(&startgpudgemm) == 0 && exitgpudgemm == 0)
	{
		fprintfdvv(stderr, "GPU DGEMM Thread Running\n");
		cal_dgemm.RunCALDGEMM( gpudgemmparams.A, gpudgemmparams.B, gpudgemmparams.C, gpudgemmparams.ALPHA, gpudgemmparams.BETA, gpudgemmparams.M, gpudgemmparams.K, gpudgemmparams.N, gpudgemmparams.LDA, gpudgemmparams.LDB,
			gpudgemmparams.LDC, gpudgemmparams.ORDER, gpudgemmparams.TRANSA, gpudgemmparams.TRANSB, gpudgemmparams.LinpackCallbacks );
		pthread_mutex_unlock(&gpudgemmdone);
	}
	fprintfdvv(stderr, "GPU DGEMM Thread Terminating\n");
	pthread_mutex_unlock(&gpudgemmdone);
	
	return(NULL);
}

void funneled_factorize_wrapper()
{
	fprintfdvv(stderr, "Factorize Funneled Wrapper\n");
	pthread_mutex_unlock(&startfactorize);
	pthread_mutex_lock(&factorizedone);
	fprintfdvv(stderr, "Factorize Funneled Wrapper Ended\n");
}

void funneled_broadcast_wrapper()
{
        fprintfdvv(stderr, "Broadcast Funneled Wrapper\n");
	pthread_mutex_unlock(&startbroadcast);
	pthread_mutex_lock(&broadcastdone);
        fprintfdvv(stderr, "Broadcast Funneled Wrapper Ended\n");
}

#endif

int CALDGEMM_Init()
{
#ifdef HPL_GPU_VERIFY
	cal_info.Verify = CAL_TRUE;
#endif
	//cal_info.Disassemble = CAL_FALSE;
	cal_info.TabularTiming = CAL_TRUE;
#if defined(TRACE_CALLS) | defined(HPL_GPU_NOT_QUIET)
	cal_info.Quiet = CAL_FALSE;
#else
	cal_info.Quiet = CAL_TRUE;
#endif
#ifdef HPL_GPU_TIMING
	cal_info.DisplayTiming = CAL_TRUE;
#endif
	//cal_info.DeviceNum = 0;
	cal_info.Width = 1024; //k for matrix multiply
	//cal_info.Height = 4096;
	//cal_info.AutoHeight = CAL_TRUE;
	//cal_info.Iterations = 1;
	//cal_info.DstMemory = 'c';
	//cal_info.VerboseTiming = CAL_FALSE;
	//cal_info.Debug = CAL_TRUE;
	//cal_info.MultiThread = CAL_TRUE;
	//cal_info.UseGPU = CAL_TRUE;
	//cal_info.UseCPU = CAL_TRUE;
	//cal_info.GPURatio = -1;
	//cal_info.DynamicSched = CAL_TRUE;
	//cal_info.MemPolicy = CAL_TRUE;
	//cal_info.DumpMatrix = CAL_FALSE;
#ifdef HPL_NO_PERFORMANCE_WARNINGS
	cal_info.NoPerformanceWarnings = CAL_TRUE;
#else
	cal_info.NoPerformanceWarnings = CAL_TRUE;
#endif
	//cal_info.AsyncTiming = (CALboolean) !cal_info.NoPerformanceWarnings;
	cal_info.KeepBuffersMapped = CAL_TRUE;

#ifdef HPL_MPI_FUNNELED_THREADING
	cal_info.linpack_factorize_function = funneled_factorize_wrapper;
	cal_info.linpack_broadcast_function = funneled_broadcast_wrapper;
#else	
	cal_info.linpack_factorize_function = HPL_CALLDGEMM_wrapper_factorize;
	cal_info.linpack_broadcast_function = HPL_CALLDGEMM_wrapper_broadcast;
#endif

	int retVal = cal_dgemm.InitCALDGEMM( &cal_info );
	
#ifdef HPL_MPI_FUNNELED_THREADING
	pthread_mutex_init(&startgpudgemm, NULL);
	pthread_mutex_init(&gpudgemmdone, NULL);
	pthread_mutex_init(&startfactorize, NULL);
	pthread_mutex_init(&factorizedone, NULL);
	pthread_mutex_init(&startbroadcast, NULL);
	pthread_mutex_init(&broadcastdone, NULL);
	pthread_t thr;
	fprintfdvv(stderr, "Starting wrapper thread\n");
	pthread_create(&thr, NULL, gpudgemm_wrapper, NULL);
	while (pthread_mutex_trylock(&startgpudgemm) != EBUSY) pthread_mutex_unlock(&startgpudgemm);
	fprintfdvv(stderr, "Wrapper started sucessfully\n");
	pthread_mutex_lock(&startfactorize);
	pthread_mutex_lock(&startbroadcast);
	pthread_mutex_lock(&factorizedone);
	pthread_mutex_lock(&broadcastdone);
	pthread_mutex_lock(&gpudgemmdone);
#endif
	
	return(retVal);
}

void CALDGEMM_Shutdown()
{
#ifdef HPL_MPI_FUNNELED_THREADING
	exitgpudgemm = 1;
	pthread_mutex_unlock(&startgpudgemm);
	pthread_mutex_lock(&gpudgemmdone);
	pthread_mutex_unlock(&startgpudgemm);
	pthread_mutex_unlock(&gpudgemmdone);
	pthread_mutex_unlock(&startfactorize);
	pthread_mutex_unlock(&factorizedone);
	pthread_mutex_unlock(&startbroadcast);
	pthread_mutex_unlock(&broadcastdone);
	pthread_mutex_destroy(&startgpudgemm);
	pthread_mutex_destroy(&gpudgemmdone);
	pthread_mutex_destroy(&startfactorize);
	pthread_mutex_destroy(&factorizedone);
	pthread_mutex_destroy(&startbroadcast);
	pthread_mutex_destroy(&broadcastdone);
#endif
	
	cal_dgemm.ExitCALDGEMM();
}

void CALDGEMM_dgemm( const enum CBLAS_ORDER ORDER, const enum CBLAS_TRANSPOSE TRANSA,
                     const enum CBLAS_TRANSPOSE TRANSB, const int M, const int N,
                     const int K, const double ALPHA, const double * A, const int LDA,
                     const double * B, const int LDB, const double BETA, double * C,
                     const int LDC, int LinpackCallbacks )
{
	if (!LinpackCallbacks && (M == 0 || N == 0 || K == 0)) return;
        else if(LinpackCallbacks || ( M >= 2048 && N >= 2048 && K >= 512 ))
        {
#ifdef HPL_MPI_FUNNELED_THREADING
	    gpudgemmparams.A = (double*) A;
	    gpudgemmparams.B = (double*) B;
	    gpudgemmparams.C = (double*) C;
	    gpudgemmparams.ALPHA = (double) ALPHA;
	    gpudgemmparams.BETA = (double) BETA;
	    gpudgemmparams.M = (int) M;
	    gpudgemmparams.K = (int) K;
	    gpudgemmparams.N = (int) N;
	    gpudgemmparams.LDA = (int) LDA;
	    gpudgemmparams.LDB = (int) LDB;
	    gpudgemmparams.LDC = (int) LDC;
	    gpudgemmparams.ORDER = ORDER;
	    gpudgemmparams.TRANSA = TRANSA;
	    gpudgemmparams.TRANSB = TRANSB;
	    gpudgemmparams.LinpackCallbacks = LinpackCallbacks;
	    fprintfdvv(stderr, "Running GPU dgemm\n");
	    if (pthread_mutex_unlock(&startgpudgemm)) fprintf(STD_OUT, "Mutex Error: %s - %d\n", __FILE__, __LINE__);
	    if (LinpackCallbacks)
	    {
		    fprintfdvv(stderr, "Waiting to factorize\n");
		    pthread_mutex_lock(&startfactorize);
		    HPL_CALLDGEMM_wrapper_factorize();
		    pthread_mutex_unlock(&factorizedone);
		    
		    if (cal_info.LinpackNodes > 1)
		    {
			    cpu_set_t old_mask, linpack_mask;
			    CPU_ZERO(&linpack_mask);
			    CPU_SET(cal_dgemm.broadcastcore(), &linpack_mask);
			    sched_getaffinity(0, sizeof(cpu_set_t), &old_mask);
			    sched_setaffinity(0, sizeof(cpu_set_t), &linpack_mask);		    
		    
			    fprintfdvv(stderr, "Waiting to broadcast\n");
			    pthread_mutex_lock(&startbroadcast);
			    HPL_CALLDGEMM_wrapper_broadcast();
			    pthread_mutex_unlock(&broadcastdone);
			    sched_setaffinity(0, sizeof(cpu_set_t), &old_mask);
		    }
		    fprintfdvv(stderr, "Factorize and broadcast done\n");
	    }
	    
	    pthread_mutex_lock(&gpudgemmdone);
#else
	    if (cal_dgemm.RunCALDGEMM( (double*) A, (double*) B, C, (double) ALPHA, (double) BETA, (int) M, (int) K, (int) N, (int) LDA, (int) LDB, (int) LDC, ORDER, TRANSA, TRANSB, LinpackCallbacks ))
	    {
		printf("Error in CALDGEMM Run, aborting HPL Run\n");
		exit(1);
	    }
#endif
        }
        else
        {
            // Use plain cblas
            cblas_dgemm( ORDER, TRANSA, TRANSB, M, N, K, ALPHA, (double*) A, LDA, (double*) B, LDB, BETA, C, LDC );
        }

}

void* CALDGEMM_alloc(size_t size)
{
    if (size % sizeof(double)) size += sizeof(double);
#ifdef HPL_PAGELOCKED_MEM
    bool page_locked = true;
#else
    bool page_locked = false;
#endif

#ifdef HPL_HUGE_TABLES
    bool huge_tables = true;
#else
    bool huge_tables = false;
#endif
    return((void*) cal_dgemm.AllocMemory(size / sizeof(double), page_locked, huge_tables));
}

void CALDGEMM_free(void* ptr)
{
    cal_dgemm.FreeMemory((double*) ptr);
}

void CALDGEMM_set_num_nodes(int num)
{
    cal_info.LinpackNodes = num;
}
