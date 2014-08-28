/**
 * C wrapper for CALDGEMM
 *
 * Copyright 2010:
 *  - David Rohr (drohr@jwdt.org)
 *  - Matthias Bach (bach@compeng.uni-frankfurt.de)
 *  - Matthias Kretz (kretz@compeng.uni-frankfurt.de)
 *
 * This file is part of HPL-GPU.
 *
 * HPL-GPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HPL-GPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HPL-GPU.  If not, see <http://www.gnu.org/licenses/>.
 *
 * In addition to the rules layed out by the GNU General Public License
 * the following exception is granted:
 *
 * Use with the Original BSD License.
 *
 * Notwithstanding any other provision of the GNU General Public License
 * Version 3, you have permission to link or combine any covered work with
 * a work licensed under the 4-clause BSD license into a single combined
 * work, and to convey the resulting work.  The terms of this License will
 * continue to apply to the part which is the covered work, but the special
 * requirements of the 4-clause BSD license, clause 3, concerning the
 * requirement of acknowledgement in advertising materials will apply to
 * the combination as such.
 */
 
#ifdef HPL_CALL_CALDGEMM

#define mcat(a, b) a ## b
#define mxcat(a, b) mcat(a, b)

#ifndef HPL_CALDGEMM_BACKEND
#define HPL_CALDGEMM_BACKEND cal
#endif

#define CALDGEMM_IMPL mxcat(caldgemm_, HPL_CALDGEMM_BACKEND)

#define CALDGEMM_HEADER < CALDGEMM_IMPL.h >

#include CALDGEMM_HEADER
extern "C"
{

	typedef unsigned int blasint;
#include <cblas.h>
	int max_gpu_nb = 1024;
	extern int HPL_CALDGEMM_gpu_height;
}
#include <pthread.h>
#include <errno.h>
#include "util_cal.h"

#define fprintfdvv( a, b )	//Disable verbose verbose debug output
//#define fprintfdvv fprintf


extern "C"
{
	extern volatile size_t HPL_CALDGEMM_swap_current_n;
	extern void HPL_CALDGEMM_wrapper_factorize();
	extern void HPL_CALDGEMM_wrapper_broadcast();
	extern void HPL_CALDGEMM_wrapper_swap();
}

static caldgemm::caldgemm_config cal_info;
static caldgemm* cal_dgemm;
char PreOutput[64] = "";

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
	fprintfdvv(STD_OUT, "GPU DGEMM Thread Started\n");
	pthread_mutex_lock(&startgpudgemm);
	fprintfdvv(STD_OUT, "GPU DGEMM Thread waiting for commands\n");
	while (pthread_mutex_lock(&startgpudgemm) == 0 && exitgpudgemm == 0)
	{
		fprintfdvv(STD_OUT, "GPU DGEMM Thread Running\n");
		cal_dgemm->RunCALDGEMM( gpudgemmparams.A, gpudgemmparams.B, gpudgemmparams.C, gpudgemmparams.ALPHA, gpudgemmparams.BETA, gpudgemmparams.M, gpudgemmparams.K, gpudgemmparams.N, gpudgemmparams.LDA, gpudgemmparams.LDB,
			gpudgemmparams.LDC, gpudgemmparams.ORDER == CblasColMajor, gpudgemmparams.TRANSA == CblasTrans, gpudgemmparams.TRANSB == CblasTrans, gpudgemmparams.LinpackCallbacks );
		pthread_mutex_unlock(&gpudgemmdone);
	}
	fprintfdvv(STD_OUT, "GPU DGEMM Thread Terminating\n");
	pthread_mutex_unlock(&gpudgemmdone);
	
	return(NULL);
}

void funneled_factorize_wrapper()
{
	fprintfdvv(STD_OUT, "Factorize Funneled Wrapper\n");
	pthread_mutex_unlock(&startfactorize);
	pthread_mutex_lock(&factorizedone);
	fprintfdvv(STD_OUT, "Factorize Funneled Wrapper Ended\n");
}

void funneled_broadcast_wrapper()
{
        fprintfdvv(STD_OUT, "Broadcast Funneled Wrapper\n");
	pthread_mutex_unlock(&startbroadcast);
	pthread_mutex_lock(&broadcastdone);
        fprintfdvv(STD_OUT, "Broadcast Funneled Wrapper Ended\n");
}

#endif

void* CALDGEMM_GetObject()
{
	return(cal_dgemm);
}

#ifdef HPL_EMULATE_MULTINODE
extern "C" int HPL_CALDGEMM_wrapper_n;

float* broadcast_fake_array;

void multinode_broadcast_fake()
{
	for (int j = 0;j < HPL_CALDGEMM_wrapper_n / 32;j++)
	{
		for (int i = 0;i < 1024 * 1024 * 2;i++)
		{
			broadcast_fake_array[i] += 1.;
		}
	}
}
#endif

#ifdef HPL_RESTRICT_CALLBACK
int HPL_Restrict_Callback_Function(int matrix_n)
{
	return(HPL_RESTRICT_CALLBACK(matrix_n));
}
#endif

int CALDGEMM_Init()
{
#ifdef HPL_GPU_VERIFY
	cal_info.Verify = true;
#endif

	//cal_info.Disassemble = false;
	cal_info.TabularTiming = true;

#if defined(TRACE_CALLS) | defined(HPL_GPU_NOT_QUIET)
	cal_info.Quiet = false;
#else
	cal_info.Quiet = true;
#endif

#ifdef HPL_GPU_TIMING
	cal_info.DisplayTiming = true;
#endif

#ifdef HPL_MULTI_GPU
	cal_info.DeviceNum = -1;
	cal_info.Height = 3072;
	cal_info.ImprovedScheduler = true;
	cal_info.MultiThreadDivide = true;
#else
	cal_info.DeviceNum = 0;
	//cal_info.Height = 4096;
#endif

#ifdef HPL_GPU_THREADSAVE_DRIVER
	cal_info.ThreadSaveDriver = 1;
#elif defined(HPL_GPU_GLOBAL_DRIVER_MUTEX)
	cal_info.ThreadSaveDriver = -1;
#endif

#ifdef HPL_GPU_MAPPING
	const int mapping[] = HPL_GPU_MAPPING;
	for (unsigned int i = 0;i < sizeof(mapping) / sizeof(int);i++)
	{
	    cal_info.GPUMapping[i] = mapping[i];
	}
#endif
#ifdef HPL_GPU_POSTPROCESS_MAPPING
	const int postprocess_mapping[] = HPL_GPU_POSTPROCESS_MAPPING;
	for (unsigned int i = 0;i < sizeof(postprocess_mapping) / sizeof(int);i++)
	{
	    cal_info.PostprocessMapping[i] = postprocess_mapping[i];
	}
#endif
#ifdef HPL_GPU_ALLOC_MAPPING
	const int alloc_mapping[] = HPL_GPU_ALLOC_MAPPING;
	for (unsigned int i = 0;i < sizeof(alloc_mapping) / sizeof(int);i++)
	{
	    cal_info.AllocMapping[i] = alloc_mapping[i];
	}
#endif
#ifdef HPL_GPU_DMA_MAPPING
	const int dma_mapping[] = HPL_GPU_DMA_MAPPING;
	for (unsigned int i = 0;i < sizeof(dma_mapping) / sizeof(int);i++)
	{
	    cal_info.DMAMapping[i] = dma_mapping[i];
	}
#endif
#ifdef HPL_GPU_EXCLUDE_CORES
	static const int exclude_cores[] = HPL_GPU_EXCLUDE_CORES;
	cal_info.nExcludeCPUCores = sizeof(exclude_cores) / sizeof(int);
	cal_info.ExcludeCPUCores = (int*) exclude_cores;
#endif
#ifdef HPL_GPU_DEVICE_IDS
	const int device_ids[] = HPL_GPU_DEVICE_IDS;
	for (unsigned int i = 0;i < sizeof(device_ids) / sizeof(int);i++)
	{
	    cal_info.DeviceNums[i] = device_ids[i];
	}
	cal_info.NumDevices = sizeof(device_ids) / sizeof(int);
#endif

#ifdef HPL_GPU_PIN_MAIN
	cal_info.PinMainThread = HPL_GPU_PIN_MAIN;
#endif
#ifdef HPL_MPI_AFFINITY
	const int mpi_affinity[] = HPL_MPI_AFFINITY;
	cal_info.PinBroadcastThread = mpi_affinity[0];
#endif

#ifdef HPL_GPU_OG
	cal_info.DstMemory = 'g';
	cal_info.ImplicitDriverSync = 1;
	cal_info.KeepBuffersMapped = false;
#else
#ifdef HPL_NO_HACKED_LIB
	cal_info.KeepBuffersMapped = false;
#else
	cal_info.KeepBuffersMapped = true;
#endif
#endif

#ifndef HPL_GPU_MAX_NB
	cal_info.Width = max_gpu_nb; //k for matrix multiply
#else
	cal_info.Width = HPL_GPU_MAX_NB;
#endif

#ifdef HPL_FAST_GPU
	cal_info.SmallTiles = 1;
	cal_info.GPURatio = 1.0;
	cal_info.DynamicSched = false;
#endif

#ifdef HPL_SLOW_CPU
	cal_info.SlowCPU = true;
#endif

#ifndef HPL_RESTRICT_CPUS
#define HPL_RESTRICT_CPUS 2
#endif

#ifdef HPL_INITIAL_GPU_RATIO
	cal_info.GPURatio = -HPL_INITIAL_GPU_RATIO;
#endif

#ifdef HPL_ALTERNATE_LOOKAHEAD
	cal_info.AlternateLookahead = HPL_ALTERNATE_LOOKAHEAD;
#endif

	cal_info.HPLFactorizeRestrictCPUs = HPL_RESTRICT_CPUS;
#ifdef HPL_RESTRICT_CALLBACK
	cal_info.HPLFactorizeRestrictCallback = HPL_Restrict_Callback_Function;
#endif
	
	//cal_info.AutoHeight = true;
	//cal_info.Iterations = 1;
	//cal_info.VerboseTiming = false;
	//cal_info.Debug = true;
	//cal_info.MultiThread = true;
	//cal_info.UseGPU = true;
	//cal_info.UseCPU = true;
	//cal_info.GPURatio = -1;
	//cal_info.DynamicSched = false;
	//cal_info.DumpMatrix = false;
	
#ifdef HPL_INTERLEAVE_MEMORY
	cal_info.MemPolicy = true;
#else
	cal_info.MemPolicy = false;
#endif

#ifdef HPL_GPU_PERFORMANCE_WARNINGS
	cal_info.NoPerformanceWarnings = false;
#else
	cal_info.NoPerformanceWarnings = true;
#endif
	//cal_info.AsyncTiming = (CALboolean) !cal_info.NoPerformanceWarnings;
	
	cal_info.linpack_swap_function = HPL_CALDGEMM_wrapper_swap;
	cal_info.PreOut = PreOutput;

#ifdef HPL_MPI_FUNNELED_THREADING
	cal_info.linpack_factorize_function = funneled_factorize_wrapper;
	cal_info.linpack_broadcast_function = funneled_broadcast_wrapper;
#else	
	cal_info.linpack_factorize_function = HPL_CALDGEMM_wrapper_factorize;
#ifdef HPL_EMULATE_MULTINODE
	cal_info.linpack_broadcast_function = multinode_broadcast_fake;
	broadcast_fake_array = new float[1024 * 1024 * 2];
#elif !defined(HPL_NO_MPI_LIB)
	cal_info.linpack_broadcast_function = HPL_CALDGEMM_wrapper_broadcast;
#endif
#endif

#ifdef HPL_PRINT_THROTTLING_NODES
	cal_info.GPUClock = HPL_PRINT_THROTTLING_NODES;
#endif

#ifdef HPL_GPU_EXTRA_CALDGEMM_OPTIONS
HPL_GPU_EXTRA_CALDGEMM_OPTIONS
#endif

#ifdef HPL_GPU_FACTORIZE
	bool a = true;
#else
	bool a = false;
#endif
	cal_dgemm = new CALDGEMM_IMPL;
	if (cal_dgemm == NULL) return(1);
	cal_info.config_backend = cal_dgemm->create_caldgemm_config_backend();
#ifdef HPL_GPU_EXTRA_CALDGEMM_BACKEND_OPTIONS
	HPL_GPU_EXTRA_CALDGEMM_BACKEND_OPTIONS
#endif
	int retVal = cal_dgemm->InitCALDGEMM( &cal_info, a );
	
#ifdef HPL_MPI_FUNNELED_THREADING
	pthread_mutex_init(&startgpudgemm, NULL);
	pthread_mutex_init(&gpudgemmdone, NULL);
	pthread_mutex_init(&startfactorize, NULL);
	pthread_mutex_init(&factorizedone, NULL);
	pthread_mutex_init(&startbroadcast, NULL);
	pthread_mutex_init(&broadcastdone, NULL);
	pthread_t thr;
	fprintfdvv(STD_OUT, "Starting wrapper thread\n");
	pthread_create(&thr, NULL, gpudgemm_wrapper, NULL);
	while (pthread_mutex_trylock(&startgpudgemm) != EBUSY) pthread_mutex_unlock(&startgpudgemm);
	fprintfdvv(STD_OUT, "Wrapper started sucessfully\n");
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
	cal_dgemm->ExitCALDGEMM();
	delete cal_dgemm;
}

void CALDGEMM_dgemm( const enum CBLAS_ORDER ORDER, const enum CBLAS_TRANSPOSE TRANSA,
                     const enum CBLAS_TRANSPOSE TRANSB, const int M, const int N,
                     const int K, const double ALPHA, const double * A, const int LDA,
                     const double * B, const int LDB, const double BETA, double * C,
                     const int LDC, int LinpackCallbacks )
{
	static int LinpackIteration = 0;
	if (M == 0 || N == 0 || K == 0)
	{
            if (cal_info.LinpackSwapN != NULL)
            {
        	HPL_CALDGEMM_gpu_height = 0;
        	HPL_CALDGEMM_wrapper_swap();
            }
            if (LinpackCallbacks)
            {
        	HPL_CALDGEMM_wrapper_factorize();
        	if (cal_info.LinpackNodes > 1) HPL_CALDGEMM_wrapper_broadcast();
            }
	    return;
	}
        else if (K >= 512 && (M >= 2048 || N >= 2048))
        {
            sprintf(PreOutput, "#(%-3d,%4d) ", cal_info.MPIRank, LinpackIteration++);
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
	    fprintfdvv(STD_OUT, "Running GPU dgemm\n");
	    if (pthread_mutex_unlock(&startgpudgemm)) fprintf(STD_OUT, "Mutex Error: %s - %d\n", __FILE__, __LINE__);
	    if (LinpackCallbacks)
	    {
		    fprintfdvv(STD_OUT, "Waiting to factorize\n");
		    pthread_mutex_lock(&startfactorize);
		    HPL_CALDGEMM_wrapper_factorize();
		    pthread_mutex_unlock(&factorizedone);
		    
		    if (cal_info.LinpackNodes > 1)
		    {
			    cpu_set_t old_mask, linpack_mask;
			    CPU_ZERO(&linpack_mask);
			    CPU_SET(cal_dgemm->broadcastcore(), &linpack_mask);
			    sched_getaffinity(0, sizeof(cpu_set_t), &old_mask);
			    sched_setaffinity(0, sizeof(cpu_set_t), &linpack_mask);		    
		    
			    fprintfdvv(STD_OUT, "Waiting to broadcast\n");
			    pthread_mutex_lock(&startbroadcast);
			    HPL_CALDGEMM_wrapper_broadcast();
			    pthread_mutex_unlock(&broadcastdone);
			    sched_setaffinity(0, sizeof(cpu_set_t), &old_mask);
		    }
		    fprintfdvv(STD_OUT, "Factorize and broadcast done\n");
	    }
	    
	    pthread_mutex_lock(&gpudgemmdone);
#else
	    if (cal_dgemm->RunCALDGEMM( (double*) A, (double*) B, C, (double) ALPHA, (double) BETA, (int) M, (int) K, (int) N, (int) LDA, (int) LDB, (int) LDC, ORDER == CblasColMajor, TRANSA == CblasTrans, TRANSB == CblasTrans, LinpackCallbacks ))
	    {
		printf("Error in CALDGEMM Run, aborting HPL Run\n");
		exit(1);
	    }
#endif
        }
        else
        {
            // Use plain cblas
            if (cal_info.LinpackSwapN != NULL)
            {
        	HPL_CALDGEMM_gpu_height = 0;
        	HPL_CALDGEMM_wrapper_swap();
            }
            cblas_dgemm( ORDER, TRANSA, TRANSB, M, N, K, ALPHA, (double*) A, LDA, (double*) B, LDB, BETA, C, LDC );
            if (LinpackCallbacks)
            {
        	HPL_CALDGEMM_wrapper_factorize();
        	if (cal_info.LinpackNodes > 1) HPL_CALDGEMM_wrapper_broadcast();
            }
        }
}

void* CALDGEMM_alloc(size_t size, int interleaved)
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
    return((void*) cal_dgemm->AllocMemory(size / sizeof(double), page_locked, huge_tables, false, interleaved));
}

void CALDGEMM_free(void* ptr)
{
    cal_dgemm->FreeMemory((double*) ptr);
}

void CALDGEMM_set_num_nodes(int num, int rank)
{
    cal_info.LinpackNodes = num;
    cal_info.MPIRank = rank;
}
void CALDGEMM_enable_async_laswp(int enable)
{
	cal_info.LinpackSwapN = enable ? &HPL_CALDGEMM_swap_current_n : NULL;
}
#else
#include <stdlib.h>
extern "C" void CALDGEMM_free(void* ptr);
extern "C" void* CALDGEMM_alloc(size_t size, int interleaved);

void CALDGEMM_free(void* ptr)
{
    free(ptr);
}
void* CALDGEMM_alloc(size_t size, int interleaved)
{
    return malloc(size);
}
#endif
