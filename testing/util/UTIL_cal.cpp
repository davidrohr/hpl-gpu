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

#define xstr(s) str(s)
#define str(s) #s

#ifndef HPL_CALDGEMM_BACKEND
#define HPL_CALDGEMM_BACKEND cal
#endif

#define CALDGEMM_IMPL mxcat(caldgemm_, HPL_CALDGEMM_BACKEND)

#define CALDGEMM_HEADER < CALDGEMM_IMPL.h >

#include CALDGEMM_HEADER
#define HPL_DIAG CBLAS_DIAG
#define HPL_SIDE CBLAS_SIDE
#define HPL_UPLO CBLAS_UPLO
#define HPL_TRANS CBLAS_TRANSPOSE
#define HPL_ORDER CBLAS_ORDER
extern "C"
{
	int max_gpu_nb = 1024;
	extern int global_m_remain;
	extern int HPL_CALDGEMM_gpu_height;
	
#ifdef HPL_CPUFREQ
	extern int curcpufreq;
	void setcpufreq(int freq, int dgemmfreq);
#endif
}
#include <pthread.h>
#include <errno.h>
#define enum
#include "util_cal.h"
#undef enum

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

	cal_info.linpack_factorize_function = HPL_CALDGEMM_wrapper_factorize;
#ifdef HPL_EMULATE_MULTINODE
	cal_info.linpack_broadcast_function = multinode_broadcast_fake;
	broadcast_fake_array = new float[1024 * 1024 * 2];
#elif !defined(HPL_NO_MPI_LIB)
	cal_info.linpack_broadcast_function = HPL_CALDGEMM_wrapper_broadcast;
#endif

#ifdef HPL_PRINT_THROTTLING_NODES
	cal_info.GPUClock = HPL_PRINT_THROTTLING_NODES;
#endif

#ifdef HPL_GPU_EXTRA_CALDGEMM_OPTIONS
	HPL_GPU_EXTRA_CALDGEMM_OPTIONS
#endif

	cal_dgemm = new CALDGEMM_IMPL;
	if (cal_dgemm == NULL) return(1);
	cal_info.config_backend = cal_dgemm->create_caldgemm_config_backend();
#ifdef HPL_GPU_EXTRA_CALDGEMM_BACKEND_OPTIONS
	HPL_GPU_EXTRA_CALDGEMM_BACKEND_OPTIONS
#endif

#ifdef HPL_CALDGEMM_PARAM
	if (cal_dgemm->ParseParameters(str(HPL_CALDGEMM_PARAM))) return(1);
#endif

	int retVal = cal_dgemm->InitCALDGEMM(&cal_info);

	return(retVal);
}

void CALDGEMM_Shutdown()
{
	cal_dgemm->ExitCALDGEMM();
	delete cal_dgemm;
}

void CALDGEMM_async_dtrsm(const HPL_ORDER ORDER, const HPL_SIDE SIDE, const HPL_UPLO UPLO, const HPL_TRANS TRANS, const HPL_DIAG DIAG, const int M, const int N,
   const double ALPHA, const double *A, const int LDA, double *B, const int LDB)
{
#ifdef HPL_CALDGEMM_ASYNC_FACT_DTRSM
	if (global_m_remain <= HPL_CALDGEMM_ASYNC_FACT_DTRSM)
	{
		if (cal_dgemm->RunAsyncSingleTileDTRSM(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB))
		{
			printf("Error in async CALDGEMM DTRSM Run, aborting HPL Run\n");
			exit(1);
		}
	}
	else
#endif
	{
		cblas_dtrsm(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB);
	}
}

void CALDGEMM_async_dtrsm2(const HPL_ORDER ORDER, const HPL_SIDE SIDE, const HPL_UPLO UPLO, const HPL_TRANS TRANS, const HPL_DIAG DIAG, const int M, const int N,
   const double ALPHA, const double *A, const int LDA, double *B, const int LDB)
{
#ifdef HPL_CALDGEMM_ASYNC_DTRSM
	if (global_m_remain <= HPL_CALDGEMM_ASYNC_DTRSM)
	{
		if (cal_dgemm->RunAsyncSingleTileDTRSM(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB))
		{
			printf("Error in async CALDGEMM DTRSM Run, aborting HPL Run\n");
			exit(1);
		}
	}
	else
#endif
	{
		cblas_dtrsm(ORDER, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B, LDB);
	}
}


void CALDGEMM_async_dgemm( const HPL_ORDER ORDER, const HPL_TRANS TRANSA,
	const HPL_TRANS TRANSB, const int M, const int N,
	const int K, const double ALPHA, const double * A, const int LDA,
	const double * B, const int LDB, const double BETA, double * C,
	const int LDC)
{
#ifdef HPL_CALDGEMM_ASYNC_FACT_DGEMM
	if (global_m_remain <= HPL_CALDGEMM_ASYNC_FACT_DGEMM)
	{
		if (cal_dgemm->RunAsyncSingleTileDGEMM( (double*) A, (double*) B, C, (double) ALPHA, (double) BETA, (int) M, (int) K, (int) N, (int) LDA, (int) LDB, (int) LDC, ORDER == CblasColMajor, TRANSA == CblasTrans, TRANSB == CblasTrans))
		{
			printf("Error in async CALDGEMM Run, aborting HPL Run\n");
			exit(1);
		}
	}
	else
#endif
	{
		cblas_dgemm( ORDER, TRANSA, TRANSB, M, N, K, ALPHA, (double*) A, LDA, (double*) B, LDB, BETA, C, LDC );
	}
}

void CALDGEMM_dgemm( const HPL_ORDER ORDER, const HPL_TRANS TRANSA,
	const HPL_TRANS TRANSB, const int M, const int N,
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

#ifdef HPL_CUSTOM_PARAMETER_CHANGE_CALDGEMM
		HPL_CUSTOM_PARAMETER_CHANGE_CALDGEMM
#endif

		if (cal_dgemm->RunCALDGEMM( (double*) A, (double*) B, C, (double) ALPHA, (double) BETA, (int) M, (int) K, (int) N, (int) LDA, (int) LDB, (int) LDC, ORDER == CblasColMajor, TRANSA == CblasTrans, TRANSB == CblasTrans, LinpackCallbacks ))
		{
			printf("Error in CALDGEMM Run, aborting HPL Run\n");
			exit(1);
		}
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

#ifdef HPL_REGISTER_MEMORY
	bool gpu_access = true;
#else
	bool gpu_access = false;
#endif
	return((void*) cal_dgemm->AllocMemory(size / sizeof(double), page_locked, huge_tables, gpu_access, interleaved));
}

void CALDGEMM_free(void* ptr)
{
#ifdef HPL_REGISTER_MEMORY
	bool gpu_access = true;
#else
	bool gpu_access = false;
#endif
	cal_dgemm->FreeMemory((double*) ptr, gpu_access);
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
