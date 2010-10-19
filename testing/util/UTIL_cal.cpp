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

#include "util_cal.h"

extern "C"
{
    extern void HPL_CALLDGEMM_wrapper_factorize();
    extern void HPL_CALLDGEMM_wrapper_broadcast();
}

static caldgemm::SampleInfo cal_info;
static caldgemm cal_dgemm;

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
	cal_info.NoPerformanceWarnings = CAL_FALSE;
	cal_info.KeepBuffersMapped = CAL_TRUE;
	
	cal_info.linpack_factorize_function = HPL_CALLDGEMM_wrapper_factorize;
	cal_info.linpack_broadcast_function = HPL_CALLDGEMM_wrapper_broadcast;

	return(cal_dgemm.InitCALDGEMM( &cal_info ));
}

void CALDGEMM_Shutdown()
{
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
	    if (cal_dgemm.RunCALDGEMM( (double*) A, (double*) B, C, (double) ALPHA, (double) BETA, (int) M, (int) K, (int) N, (int) LDA, (int) LDB, (int) LDC, ORDER, TRANSA, TRANSB, LinpackCallbacks ))
	    {
		printf("Error in CALDGEMM Run, aborting HPL Run\n");
		exit(1);
	    }
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
