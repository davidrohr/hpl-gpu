/**
 * A utility for timing stuff
 *
 * (c) Matthias Bach <bach@compeng.uni-frankfurt.de> 2010
 */

#include <caldgemm.h>

#include "util_cal.h"

static caldgemm::SampleInfo cal_info;
static caldgemm cal_dgemm;

void CALDGEMM_Init()
{
	cal_info.Pin = -3;
	cal_info.Verify = CAL_FALSE;
	cal_info.PrintIL = CAL_FALSE;
	cal_info.Disassemble = CAL_FALSE;
	cal_info.Quiet = CAL_FALSE;
	cal_info.DeviceNum = 0;
	cal_info.Width = 1024; //k for matrix multiply
	cal_info.Height = 2048;
	cal_info.Iterations = 1;
	cal_info.DstMemory = 'c';
	cal_info.DstCacheable = CAL_TRUE;
	cal_info.VerboseTiming = CAL_FALSE;
	cal_info.Debug = CAL_FALSE;
	cal_info.MultiThread = CAL_TRUE;
	cal_info.UseGPU = CAL_TRUE;
	cal_info.UseCPU = CAL_TRUE;
	cal_info.GPURatio = 0.66;

	cal_dgemm.InitCALDGEMM( &cal_info );
}

void CALDGEMM_Shutdown()
{
	cal_dgemm.ExitCALDGEMM();
}

void CALDGEMM_dgemm( const enum CBLAS_ORDER ORDER, const enum CBLAS_TRANSPOSE TRANSA,
                     const enum CBLAS_TRANSPOSE TRANSB, const int M, const int N,
                     const int K, const double ALPHA, const double * A, const int LDA,
                     const double * B, const int LDB, const double BETA, double * C,
                     const int LDC )
{
	// check wether K matches the restrictions of calblas
	if( K == cal_info.Width && ORDER == CblasColMajor && TRANSA == CblasNoTrans && TRANSB == CblasNoTrans )
	{
		// cal_dgemm assumes row major order but linpack is always column major
		// therefore it needs to work on C^T -> we need to switch A and B
		// TODO insert a check, though redundent
		cal_dgemm.ResetTimers();
		cal_dgemm.RunCALDGEMM( (double*) B, (double*) A, C, (double) ALPHA, (double) BETA, (int) N, (int) M, (int) LDB, (int) LDA, (int) LDC );
	}
	else
	{
		// Use plain cblas
		cblas_dgemm( ORDER, TRANSA, TRANSB, M, N, K, ALPHA, (double*) A, LDA, (double*) B, LDB,
		             BETA, C, LDC );
	}
}

