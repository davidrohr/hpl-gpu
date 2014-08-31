#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

extern "C" {
#include "hpl_auxil.h"
#include "hpl_pauxil.h"
}
#include "../../caldgemm/cmodules/timer.h"
#include "../../caldgemm/cmodules/qmalloc.h"

#include "../../testing/ptest/fastmatgen.h"

#define quit(...) {fprintf(stderr, __VA_ARGS__); exit(1);}

int main(int argc, char** argv)
{
	HighResTimer timer;
	FILE* fp;
	//dlacpy.24960.01920.28872.24960. 0.0808s.dat  dlaswp00N.1920.87169.89160. 0.3553s.dat  dlaswp01T.2496.23809.29640.23816. 0.4362s.dat  dlaswp06T.0590.23809.29640.23816. 0.1770s.dat  dlaswp10N.23040.01920.23048. 0.0441s.dat  dlatcpy.1920.23809.23816.29640. 0.0516s.dat
	
	fprintf(stderr, "Allocating swap vector memory\n");
	
	int benchmark_warmup = 10;
	int benchmark_iterations = 50;
	
	int dlacpy_M = 24960;
	int dlacpy_N = 1920;
	int dlacpy_LDA = 28872;
	int dlacpy_LDB = 24960;
	
	int dlatcpy_M = 1920;
	int dlatcpy_N = 23809;
	int dlatcpy_LDA = 23816;
	int dlatcpy_LDB = 29640;
	
	int dlaswp00N_M = 1920;
	int dlaswp00N_N = 87169;
	int dlaswp00N_LDA = 89160;
	int* dlaswp00N_IPIV = new int[dlaswp00N_M];
	const char* dlaswp00N_FILE = "dlaswp00N.1920.87169.89160. 0.3553s.dat";
	
	int dlaswp01T_M = 2496;
	int dlaswp01T_N = 23809;
	int dlaswp01T_LDA = 29640;
	int dlaswp01T_LDU = 23816;
	int* dlaswp01T_LINDXA = new int[dlaswp01T_M];
	int* dlaswp01T_LINDXAU = new int[dlaswp01T_M];
	const char* dlaswp01T_FILE = "dlaswp01T.2496.23809.29640.23816. 0.4362s.dat";
	
	int dlaswp06T_M = 590;
	int dlaswp06T_N = 23809;
	int dlaswp06T_LDA = 29640;
	int dlaswp06T_LDU = 23816;
	int* dlaswp06T_LINDXA = new int[dlaswp06T_M];
	const char* dlaswp06T_FILE = "dlaswp06T.0590.23809.29640.23816. 0.1770s.dat";
	
	int dlaswp10N_M = 23040;
	int dlaswp10N_N = 1920;
	int dlaswp10N_LDA = 23048;
	int* dlaswp10N_IPIV = new int[dlaswp10N_N];
	const char* dlaswp10N_FILE = "dlaswp10N.23040.01920.23048. 0.0441s.dat";
	
	if (dlaswp00N_IPIV == NULL || dlaswp01T_LINDXA == NULL || dlaswp01T_LINDXAU == NULL || dlaswp06T_LINDXA == NULL || dlaswp10N_IPIV == NULL) quit("Memmory allocation error (swap vectors)\n");
	
	fprintf(stderr, "Opening input files, reading swap vectors\n");

	if ((fp = fopen(dlaswp00N_FILE, "rb")) == NULL) quit("Error opening file %s\n", dlaswp00N_FILE);
	fread(dlaswp00N_IPIV, sizeof(dlaswp00N_IPIV[0]), dlaswp00N_M, fp);
	fclose(fp);

	if ((fp = fopen(dlaswp01T_FILE, "rb")) == NULL) quit("Error opening file %s\n", dlaswp01T_FILE);
	fread(dlaswp01T_LINDXA, sizeof(dlaswp01T_LINDXA[0]), dlaswp01T_M, fp);
	fread(dlaswp01T_LINDXAU, sizeof(dlaswp01T_LINDXAU[0]), dlaswp01T_M, fp);
	fclose(fp);

	if ((fp = fopen(dlaswp06T_FILE, "rb")) == NULL) quit("Error opening file %s\n", dlaswp06T_FILE);
	fread(dlaswp06T_LINDXA, sizeof(dlaswp06T_LINDXA[0]), dlaswp06T_M, fp);
	fclose(fp);

	if ((fp = fopen(dlaswp10N_FILE, "rb")) == NULL) quit("Error opening file %s\n", dlaswp10N_FILE);
	fread(dlaswp10N_IPIV, sizeof(dlaswp10N_IPIV[0]), dlaswp10N_N, fp);
	fclose(fp);

	int max_size = std::max({dlacpy_M, dlacpy_N, dlacpy_LDA, dlacpy_LDB, dlatcpy_M, dlatcpy_N, dlatcpy_LDA, dlatcpy_LDB, dlaswp00N_M, dlaswp00N_N, dlaswp00N_LDA, dlaswp01T_M, dlaswp01T_N, dlaswp01T_LDA, dlaswp01T_LDU, dlaswp06T_M, dlaswp06T_N, dlaswp06T_LDA, dlaswp06T_LDU, dlaswp10N_M, dlaswp10N_N, dlaswp10N_LDA});
	size_t matrix_size = (size_t) max_size * (size_t) max_size * sizeof(double);
	
	fprintf(stderr, "Allocating matrix memory (size %lld KB)", (long long int) matrix_size / 1024);
	
	double* matrix_1 = (double*) qmalloc::qMalloc(matrix_size, false, false, true, NULL, true);
	fprintf(stderr, ".");
	double* matrix_2 = (double*) qmalloc::qMalloc(matrix_size, false, false, true, NULL, true);
	fprintf(stderr, ".\n");
	if (matrix_1 == NULL || matrix_2 == NULL) quit("Memory allocation error (matrices)\n");
	
	fprintf(stderr, "Filling matrices with random data\n");
	fastmatgen(1, matrix_1, matrix_size / sizeof(double));
	fastmatgen(2, matrix_2, matrix_size / sizeof(double));
	
	fprintf(stderr, "Initialization done, running benchmarks\n");
	
	for (int i = 0;i < benchmark_warmup + benchmark_iterations;i++)
	{
		HPL_dlacpy(dlacpy_M, dlacpy_N, matrix_1, dlacpy_LDA, matrix_2, dlacpy_LDB, 1);
		if (i == benchmark_warmup) timer.ResetStart();
	}
	fprintf(stderr, "Time dlacpy: %fs\n", timer.GetCurrentElapsedTime() / (double) benchmark_iterations);

	for (int i = 0;i < benchmark_warmup + benchmark_iterations;i++)
	{
		HPL_dlatcpy(dlatcpy_M, dlatcpy_N, matrix_1, dlatcpy_LDA, matrix_2, dlatcpy_LDB);
		if (i == benchmark_warmup) timer.ResetStart();
	}
	fprintf(stderr, "Time dlatcpy: %fs\n", timer.GetCurrentElapsedTime() / (double) benchmark_iterations);

	for (int i = 0;i < benchmark_warmup + benchmark_iterations;i++)
	{
		HPL_dlaswp00N(dlaswp00N_M, dlaswp00N_N, matrix_1, dlaswp00N_LDA, dlaswp00N_IPIV);
		if (i == benchmark_warmup) timer.ResetStart();
	}
	fprintf(stderr, "Time dlaswp00N: %fs\n", timer.GetCurrentElapsedTime() / (double) benchmark_iterations);

	for (int i = 0;i < benchmark_warmup + benchmark_iterations;i++)
	{
		HPL_dlaswp01T(dlaswp01T_M, dlaswp01T_N, matrix_1, dlaswp01T_LDA, matrix_2, dlaswp01T_LDU, dlaswp01T_LINDXA, dlaswp01T_LINDXAU);
		if (i == benchmark_warmup) timer.ResetStart();
	}
	fprintf(stderr, "Time dlaswp01T: %fs\n", timer.GetCurrentElapsedTime() / (double) benchmark_iterations);

	for (int i = 0;i < benchmark_warmup + benchmark_iterations;i++)
	{
		HPL_dlaswp06T(dlaswp06T_M, dlaswp06T_N, matrix_1, dlaswp06T_LDA, matrix_2, dlaswp06T_LDU, dlaswp06T_LINDXA);
		if (i == benchmark_warmup) timer.ResetStart();
	}
	fprintf(stderr, "Time dlaswp06T: %fs\n", timer.GetCurrentElapsedTime() / (double) benchmark_iterations);

	for (int i = 0;i < benchmark_warmup + benchmark_iterations;i++)
	{
		HPL_dlaswp10N(dlaswp10N_M, dlaswp10N_N, matrix_1, dlaswp10N_LDA, dlaswp10N_IPIV);
		if (i == benchmark_warmup) timer.ResetStart();
	}
	fprintf(stderr, "Time dlaswp10N: %fs\n", timer.GetCurrentElapsedTime() / (double) benchmark_iterations);
	
	fprintf(stderr, "Benchmarks done, cleaning up\n");
	
	delete[] dlaswp00N_IPIV;
	delete[] dlaswp01T_LINDXA;
	delete[] dlaswp01T_LINDXAU;
	delete[] dlaswp06T_LINDXA;
	delete[] dlaswp10N_IPIV;
	
	qmalloc::qFree(matrix_1);
	qmalloc::qFree(matrix_2);
	
	fprintf(stderr, "All done, exiting\n");
	
	return(0);
}
