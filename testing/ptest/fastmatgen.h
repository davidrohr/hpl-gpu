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

#include <pthread.h>
#include <unistd.h>

#define FASTRAND_THREADS nfastmatthreads
#define FASTRAND_THREADS_MAX 64
int fastrand_seed;
volatile int* fastrand_done;
double* fastrand_A;
size_t fastrand_size;

static int nfastmatthreads;

void* fastmatgen_slave(void* arg)
{
	int num = (int) (size_t) arg;
	
	cpu_set_t mask;
	CPU_ZERO(&mask);
	CPU_SET(num, &mask);
	sched_setaffinity(0, sizeof(cpu_set_t), &mask);

	size_t fastrand_num = fastrand_seed + 65537 * num;
	const size_t fastrand_mul = 84937482743;
	const size_t fastrand_add = 138493846343;
	const size_t fastrand_mod = 538948374763;

	size_t sizeperthread = fastrand_size / FASTRAND_THREADS;

	double* A = fastrand_A + num * sizeperthread;
	size_t size = (num == FASTRAND_THREADS - 1) ? (fastrand_size - (FASTRAND_THREADS - 1) * sizeperthread) : sizeperthread;

	for (size_t i = 0;i < size;i++)
	{
		fastrand_num = (fastrand_num * fastrand_mul + fastrand_add) % fastrand_mod;
		A[i] = (double) -0.5 + (double) fastrand_num / (double)fastrand_mod;
	}

	fastrand_done[num] = 1;
	return(NULL);
}

void fastmatgen(int SEED, double* A, size_t size)
{
	fastrand_seed = SEED;
	fastrand_A = A;
	fastrand_size = size;

	nfastmatthreads = sysconf(_SC_NPROCESSORS_CONF);
	if (nfastmatthreads < 1) nfastmatthreads = 1;
	if (nfastmatthreads > FASTRAND_THREADS_MAX) nfastmatthreads = FASTRAND_THREADS_MAX;

	fastrand_done = (volatile int*) malloc(FASTRAND_THREADS * sizeof(int));

	memset((void*) fastrand_done, 0, FASTRAND_THREADS * sizeof(int));

	cpu_set_t oldmask;
	sched_getaffinity(0, sizeof(cpu_set_t), &oldmask);

	for (int i = 0;i < FASTRAND_THREADS - 1;i++)
	{
		pthread_t thr;
		pthread_create(&thr, NULL, fastmatgen_slave, (void*) (size_t) i);
	}
	fastmatgen_slave((void*) (size_t) (FASTRAND_THREADS - 1));

	for (int i = 0;i < FASTRAND_THREADS;i++)
	{
		while (fastrand_done[i] == 0);
	}
	free((void*) fastrand_done);
	sched_setaffinity(0, sizeof(cpu_set_t), &oldmask);
}

