/*
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

#ifndef USE_ORIGINAL_LASWP

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <cstdlib>
#include <cstdio>
#include <sched.h>

#include "glibc_hacks.h"

#ifdef HPL_CALL_CALDGEMM
#include "../../caldgemm/cmodules/affinity.h"
#endif

#ifdef HPL_CALL_CALDGEMM
#define USE_DIES 4
#define USE_CORES 3
#define CORE_COUNT 6

extern "C" int get_num_procs();

#ifndef HPL_USE_ALL_CORES_FOR_LASWP
#define RESTRICT_CORES
#endif
#include <unistd.h>

#include "../../caldgemm/caldgemm.h"
#else

static inline int get_num_procs()
{
    return(sysconf(_SC_NPROCESSORS_ONLN));
}

#endif


namespace
{
    extern "C" int HPL_init_laswp(void* ptr);

    class HPL_init_laswp_foo
    {
		public:
		    void operator()(const tbb::blocked_range<size_t> &) const {
#ifdef HPL_CALL_CALDGEMM
			if (gettid() != getpid()) setThreadName("LASWP");
#endif
		    }
    };


    int HPL_init_laswp(void* ptr)
    {
		int num_threads = tbb::task_scheduler_init::automatic;
		
		// before we init TBB we must make sure that we are not pinned to a CPU
		cpu_set_t oldmask;
		sched_getaffinity(0, sizeof(cpu_set_t), &oldmask);
		cpu_set_t fullMask;
		CPU_ZERO(&fullMask);

#ifndef RESTRICT_CORES
		/*const char *num_threads_string = getenv("LASWP_NUM_THREADS");
		if (num_threads_string) {
		    num_threads = atoi(num_threads_string);
		    //fprintf(stderr, "TBB initialized with %d threads\n", num_threads);
		} else {
		    //fprintf(stderr, "TBB initialized: %d threads available\n", tbb::task_scheduler_init::default_num_threads());
		}*/
		for (int i = 0; i < get_num_procs(); ++i) {
		    CPU_SET(i, &fullMask);
		}
		if (CPU_COUNT(&oldmask) == 1
				&& std::getenv("LASWP_PIN_WORKER_THREADS")
				&& std::strcmp("1", std::getenv("LASWP_PIN_WORKER_THREADS")) == 0) {
		    CPU_XOR(&fullMask, &fullMask, &oldmask);
		}
		num_threads = get_num_procs();
#else
		caldgemm* cal_dgemm = (caldgemm*) ptr;
		int num_procs = get_num_procs();
		for (int i = 0;i < num_procs;i++)
		{
			if (cal_dgemm->cpuUsed(i) == false) CPU_SET(i, &fullMask);
		}
		num_threads = CPU_COUNT(&fullMask);
		printf("Using %d threads for LASWP ( ", num_threads);
		for (int i = 0;i < num_procs;i++)
		{
			if (CPU_ISSET(i, &fullMask)) printf("%d ", i);
		}
		printf(")\n");

#endif
		sched_setaffinity(0, sizeof(cpu_set_t), &fullMask);
		sched_getaffinity(0, sizeof(cpu_set_t), &fullMask);

		//fprintf(stderr, "Pin TBB worker threads to core(s) 0x%016lX\n", fullMask.__bits[0]);
		static tbb::task_scheduler_init init(num_threads);
		tbb::parallel_for (tbb::blocked_range<size_t>(0, 100), HPL_init_laswp_foo());

		//fprintf(stderr, "       Pin main thread to core(s) 0x%016lX\n", oldmask.__bits[0]);
		sched_setaffinity(0, sizeof(cpu_set_t), &oldmask);

		return 0;
    }
}

#endif
