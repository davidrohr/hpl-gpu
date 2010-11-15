/*
    Copyright (C) 2010 Matthias Kretz <kretz@kde.org>
    Copyright (C) 2010 Frankfurt Institute for Advanced Studies (FIAS)

    The source code is property of the Frankfurt Institute for Advanced Studies
    (FIAS). None of the material may be copied, reproduced, distributed,
    republished, downloaded, displayed, posted or transmitted in any form or by
    any means, including, but not limited to, electronic, mechanical,
    photocopying, recording, or otherwise, without the prior written permission
    of FIAS.

*/

#ifndef USE_ORIGINAL_LASWP

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <cstdlib>
#include <cstdio>
#include <sched.h>

#include "glibc_hacks.h"

#define USE_DIES 4
#define USE_CORES 3
#define CORE_COUNT 6

#define RESTRICT_CORES

namespace
{
    class HPL_init_laswp_foo
    {
        public:
            void operator()(const tbb::blocked_range<size_t> &) const {}
    };

    int HPL_init_laswp()
    {
        const char *num_threads_string = getenv("LASWP_NUM_THREADS");
        int num_threads = tbb::task_scheduler_init::automatic;
        
#ifndef RESTRICT_CORES
        if (num_threads_string) {
            num_threads = atoi(num_threads_string);
            //fprintf(stderr, "TBB initialized with %d threads\n", num_threads);
        } else {
            //fprintf(stderr, "TBB initialized: %d threads available\n", tbb::task_scheduler_init::default_num_threads());
        }
#else
        num_threads = USE_DIES * USE_CORES;
#endif

        // before we init TBB we must make sure that we are not pinned to a CPU
        cpu_set_t oldmask;
        sched_getaffinity(0, sizeof(cpu_set_t), &oldmask);
        cpu_set_t fullMask;
        CPU_ZERO(&fullMask);

#ifndef RESTRICT_CORES
        for (int i = 0; i < 64; ++i) {
            CPU_SET(i, &fullMask);
        }
        if (CPU_COUNT(&oldmask) == 1
                && std::getenv("LASWP_PIN_WORKER_THREADS")
                && std::strcmp("1", std::getenv("LASWP_PIN_WORKER_THREADS")) == 0) {
            CPU_XOR(&fullMask, &fullMask, &oldmask);
        }
#else

#include "caldgemm/caldgemm_config.h"
#ifdef HPL_NO_HACKED_LIB
	for (int i = CALDGEMM_OUTPUT_THREADS_SLOW + 2;i < 64;i+=2)
#else
	for (int i = CALDGEMM_OUTPUT_THREADS + 2;i < 64;i+=2)
#endif
	{
		CPU_SET(i, &fullmask);
	}

/*        for (int i = 0;i < USE_DIES;i++)
        {
    	    for (int j = 0;j < USE_CORES;j++)
    	    {
    		CPU_SET(i * CORE_COUNT + CORE_COUNT - j - 1, &fullMask);
    	    }
    	}*/
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

    int _HPL_init_laswp = HPL_init_laswp();
}

#endif
