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
        if (num_threads_string) {
            num_threads = atoi(num_threads_string);
            fprintf(stderr, "TBB initialized with %d threads\n", num_threads);
        } else {
            fprintf(stderr, "TBB initialized: %d threads available\n", tbb::task_scheduler_init::default_num_threads());
        }

        // before we init TBB we must make sure that we are not pinned to a CPU
        cpu_set_t oldmask;
        sched_getaffinity(0, sizeof(cpu_set_t), &oldmask);
        cpu_set_t fullMask;
        CPU_ZERO(&fullMask);
        for (int i = 0; i < 24; ++i) {
            CPU_SET(i, &fullMask);
        }
        sched_setaffinity(0, sizeof(cpu_set_t), &fullMask);

        static tbb::task_scheduler_init init(num_threads);
        tbb::parallel_for (tbb::blocked_range<size_t>(0, 100), HPL_init_laswp_foo());

        sched_setaffinity(0, sizeof(cpu_set_t), &oldmask);

        return 0;
    }

    int _HPL_init_laswp = HPL_init_laswp();
}

#endif
