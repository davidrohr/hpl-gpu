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

#include <tbb/task_scheduler_init.h>
#include <cstdlib>
#include <cstdio>

namespace
{
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
        static tbb::task_scheduler_init init(num_threads);
        return 0;
    }

    int _HPL_init_laswp = HPL_init_laswp();
}
