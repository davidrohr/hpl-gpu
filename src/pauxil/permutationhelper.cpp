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

#ifndef HPL_USE_ORIGINAL_LASWP

#include "permutationhelper.h"

#include <cstdlib>
#include <cstdio>

PermutationHelper *__restrict__ PermutationHelper::s_instance = 0;

PermutationHelper::~PermutationHelper()
{
    free(m_data);
}

void PermutationHelper::resize(size_t size)
{
    m_size = size;
    free(m_data);
    void *tmp;
    int err = posix_memalign(&tmp, 4096, m_size * sizeof(Data));
    m_data = static_cast<Data *>(tmp);
    if (0 != err) {
        fprintf(stderr, "posix_memalign failed to allocate %ld Bytes\n", m_size * sizeof(Data));
        abort();
    }
}

int PermutationHelper::_init()
{
    static PermutationHelper tmp;
    s_instance = &tmp;
    s_instance->resize(1024);

    return 1;
}

namespace
{
    int init = PermutationHelper::_init();
} // namespace 

#endif

// vim: sw=4 sts=4 et tw=100
