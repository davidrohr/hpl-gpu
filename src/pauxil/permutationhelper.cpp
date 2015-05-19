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
#pragma GCC diagnostic ignored "-Wunused-variable"
	int init = PermutationHelper::_init();
} // namespace 

#endif

// vim: sw=4 sts=4 et tw=100
