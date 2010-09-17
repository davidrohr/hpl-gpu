/*
    Copyright (C) 2010 Matthias Kretz <kretz@kde.org>

    Vc is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of
    the License, or (at your option) any later version.

    Vc is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Vc.  If not, see <http://www.gnu.org/licenses/>.

*/

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

// vim: sw=4 sts=4 et tw=100
