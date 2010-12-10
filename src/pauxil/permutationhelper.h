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

#ifndef PERMUTATIONHELPER_H
#define PERMUTATIONHELPER_H

#include <cstddef>

class PermutationHelper
{
    public:
        static PermutationHelper &instance() { return *s_instance; }

        void resize(size_t size);
        inline void ensureSize(size_t size) {
            if (__builtin_expect(size > m_size, 0)) {
                resize(size);
            }
        }

        struct Data
        {
            int a, b;
        };

        Data &__restrict__ operator[](int i) { return m_data[i]; }
        const Data &__restrict__ operator[](int i) const { return m_data[i]; }

        //private:
        static int _init();

    private:
        inline PermutationHelper() : m_data(0) {}
        ~PermutationHelper();

        Data *__restrict__ m_data;
        size_t m_size;

        static PermutationHelper *__restrict__ s_instance;
};

#endif // PERMUTATIONHELPER_H
