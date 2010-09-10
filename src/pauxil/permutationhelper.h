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

#ifndef PERMUTATIONHELPER_H
#define PERMUTATIONHELPER_H

#include <cstddef>

class PermutationHelper
{
    public:
        static PermutationHelper &instance() { return *s_instance; }

        void resize(size_t size);
        inline void ensureSize(size_t size) {
            if (size > m_size) {
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
