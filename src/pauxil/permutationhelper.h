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
