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

#include <cstdlib>
#include <cstring>
#include <signal.h>
#include <cstdio>

class Matrix
{
    public:
        Matrix(size_t rows, size_t cols);
        Matrix(const Matrix &rhs);
        ~Matrix();

        bool operator==(const Matrix &rhs) const;

        size_t LDA() const { return m_rows; }

        double *at(size_t row, size_t col) {
            return &m_data[row + m_rows * col];
        }

    private:
        Matrix &operator=(const Matrix &);

        double *m_data;
        size_t m_rows, m_cols;
};

inline Matrix::Matrix(size_t rows, size_t cols)
    : m_rows(rows), m_cols(cols)
{
    void *tmp;
    posix_memalign(&tmp, 64, m_rows * m_cols * sizeof(double));
    m_data = static_cast<double *>(tmp);
    for (size_t i = 0; i < m_rows * m_cols; ++i) {
        m_data[i] = i;
    }
}

inline Matrix::Matrix(const Matrix &rhs)
    : m_rows(rhs.m_rows), m_cols(rhs.m_cols)
{
    void *tmp;
    posix_memalign(&tmp, 64, m_rows * m_cols * sizeof(double));
    m_data = static_cast<double *>(tmp);
    std::memcpy(m_data, rhs.m_data, m_rows * m_cols * sizeof(double));
}

inline Matrix::~Matrix()
{
    free(m_data);
}

inline bool Matrix::operator==(const Matrix &rhs) const
{
    return 0 == std::memcmp(m_data, rhs.m_data, m_rows * m_cols * sizeof(double));
}

#define COMPARE(a, b) \
    if ((a) == (b)) {} else { \
        fprintf(stderr, "%s:%d: %s != %s\n", __FILE__, __LINE__, #a, #b); \
        raise(SIGTRAP); \
    }
