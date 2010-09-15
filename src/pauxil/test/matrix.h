/*  This file is part of the Vc library.

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
