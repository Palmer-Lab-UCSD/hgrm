
#include "../include/Matrix.h"
#include <gtest/gtest.h>


TEST(TestMatrix, initialize) {
    size_t n_row { 3 };
    size_t m_col { 2 };

    Matrix a { n_row, m_col };

    std::array<size_t, 2> dims { a.dims() };
    EXPECT_EQ(dims[0], n_row);
    EXPECT_EQ(dims[1], m_col);

    for (int i = 0; i < n_row; i++)
        for (int j = 0; j < m_col; j++)
            EXPECT_EQ(a(i, j), 0);

    n_row = 0;

    EXPECT_THROW({
            size_t n_row = 0;
            size_t m_col = 2;
            Matrix b(n_row, m_col);
            },
            std::runtime_error);

    EXPECT_ANY_THROW({
            size_t n_row = -1;
            size_t m_col = 2;
            Matrix b(n_row, m_col);
            });

    EXPECT_THROW({Matrix b(1, 0);}, std::runtime_error);
    EXPECT_ANY_THROW({Matrix b(1, -1);});
}



TEST(TestMatrix, val) {
    size_t n_row { 3 };
    size_t m_col { 5 };

    Matrix a { n_row, m_col };
    
    std::array<size_t, 2> dims { a.dims() };

    double x = { 1 };
    for (int i = 0; i < dims[0]; i++)
        for (int j = 0; j < dims[1]; j++)
            a(i, j) = x++;

    x = 1;
    for (int i = 0; i < dims[0]; i++)
        for (int j = 0; j < dims[1]; j++)
            EXPECT_FLOAT_EQ(a(i, j), x++);
}


TEST(TestMatrix, out_of_bounds) {
    size_t n_row { 3 };
    size_t m_col { 5 };
    
    Matrix a { n_row, m_col };

    EXPECT_THROW({ double tmp = a(4, 3); }, std::runtime_error);
    EXPECT_THROW({ double tmp = a(3, 5); }, std::runtime_error);
    EXPECT_THROW({ double tmp = a(3, 2); }, std::runtime_error);
    EXPECT_THROW({ double tmp = a(2, 5); }, std::runtime_error);
    EXPECT_THROW({ double tmp = a(-2, 4); }, std::runtime_error);

}


TEST(TestMatrix, dim_and_size) {
    size_t n_row { 3 };
    size_t m_col { 5 };
    
    Matrix a { n_row, m_col };
    
    std::array<size_t, 2> dims { a.dims() };
    EXPECT_EQ(dims[0], n_row);
    EXPECT_EQ(dims[1], m_col);

    EXPECT_EQ(a.size(), n_row * m_col);
}
