
#include <gtest/gtest.h>
#include <grm.h>
#include <cstddef>


TEST(TestGrm, Init) {
    size_t n_row { 3 };
    size_t m_col { 2 };
    Grm a { n_row, m_col };
    std::array<size_t, 2> dims { a.dims() };
    EXPECT_EQ(dims[0], n_row);
    EXPECT_EQ(dims[1], m_col);

    for (size_t i = 0; i < n_row; i++)
        for (size_t j = 0; j < m_col; j++)
            EXPECT_EQ(a(i, j), 0);

    n_row = 0;

    EXPECT_THROW({
            size_t n_row = 0;
            size_t m_col = 2;
            Grm b(n_row, m_col);
            },
            std::runtime_error);

    EXPECT_ANY_THROW({
            size_t n_row = -1;
            size_t m_col = 2;
            Grm b(n_row, m_col);
            });

    EXPECT_THROW({Grm b(1, 0);}, std::runtime_error);
    EXPECT_ANY_THROW({Grm b(1, -1);});
}



TEST(TestGrm, Vals) {
    size_t n_row { 3 };
    size_t m_col { 5 };

    Grm a { n_row, m_col };
    
    std::array<size_t, 2> dims { a.dims() };

    double x = { 1 };
    for (size_t i = 0; i < dims[0]; i++)
        for (size_t j = 0; j < dims[1]; j++)
            a(i, j) = x++;

    x = 1;
    for (size_t i = 0; i < dims[0]; i++)
        for (size_t j = 0; j < dims[1]; j++)
            EXPECT_FLOAT_EQ(a(i, j), x++);
}


TEST(TestGrm, OutOfBounds) {
    size_t n_row { 3 };
    size_t m_col { 5 };
    
    Grm a { n_row, m_col };

    EXPECT_THROW({ a(4, 3); }, std::runtime_error);
    EXPECT_THROW({ a(3, 5); }, std::runtime_error);
    EXPECT_THROW({ a(3, 2); }, std::runtime_error);
    EXPECT_THROW({ a(2, 5); }, std::runtime_error);
    EXPECT_THROW({ a(-2, 4); }, std::runtime_error);

}


TEST(TestGrm, DimAndSize) {
    size_t n_row { 3 };
    size_t m_col { 5 };
    
    Grm a { n_row, m_col };
    
    std::array<size_t, 2> dims { a.dims() };
    EXPECT_EQ(dims[0], n_row);
    EXPECT_EQ(dims[1], m_col);

    EXPECT_EQ(a.size(), n_row * m_col);
}
