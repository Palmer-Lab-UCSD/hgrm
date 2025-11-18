//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
//
// Acknowledgment
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Sonnet, the AI assistant from Anthropic
// (Jan 2025), with minor recommendations incorporated.
//
//
#ifndef HEADER_MATRIX_H
#define HEADER_MATRIX_H

#include <cstddef>
#include <stdexcept>
#include <memory>
#include <array>
#include <utility>

class Matrix
{
public:
    Matrix(size_t, size_t);                         // constructorconstructor
    Matrix(const Matrix&);                          // copy constructor
    Matrix(Matrix&&);                               // move constructor
    Matrix& operator=(const Matrix&)=delete;        // copy assignment
    Matrix& operator=(Matrix&&)=delete;             // move assignment
                                            

    double operator()(const size_t&, const size_t&) const;
    double& operator()(const size_t&, const size_t&);

    size_t size() const;
    std::array<size_t,2> dims() const;

private:
    const size_t nrow_;
    const size_t mcol_;
    std::unique_ptr<double[]> data_;
    size_t mat_idx_to_array_(const size_t&, const size_t&) const;
};

#endif
