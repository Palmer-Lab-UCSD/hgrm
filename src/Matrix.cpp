// MAtrix
//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-10
//
//
// Acknowledgment
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Sonnet, the AI assistant from Anthropic
// (Jan 2025), with minor recommendations incorporated.
//
//

#include "Matrix.h"

// default constructor
Matrix::Matrix(size_t nrow, size_t mcol)
    : nrow_(nrow), mcol_(mcol), 
    data_(nrow_ > 0 &&  mcol_ > 0 ? std::make_unique<double[]>(size()) : nullptr) {
    
        if (nrow_ <= 0 || mcol_ <= 0)
            throw std::runtime_error("Matrix must have minimum size of 1");

        // set default values to zero
        for (size_t i = 0; i < size(); i++)
            data_[i] = 0;
    };


// copy constructor
//
Matrix::Matrix(const Matrix& other) 
    : nrow_(other.nrow_), mcol_(other.mcol_),
    data_(std::make_unique<double[]>(other.size())) {

        // Matrix values have already been validated
        for (size_t i = 0; i < size(); i++)
            data_[i] = other.data_[i];
}


Matrix::Matrix(Matrix&& other) 
    : nrow_(other.nrow_), mcol_(other.mcol_), data_(std::move(other.data_)) {};


double Matrix::operator()(const size_t& i, const size_t& j) const {
    return data_[mat_idx_to_array_(i, j)];
}

double& Matrix::operator()(const size_t& i, const size_t& j) {
    return data_[mat_idx_to_array_(i, j)];
}

std::array<size_t,2> Matrix::dims() const {
    return {nrow_, mcol_};
}


size_t Matrix::mat_idx_to_array_(const size_t& i, const size_t& j) const {
    if (i >= nrow_ || j >= mcol_)
        throw std::runtime_error("Indices must be postive integers or zero.");

    return i*mcol_ + j;
}


size_t Matrix::size() const { return nrow_ * mcol_; };
