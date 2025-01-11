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


#include "HaplotypeVcfParser.h"

Matrix::Matrix(size_t nrow, size_t mcol)
    : nrow_(nrow), mcol_(mcol), 
    data_(nrow_ > 0 &&  mcol_ > 0 ? std::make_unique<double[]>(size()) : nullptr) {
    
        if (nrow_ <= 0 || mcol_ <= 0)
            throw std::runtime_error("Matrix must have minimum size of 1");

        // set default values to zero
        for (int i = 0; i < size(); i++)
            data_[i] = 0;
    };


double Matrix::operator()(size_t i, size_t j) const {
    return data_[mat_idx_to_array_(i, j)]
}

double& Matrix::operator()(size_t i, size_t j) {
    return data_[mat_idx_to_array_(i, j)]
}

std::array<size_t,2> Matrix::dims() {
    return {nrow_, mcol_};
}


size_t Matrix::mat_idx_to_array_(size_t i, size_t j) {
    if (i >= nrow_ || j >= mcol_)
        throw std::runtime_error("Indices must be postive integers or zero.");

    array_idx = i*mcol_ + j;

    return array_idx
}


size_t Matrix::size() const { return nrow_ * mcol_ };
