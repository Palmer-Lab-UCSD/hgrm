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

#include <grm.h>

// default constructor
Grm::Grm(const size_t nrow, const size_t mcol)
    : nrow_(nrow), mcol_(mcol), 
    data_(nrow_ > 0 &&  mcol_ > 0 ? std::make_unique<double[]>(size()) : nullptr) {
    
        if (nrow_ == 0 || mcol_ == 0)
            throw std::runtime_error("Grm must have minimum size of 1");

        // set default values to zero
        for (size_t i = 0; i < size(); i++)
            data_[i] = 0;
    };


// copy constructor
//
Grm::Grm(const Grm& other) 
    : nrow_(other.nrow_), mcol_(other.mcol_),
    data_(std::make_unique<double[]>(other.size())) {

        // Grm values have already been validated
        for (size_t i = 0; i < size(); i++)
            data_[i] = other.data_[i];
}

// TODO: check this.
Grm::Grm(Grm&& other) 
    : nrow_(other.nrow_), mcol_(other.mcol_), data_(std::move(other.data_)) {};


double Grm::operator()(const size_t& i, const size_t& j) const {
    return data_[mat_idx_to_array_(i, j)];
}

double& Grm::operator()(const size_t& i, const size_t& j) {
    return data_[mat_idx_to_array_(i, j)];
}

std::array<size_t,2> Grm::dims() const {
    return {nrow_, mcol_};
}


size_t Grm::mat_idx_to_array_(const size_t& i, const size_t& j) const {
    if (i >= nrow_ || j >= mcol_)
        throw std::runtime_error("Indices must be postive integers or zero.");

    return i*mcol_ + j;
}


size_t Grm::size() const { return nrow_ * mcol_; };

int Grm::write(const std::string& filename) const {
    return -1;
}

int Grm::write(const char *filename) const {
    return -1;
}

int Grm::read(const char *filename, Grm *grm) {
    return -1;
}
