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
#ifndef HEADER_GRM_H
#define HEADER_GRM_H

#include <cstdio>
#include <cstddef>
#include <stdexcept>
#include <memory>
#include <array>
#include <utility>


class Grm {
public:
    Grm(const size_t, const size_t);
    Grm(const Grm&);                          // copy constructor
    Grm(Grm&&);                               // move constructor
    Grm& operator=(const Grm&)=delete;        // copy assignment
    Grm& operator=(Grm&&)=delete;             // move assignment
                                            

    double operator()(const size_t&, const size_t&) const;
    double& operator()(const size_t&, const size_t&);

    size_t size() const;
    std::array<size_t,2> dims() const;

    int write(const char *filename) const;
    int write(const std::string& filename) const;
    static int read(const char *filename, Grm *grm);
    static int read(const std::string& filename, Grm *grm);

private:
    const size_t nrow_;
    const size_t mcol_;
    std::unique_ptr<double[]> data_;
    size_t mat_idx_to_array_(const size_t&, const size_t&) const;
};


#endif
