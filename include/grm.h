// Palmer Lab at UCSD
//
// This library provides the data structure of a genetic relationship matrix
// and functions for file I/O.
// 
// GRM BINARY FILE SPECIFICATION
//
// A computed GRM is stored in a custom binary format.  The extension ".grm"
// of these files is mandatory.  The file is divided into two components, a
// header with meta-data necessary to reproduce the grm calculation and the
// the computed grm values, named the payload.  
//
// The .grm file header is defined by the struct Hdr, and contains, at a
// minimum, the following information:
//      * program_version: grm program version number
//      * data_type: alt_count, expected_alt_count, expected_haplotype_count, 
//          both expected_alt_count and expected_haplotype_count.
//      * coords: Genomic coordinates used in the grm calculation.
//      * samples: list of sample id's in order of the grm
// the coords and samples are defined by their own structs with field pointers
// to heap allocated memory addresses.  Reading and writing such heap allocated
// structs make use of runtime polymorphism of function "read" and "write"
//
//
// ACKNOWLEDGMENT
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Sonnet, the AI assistant from Anthropic
// (Jan 2025), with minor recommendations incorporated.
// 
#ifndef HEADER_GRM_H
#define HEADER_GRM_H

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <memory>
#include <utility>
#include <string>

#include "io.h"


// The algorithm for getting the array idx from matrix indexes is simply
//
// i * n_samples - n_skipped_idxs + j
//
// where i is the matrix row index and j is the matrix column index.
// interesting term is n_skipped_idxs, this is the number of elements
// that referencing (i, j) skip when only storing upper triangle. For
// example, consider the following table with matrix to array indexes
//
//  i       j       num_skipped     idx     
//  0       0       0               0
//  0       5       0               5
//  1       0       0               1n - 0
//  2       0       1               2n - 1
//  3       0       3               3n - 3
//  4       0       6               4n - 6
//  
// we see that number skipped is the number of lower triangular elements
// of a matrix constructed from i rows,  (i-1) * i / 2.  Here, we see an
// obvious problem, that when i = 0 we get a negative number, which doesn't
// make sense.  This can be avoided by using the equivalent formulat
//
// n_skipped_idxs = i * (i + 1) / 2 - i
//
// making the equation above read
//  
//  i * (n_samples + 1) - i * (i+1)/2 + j
#define MATRIX_IDX_TO_ARRAY(i, j, n)   ((i) * (n + 1) - (i)*(i+1)/2 + j)


namespace grm {

enum STATUS { 
    SUCCESS, 
    FAILED, 
    UNKNOWN_FAILURE,
    ERROR_IDX_ARR_BOUNDS,
    ERROR_FOPEN,
    ERROR_EOF_NOT_REACHED,
    ERROR_ON_WRITE,
    ERROR_FILE_NOT_OPEN,
};


enum GrmType {
    EHC,        // Expected Haplotype Count
    EAC,        // Expected Alternative Allele Count
    BOTH,       // Both EHC AND EAC
    DS,         // Dosage, i.e. Called Alternative Allele Count
}; 

// @title: Store genomic coordinates used in GRM calculation
// @description: 
struct Coordinates {
    const std::string contig;
    const size_t len;
    std::unique_ptr<size_t> *pos;
};


STATUS write(io::FileIO* fio, Coordinates* coords);
STATUS read(io::FileIO* fio, Coordinates* coords);


// Storage in binary format.  Sample names are comma separated
// [size_t len][names[0],names[1],names[2]...names[len-1]\0]
struct Samples {
    Samples(const size_t len): 
        len(len), 
        names(len == 0 ? nullptr : std::make_unique<std::string>(len)){};
    const size_t len;               //number of samples
    std::unique_ptr<std::string> names;
};

STATUS write(io::FileIO *fio, Samples *samples);
STATUS read(io::FileIO *fio, Samples *samples);


// Header
struct Hdr {
    // const std::string program_version;
    const size_t n_samples;
    const GrmType grm_type;
    const Coordinates *coords;
    const Samples *samples;
};


STATUS write(io::FileIO *fio, const Hdr *hdr);
STATUS read(io::FileIO *fio, Hdr *hdr);


// Grm class manages storage and access of GRM matrix
// 
// The GRM as an n_sample by n_sample symmetric, positive semi-definite
// matrix.  Let Z represent the n_sample by m_marker data genetic data.  
// From these data the GRM is computed as GRM = ZZ^T.
//
// @param n_samples of the GRM.
//
class Grm {
public:
    // 
    Grm(const size_t n_samples);

    Grm(const Grm&)=delete;                          
    Grm(Grm&&)=delete;
    Grm& operator=(const Grm&)=delete;
    Grm& operator=(Grm&&)=delete;
                                            
    // Unchecked indexes when setting and getting of matrix values
    float operator()(const size_t i, const size_t j) const;
    float& operator()(const size_t i, const size_t j);

    // Checked indexes when setting and getting of matrix values
    STATUS set(const size_t i, const size_t j, const float val); 
    STATUS get(const size_t i, const size_t j, float *val) const; 

    size_t size() const;

private:
    const size_t n_samples_;
    std::unique_ptr<float[]> data_;
    size_t midx_to_arr_(const size_t&, const size_t&) const;
};

// @title: Write meta-data and computed grm elements to file
// @description: The binary file written contains a header and payload:
//      * Header
//          - an instance of grm::Header
//      * Payload
//          - grm data in row major order
// @param filename: name of file that the data are written
// @param hdr: an instance of grm::Hdr with important meta data
// @return grm::STATUS: 
//
STATUS write(io.FileIO *fio, const Hdr *hdr, const Grm *grmatrix) const;
STATUS read(io.FileIO *fio, const Hdr *hdr, Grm *grmatrix);


}

#endif
