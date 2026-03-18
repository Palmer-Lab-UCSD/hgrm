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
// reviewed by Claude Opus 4.6, the AI assistant from Anthropic.
// Some recommendations have been incorporated.
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

#include "constants.h"


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
    ERROR_ON_READ,
    ERROR_FILE_NOT_OPEN,
    ERROR_NULLPTR_ARG,
    ERROR_INVALID_ARG,
};


enum GrmType {
    EHC,        // Expected Haplotype Count
    EAC,        // Expected Alternative Allele Count
    BOTH,       // Both EHC AND EAC
    DS,         // Dosage, i.e. Called Alternative Allele Count
}; 

// @title: Store genomic coordinates used in GRM calculation
struct Coordinates {
    Coordinates(): contig(""), len(0), pos(nullptr) {};
    Coordinates(const char* contig, const size_t len)
        : contig(contig),
        len(len), 
        pos(std::make_unique<size_t>(len)) {};

    Coordinates(Coordinates&) = delete;
    Coordinates& operator=(Coordinates&) = delete;

    Coordinates(Coordinates&& other);
    Coordinates& operator=(Coordinates&& other);

    // Data Fields
    std::string contig;
    size_t len;
    std::unique_ptr<size_t[]> pos;
};


// Coordinates Storage Layout
//
//  type    number  description
//  --------------------------------------------------------------------
//  size_t  1       number of characters (n) in contig name
//  char    n       characters for contig name without null character
//  size_t  1       number of genomic positions (npos) 
//  size_t  npos    the positions used for computation of the grm
//
STATUS write(io::FileIO* fio, const Coordinates* coords);
STATUS read(io::FileIO* fio, Coordinates* coords);


// Samples stores sample id strings and the number of samples
//
struct Samples {
    Samples(): len(0), names(nullptr);
    Samples(size_t n_samples): 
        len(n_samples), 
        names(len == 0 ? nullptr : std::make_unique<std::string[]>(len)) {};

    Samples(const Samples&) = delete;
    Samples& operator=(const Samples&) = delete;

    Samples(Samples&& other);
    Samples& operator=(Samples&& other);

    // Data Fields
    size_t len;               //number of samples
    std::unique_ptr<std::string[]> names;

};


// Sample Storage Layout
//
//  type    number  description
//  --------------------------------------------------------------------
//  size_t  1       represents number of samples
//  size_t  1       the number of characters of longest sample id
//  size_t  1       number of characters (n_1) in first sample id
//  char    n_1     characters of sample id 1 without terminal null '\0'
//  size_t  1       number of characters (n_2) in second sample id
//  char    n_2     characters of sample id 2 without terminal null '\0'
//  ...
//  size_t  1       number of characters (n_N) in N^{th} sample id
//  char    n_N     characters of sample id N without terminal null '\0'

STATUS write(io::FileIO* fio, const Samples* samples);
STATUS read(io::FileIO* fio, Samples* samples);


// Header
struct Hdr {

    // Data Fields
    GrmType grm_type;
    Coordinates* coords;
    Samples* samples;

    const std::string version = constants::version;
};

// Header Storage Layout
//
//  type    number  description
//  --------------------------------------------------------------------
//  size_t  1       number of characters (n) in version string
//  char    n       version string without null terminator
//  GrmType 1       type of grm 
//
//  call write for coordinates
//
//  call write for samples
//
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
