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


namespace grm {

enum STATUS { 
    SUCCESS, 
    FAILED, 
    UNKNOWN_FAILURE,
    ERROR_IDX_ARR_BOUNDS,
    ERROR_FOPEN,
    ERROR_EOF_NOT_REACHED,
    ERROR_ON_WRITE,
};


namespace details {

    // @title: Count the number of non empty lines in text file
    //
    // @param fid: pointer to C file stream, i.e. that returned by fopen
    // @param num_lines: the number of lines written at this address
    // @return  -1: file I/O error as determined by ferror(fid), or
    //          -2: end of file not reached, reason undetermined, or
    //          -3: error in returning file handle to beginning of file
    //           0: success
    int num_lines_in_file(FILE *fid, size_t *num_lines);


    // int chars_to_size_t(FILE *fid, size_t *val);
    //
    
}


struct Dims {
    Dims(size_t nrow_in, size_t mcol_in): 
        nrow(nrow_in), mcol(mcol_in) {};

    const size_t nrow;
    const size_t mcol;
};


// @title: Store genomic coordinates used in GRM calculation
// @description: 
struct Coordinates {
    const std::string contig;
    const size_t len;
    std::unique_ptr<size_t> *pos;
};


STATUS write(io::FileIO *fio, Coordinates *coords);
STATUS read(io::FileIO *fio, Coordinates *coords);


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


struct Hdr {
    const std::string program_version;
    const std::string data_type;
    const Coordinates *coords;
    const Samples *samples;
};


STATUS write(FILE *fid, Hdr *hdr);
STATUS read(FILE *fid, Hdr *hdr);

class Grm {
public:
    Grm(const size_t, const size_t);
    Grm(const Grm&);                          // copy constructor
    Grm(Grm&&);                               // move constructor
    Grm& operator=(const Grm&)=delete;        // copy assignment
    Grm& operator=(Grm&&)=delete;             // move assignment
                                            
    // Unchecked indexes when setting and getting of matrix values
    float operator()(const size_t i, const size_t j) const;
    float& operator()(const size_t i, const size_t j);

    // Checked indexes when setting and getting of matrix values
    STATUS set(const size_t i, const size_t j, const float val); 
    STATUS get(const size_t i, const size_t j, float *val) const; 

    size_t size() const;
    const Dims& dims() const;

private:
    const Dims dims_;
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
STATUS write(io.FileIO *fio, Grm *grm, Hdr *hdr) const;
STATUS read(io.FileIO *fio, Grm *grm, Hdr *hdr);

}

#endif
