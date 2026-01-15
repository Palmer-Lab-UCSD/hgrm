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
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <memory>
#include <utility>
#include <string>


namespace grm {

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

}


enum class STATUS { 
    SUCCESS, 
    FAILED, 
    UNKNOWN_FAILURE,
    ERROR_IDX_ARR_BOUNDS,
    ERROR_FOPEN,
    ERROR_EOF_NOT_REACHED,
    ERROR_ON_WRITE,
};


struct Dims {
    Dims(size_t nrow_in, size_t mcol_in): 
        nrow(nrow_in), mcol(mcol_in) {};

    const size_t nrow;
    const size_t mcol;
};


// @title: Store genomic coordinates and mange binary I/O
// @description: 
struct Coordinates {

    // @param pos_filename: The name, and path if necessary, of the text
    //      file specifying variant positions on the specified contig 
    //      to be included for the grm.
    // @param contig_name: The name of the contig, e.g. chrm1
    //
    Coordinates(const char *contig_name, const size_t len);
    Coordinates(const std::string& contig_name, const size_t len);

    size_t operator[](size_t idx) const;
    size_t& operator[](size_t idx);

    STATUS write(FILE *fid);
    static STATUS read(FILE *fid, Coordinates *coords);

    const size_t len;
    const std::string contig;
    std::unique_ptr<size_t> *pos;
};


struct Samples {
    const size_t len;
    std::unique_ptr<char*> names;

    STATUS write(FILE *fid);
    static STATUS read(FILE *fid, Samples *samples);
};


STATUS load_samples(const char *filename, Samples *samples);


struct Hdr {
    const std::string program_version;
    const std::string data_type;
    const Coordinates *coords;
    const Samples *samples;

    STATUS bin_write(FILE *fid);
    static STATUS bin_read(FILE *fid, Hdr *hdr);
};


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
    STATUS write(const char *filename, const Hdr *hdr) const;

    static STATUS read(const char *filename, Grm *grm);

private:
    const Dims dims_;
    std::unique_ptr<float[]> data_;
    size_t midx_to_arr_(const size_t&, const size_t&) const;
};

}

#endif
