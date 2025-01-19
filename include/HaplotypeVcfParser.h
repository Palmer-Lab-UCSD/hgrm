// Parse STITCH vcf file
//
//
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
#ifndef HEADER_HAPLOTYPEVCFPARSER_H
#define HEADER_HAPLOTYPEVCFPARSER_H

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "Matrix.h"
#include "utils.h"


// samples are separated by white space
const char HAP_CODE[] { "HD" };
const char META_PREFIX { '#' };
const std::string MEASUREMENT_DELIM { ":" };
const std::string HAP_DELIM { "," };
const int NUM_VCF_FIELDS { 9 };
const std::string SPACE_DELIM { " \t\n\v\f\r" };
const size_t BUFFER_SIZE { 1000 };

// NOTE: in the future it may be best to test for set membership
static const char* VCF_FIELD_NAMES[NUM_VCF_FIELDS] {
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT"
};


// Move semantics, I don't want to copy data
class HaplotypeDataRecord
{
public:

    HaplotypeDataRecord()=delete;
    HaplotypeDataRecord(size_t n_samples, size_t k_founders);
    HaplotypeDataRecord(const HaplotypeDataRecord&)=delete;
    HaplotypeDataRecord(HaplotypeDataRecord&&)=delete;


    const std::string& chrom() const;
    const long pos() const;
    const std::string& id() const;
    const char ref() const;
    const char alt() const;
    const std::string& qual() const;
    const std::string& filter() const;
    const std::string& info() const;
    const std::string& format() const;

    void parse_vcf_line(const char*);
    double operator()(size_t, size_t) const;

    std::array<size_t,2> dims() const;


private:
    size_t n_samples_;
    size_t k_founders_;

    std::string chrom_ { "" };
    long pos_ { -1 };
    std::string id_ { "" };
    char ref_ { '\0' };
    char alt_ { '\0' };
    std::string qual_ { "" };
    std::string filter_ { "" };
    std::string info_ { "" };
    std::string format_ { "" };

    std::unique_ptr<Matrix> samples_ { nullptr };

    StringRecord line_parse_ { SPACE_DELIM };
    StringRecord field_parse_ { MEASUREMENT_DELIM };
    StringRecord hap_parse_ { HAP_DELIM };
};


class HaplotypeVcfParser
{
public:

    HaplotypeVcfParser()=delete;                                // default constructor
    HaplotypeVcfParser(char* filename);                         // constructor
    //HaplotypeVcfParser(std::string filename);                   // constructor
    HaplotypeVcfParser(const HaplotypeVcfParser&)=delete;       // copy constructor
    HaplotypeVcfParser(const HaplotypeVcfParser&&)=delete;       // move constructor
    HaplotypeVcfParser& operator=(const HaplotypeVcfParser&)=delete;    // copy assignment
    ~HaplotypeVcfParser();                                      // descructor

    size_t n_samples() const;
    size_t k_founders() const;

    bool load_record(HaplotypeDataRecord&);

private:
    const std::string fname_;
    FILE* file_stream_;

    size_t n_cols_ { 0 };
    size_t n_samples_ { 0 };
    size_t k_founders_ { 0 };

    char* line_buffer_ { nullptr };
    size_t line_buffer_size_ { 0 };

    void pos_(long);
    size_t get_line_num_char_();
    void set_params_();
};

#endif
