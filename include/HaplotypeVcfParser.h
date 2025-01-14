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
#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <cstdlib>
#include "Matrix.h"


// samples are separated by white space
const char HAP_CODE[] { "HD" };
const char META_PREFIX { '#' };
const char MEASUREMENT_DELIM { ':' };
const char HAP_DELIM { ',' };
const int NUM_VCF_FIELDS { 9 };
const std::array<std::string, NUM_VCF_FIELDS> VCF_FIELD_NAMES= {
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

    HaplotypeDataRecord();
    HaplotypeDataRecord(const std::string& vcf_line,
                        size_t n_cols);                 // constructor
    HaplotypeDataRecord(const HaplotypeDataRecord&);    // copy constructor
    //HaplotypeDataRecord(HaplotypeDataRecord&&);         // move constructor
    HaplotypeDataRecord& operator=(const HaplotypeDataRecord&);  //copy assignment
    //HaplotypeDataRecord& operator=(HaplotypeDataRecord&&)       //move assignment

    const std::string& chrom() const;
    const long pos() const;
    const std::string& id() const;
    const char ref() const;
    const char alt() const;
    const std::string& qual() const;
    const std::string& filter() const;
    const std::string& info() const;
    const std::string& format() const;

    void parse_vcf_line(const std::string&, size_t n_cols);
    double operator()(size_t, size_t) const;

    std::array<size_t,2> dims() const;


private:
    static constexpr size_t buffer_size_ { 1000 };
    size_t buffer_idx_ { 0 };
    std::array<char,buffer_size_> buffer_;

    static constexpr size_t hap_buffer_size_ { 100 };
    size_t hap_buffer_idx_ { 0 };
    std::array<char,hap_buffer_size_> hap_buffer_;


    std::unique_ptr<Matrix> samples_;

    std::string chrom_;
    long pos_;
    std::string id_;
    char ref_;
    char alt_;
    std::string qual_;
    std::string filter_;
    std::string info_;
    std::string format_;

};


class HaplotypeVcfParser
{
public:

    HaplotypeVcfParser(char* filename);
    HaplotypeVcfParser(std::string filename);
    ~HaplotypeVcfParser(); 

    bool get_record(HaplotypeDataRecord&);
    size_t n_samples() { return n_samples_; };

private:
    const std::string fname_;
    const size_t n_cols_;
    const size_t n_samples_;
    const size_t n_meta_lines_;
    const std::unique_ptr<std::string[]> meta_;
    const std::unique_ptr<std::string[]> vcf_standard_colnames_;
    const std::unique_ptr<std::string[]> sample_names_;

    std::ifstream file_stream_;

    static constexpr size_t buffer_size_ { 1000 };
    size_t buffer_idx_ { 0 };
    std::array<char,buffer_size_> buffer_;

    void pos_(int);
    void set_vcf_params_();
    void load_meta_();
    void load_header_();
    bool load_next_field_to_buffer_(const std::string&, size_t&);
};


