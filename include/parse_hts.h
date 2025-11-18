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
#ifndef HEADER_PARSE_HTS_H
#define HEADER_PARSE_HTS_H

#include <matrix.h>
#include <string>
#include <array>
#include <cstdlib>

namespace htslib {
extern "C" {
#include <htslib/vcf.h>
#include <htslib/hts.h>
}
}



// samples are separated by white space
const char HAP_CODE[] { "HD" };
// const char META_PREFIX { '#' };
// const char MEASUREMENT_DELIM { ':' };
// const char HAP_DELIM { ',' };
// const int NUM_VCF_FIELDS { 9 };
// const char SPACE_DELIM { '\t' };
// 
// 
// // NOTE: in the future it may be best to test for set membership
// static const char* VCF_FIELD_NAMES[NUM_VCF_FIELDS] {
//     "#CHROM",
//     "POS",
//     "ID",
//     "REF",
//     "ALT",
//     "QUAL",
//     "FILTER",
//     "INFO",
//     "FORMAT"
// };
// 
// 
// // Move semantics, I don't want to copy data
// class HaplotypeDataRecord
// {
// public:
// 
//     HaplotypeDataRecord()=delete;
//     HaplotypeDataRecord(size_t, size_t);
//     HaplotypeDataRecord(const HaplotypeDataRecord&)=delete;
//     HaplotypeDataRecord(HaplotypeDataRecord&&)=delete;
// 
// 
//     const std::string& chrom() const;
//     const long pos() const;
//     const std::string& id() const;
//     const char ref() const;
//     const char alt() const;
//     const std::string& qual() const;
//     const std::string& filter() const;
//     const std::string& info() const;
//     const std::string& format() const;
// 
//     void parse_vcf_line(const char*);
//     const double& operator()(size_t, size_t) const;
// 
//     std::array<size_t,2> dims() const;
// 
// 
// private:
//     size_t n_samples_;
//     size_t k_founders_;
// 
//     std::string chrom_ { "" };
//     long pos_ { -1 };
//     std::string id_ { "" };
//     char ref_ { '\0' };
//     char alt_ { '\0' };
//     std::string qual_ { "" };
//     std::string filter_ { "" };
//     std::string info_ { "" };
//     std::string format_ { "" };
// 
//     std::unique_ptr<Matrix> samples_ { nullptr };
// 
//     StringRecord line_parse_ { SPACE_DELIM };
//     StringRecord field_parse_ { MEASUREMENT_DELIM };
//     StringRecord hap_parse_ { HAP_DELIM };
// };

// Interface with htslib bcf tools
class ParseHtsVariantFile
{
public:
    // HaplotypeVcfParser(const char* variant_fname);
    ParseHtsVariantFile(const char *variant_fname, const char *sample_fname);
    // HaplotypeVcfParser(const std::string& variant_fname);
    // HaplotypeVcfParser(const std::string& variant_fname,
    //        const std::string& sample_fname);

    ParseHtsVariantFile()=delete; 
    ParseHtsVariantFile(const ParseHtsVariantFile&)=delete;
    ParseHtsVariantFile(const ParseHtsVariantFile&&)=delete;
    // HaplotypeVcfParser& operator=(const HaplotypeVcfParser&)=delete;

    ~ParseHtsVariantFile();

    // size_t n_samples() const;
    // size_t k_founders() const;

    // bool load_record(HaplotypeDataRecord&);

private:
    const std::string fname_;
    htslib::htsFile *fid_;
    htslib::bcf_hdr_t *hdr_;

    // size_t n_cols_ { 0 };
    // size_t n_samples_ { 0 };
    // size_t k_founders_ { 0 };
    // size_t fpos_record_one_ { 0 };


    // void pos_(size_t);
    // size_t get_line_num_char_();
    // void set_params_();
};

#endif
