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

namespace bcfio {
// @title The meta data on a BCF attribute
// @description BCF, VCF, and VCF.GZ files hold metadata in the header that
//     specify the type and format of data in records.  I call each unique
//     piece of data in a record a record attribute, e.g. an INFO column or 
//     FORMAT column of a record are attributes of that record.  HTSLIB encodes
//     attribute information an unsigned 64 bit integer, and to access any value
//     one needs to correctly implement bit shifting and masking.  This struct 
//     contains bit-fields representing each value stored in the uint64_t.
//
struct BcfHdrAttr {
    uint64_t number : 20;
    uint64_t vl_type : 4;
    uint64_t type : 4;
    uint64_t coltype : 4;
};


class BcfHeader {
public:
    htslib::bcf_hdr_t *hdr;
    
    BcfHeader(htslib::htsFile *fid): 
        hdr(fid ? htslib::bcf_hdr_read(fid) : nullptr) {};
    ~BcfHeader() { if (hdr) htslib::bcf_hdr_destroy(hdr); };

    const bool isnull() const { return hdr == nullptr; };

    // sample_names()
    const int get_format(const char *name, BcfHdrAttr *ptr) const;
    const int get_info(const char *name, BcfHdrAttr *ptr) const;
    const int get_filter(const char *name, BcfHdrAttr *ptr) const;

private:
    const int decode_hts_idinfo_(const char *name, 
            const int bcf_dt_type, 
            BcfHdrAttr *ptr) const;
};



// Interface with htslib bcf tools
class ReadBcf
{
public:
    // HaplotypeVcfParser(const char* variant_fname);
    ReadBcf(const char *variant_fname);
    ReadBcf(const char *variant_fname, const char *sample_fname);
    // HaplotypeVcfParser(const std::string& variant_fname);
    // HaplotypeVcfParser(const std::string& variant_fname,
    //        const std::string& sample_fname);

    ReadBcf()=delete; 
    ReadBcf(const ReadBcf&)=delete;
    ReadBcf(const ReadBcf&&)=delete;
    // HaplotypeVcfParser& operator=(const HaplotypeVcfParser&)=delete;

    ~ReadBcf();

    const size_t n_samples() const;
    const size_t k_founders() const;
    std::unique_ptr<std::string[]> sample_names() const;

    // bool load_record(HaplotypeDataRecord&);

private:
    const std::string fname_;
    htslib::htsFile *fid_;
    BcfHeader hdr_;

    // size_t n_cols_ { 0 };
    // size_t n_samples_ { 0 };
    // size_t k_founders_ { 0 };
    // size_t fpos_record_one_ { 0 };


    // void pos_(size_t);
    // size_t get_line_num_char_();
    // void set_params_();
};
}

#endif
