// Parse STITCH vcf file
//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
//
// Acknowledgment
//
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
// const char HAP_CODE[] { "HD" };

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
struct BcfHdrAttr { uint64_t number : 20, vl_type : 4, type : 4, coltype : 4; };


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


// @title Manage bcf record
// @description Manage the lifetime of a htslib::bcf1_t type record using 
//      htslib functions with RAII.  Provide some a simpler interface to
//      quantities of interest
struct BcfRecord {
    BcfRecord(): rec(htslib::bcf_init()) {};
    ~BcfRecord() { if (rec) htslib::bcf_destroy(rec); };

    bool is_snp() const { return htslib::bcf_is_snp(rec); }
    htslib::bcf1_t *rec;
};


// @title Interface with htslib bcf tools
// @description ReadBCF manages the lifetime of an open htslib file
//      and organizes the bcf file header and any one record for easy
//      and memory safe parsing.
// @param bcfname: the path and filename to the bcf file to be read.
// @param sample_fname: the path and filename of the text file listing the
//      samples id's of records to be retreived.  If this is not included
//      all sample records are retrieved.
class ReadBcf
{
public:
    ReadBcf(const char *bcfname);
    ReadBcf(const char *bcfname, const char *sample_fname);
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

    int next_record(BcfRecord *rec);
    // size_t fpos_record_one_ { 0 };

    // void pos_(size_t);
    // size_t get_line_num_char_();
    // void set_params_();
};
}

#endif
