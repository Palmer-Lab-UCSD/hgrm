//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
// Input argument
//    filename: vcf with haplotpye
//
//
//

#include <bcfio.h>
#include <cstdlib>

// *****************************************************************************
// class BcfHeader
// *****************************************************************************
//
const int bcfio::BcfHeader::decode_hts_idinfo_(const char *name, 
        const int bcf_dt_type, 
        bcfio::BcfHdrAttr *ptr) const {

    // BCF_DT_ID is the C macro for the ID dictionary index defined by htslib
    // see htslib/vcf.h line 86
    int idx = htslib::bcf_hdr_id2int(hdr_, BCF_DT_ID, name);

    if (idx < 0)
        return -1;

    uint64_t val = hdr_->id[BCF_DT_ID][idx].val->info[bcf_dt_type];

    ptr->number = val >> 12 & 0xfffff;
    ptr->vl_type = val >> 8 & 0xf;
    ptr->type = val >> 4 & 0xf;
    ptr->coltype = val & 0xf;

    return 0;
}

int bcfio::BcfHeader::get_format(const char *id, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(id, BCF_HL_FMT, ptr);
}

int bcfio::BcfHeader::get_info(const char *id, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(id, BCF_HL_INFO, ptr);
}

int bcfio::BcfHeader::get_filter(const char *id, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(id, BCF_HL_FLT, ptr);
}

// *****************************************************************************
// class BcfRecord
// *****************************************************************************

bcfio::BcfRecord::~BcfRecord() {
    if (rec) htslib::bcf_destroy(rec);
    if (dst_) free(dst_);
}

int bcfio::BcfRecord::load_data(bcfio::BcfHeader *hdr, const char *id) {
    int status { 0 };

    if ((status = hdr->get_format(id, &attr_)) != 0) 
        return status;

    if (attr_->type == BCF_HT_REAL)
        return htslib::bcf_get_format_values(hdr->hts_hdr(), rec, id, 
                (void**)(&fdst_), &ndst_, BCF_HT_REAL);
    else if (attr_->type == BCF_HT_STR)
        return htslib::bcf_get_format_values(hdr->hts_hdr(), rec, id, 
                (void**)(&cdst_), &ndst_, BCF_HT_STR);

    printf("ERROR: Only types string and float are currently supported."
            " contact the project maintainer if your type is not yet"
            " supported.\n");
    return -1;
}


// *****************************************************************************
// class BcfRead
// *****************************************************************************
bcfio::ReadBcf::ReadBcf(const char *bcfname)
    : fname_(bcfname),
    fid_(htslib::hts_open(bcfname, "r")),
    hdr(fid_) {};


// TODO: subset samples by those in sample_fname file
bcfio::ReadBcf::ReadBcf(const char *bcfname, const char *sample_fname)
    : fname_(bcfname),
    fid_(htslib::hts_open(bcfname, "r")),
    hdr(fid_) {

    int status { 0 };
    // Subset samples with those found in the file sample_fname 
    if (!sample_fname || *sample_fname == '\0')
        fprintf(stdout, "No file with sample names detected, retreiving"
                " records for all samples.\n");
    else
        status = htslib::bcf_hdr_set_samples(hdr.hdr, sample_fname, 1);

    if (status < 0) {
        fprintf(stderr, "Error: Couldn't read sample file\n");
        exit(EXIT_FAILURE);
    } else if (status > 0) {
        fprintf(stderr, "Error: A subset of samples in sample file are not"
                " found in the VCF,BCF, or VCF.GZ file.\n");
        exit(EXIT_FAILURE);
    }


    // get number of characters in data record for line buffer size
};

bcfio::ReadBcf::~ReadBcf() {
    if (fid_)
        htslib::hts_close(fid_);
}


const size_t bcfio::ReadBcf::k_haps() const {
    BcfHdrAttr fmt {};

    if (hdr.get_format("HD", &fmt) < 0)
        printf("error\n");

    return static_cast<const size_t>(fmt.number);
}

// Note: May be better to just return a reference?
std::unique_ptr<std::string[]> bcfio::ReadBcf::sample_names() const {

    std::unique_ptr<std::string[]> samp_names = 
        std::make_unique<std::string[]>(n_samples()); 

    for (int i = 0; i < n_samples(); i++)
        samp_names[i] = std::string(*(hdr.hdr->samples + i));

    return samp_names;
}

// title: load next record
int bcfio::ReadBcf::next_record(bcfio::BcfRecord *ptr) {
    int status = htslib::bcf_read(fid_, hdr.hdr, ptr->rec);
    if (status != 0)
        return status;

    // Unpacking options defined in htslib/vcf.h line 419
    return htslib::bcf_unpack(ptr->rec, BCF_UN_ALL);
}
