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

// decoder based upon htslib/vcf.h line 100 in the typedef struct bcf_idinfo_t. 
const int bcfio::BcfHeader::decode_hts_idinfo_(const char *name, 
        const int bcf_dt_type, 
        bcfio::BcfHdrAttr *ptr) const {

    // BCF_DT_ID is the ID dictionary index defined by htslib
    int idx = htslib::bcf_hdr_id2int(hdr, BCF_DT_ID, name);

    if (idx < 0)
        return -1;

    uint64_t val = hdr->id[BCF_DT_ID][idx].val->info[bcf_dt_type];

    ptr->number = val >> 12 & 0xfffff;
    ptr->vl_type = val >> 8 & 0xf;
    ptr->type = val >> 4 & 0xf;
    ptr->coltype = val & 0xf;

    return 0;
}

const int bcfio::BcfHeader::get_format(const char *name, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(name, BCF_HL_FMT, ptr);
}

const int bcfio::BcfHeader::get_info(const char *name, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(name, BCF_HL_INFO, ptr);
}

const int bcfio::BcfHeader::get_filter(const char *name, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(name, BCF_HL_FLT, ptr);
}


bcfio::ReadBcf::ReadBcf(const char *variant_fname)
    : fname_(variant_fname),
    fid_(htslib::hts_open(variant_fname, "r")),
    hdr_(fid_) {};


// TODO: subset samples by those in sample_fname file
bcfio::ReadBcf::ReadBcf(const char *variant_fname, const char *sample_fname)
    : fname_(variant_fname),
    fid_(htslib::hts_open(variant_fname, "r")),
    hdr_(fid_) {

    int status { 0 };
    // Subset samples with those found in the file sample_fname 
    if (!sample_fname || *sample_fname == '\0')
        fprintf(stdout, "No file with sample names detected, computing"
                "hGRM over all samples.\n");
    else
        status = htslib::bcf_hdr_set_samples(hdr_.hdr, sample_fname, 1);

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

const size_t bcfio::ReadBcf::n_samples() const {
    // See htslib/vcf.h line 649
    // Remember that n is the number of entries in the triplet of 
    // dictionaries in the VCF.  BCF_DT_SAMPLE, provides the index of n
    // that correspondes to the number of samples.
    return hdr_.hdr->n[BCF_DT_SAMPLE];
};

const size_t bcfio::ReadBcf::k_founders() const {
    BcfHdrAttr fmt {};

    if (hdr_.get_format("HD", &fmt) < 0)
        printf("errror\n");

    return static_cast<size_t>(fmt.number);
}

std::unique_ptr<std::string[]> bcfio::ReadBcf::sample_names() const {

    std::unique_ptr<std::string[]> samp_names = 
        std::make_unique<std::string[]>(n_samples()); 

    for (int i = 0; i < n_samples(); i++)
        samp_names[i] = std::string(*(hdr_.hdr->samples + i));

    return samp_names;
}
