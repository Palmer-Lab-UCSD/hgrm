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

// *****************************************************************************
// class BcfHeader
// *****************************************************************************
//
int bcfio::BcfHeader::decode_hts_idinfo_(const char *name, 
        const int bcf_dt_type, 
        bcfio::BcfHdrAttr *ptr) const {

    // BCF_DT_ID is the C macro for the ID dictionary index defined by htslib
    // see htslib/vcf.h line 86
    int idx = htslib::bcf_hdr_id2int(hdr_, BCF_DT_ID, name);

    if (idx < 0)
        return idx;

    uint64_t val = hdr_->id[BCF_DT_ID][idx].val->info[bcf_dt_type];

    ptr->number = val >> 12 & 0xfffff;
    ptr->vl_type = val >> 8 & 0xf;
    ptr->type = val >> 4 & 0xf;
    ptr->coltype = val & 0xf;

    return 0;
}

int bcfio::BcfHeader::get_format_attr(const char *id, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(id, BCF_HL_FMT, ptr);
}

int bcfio::BcfHeader::get_info_attr(const char *id, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(id, BCF_HL_INFO, ptr);
}

int bcfio::BcfHeader::get_filter_attr(const char *id, BcfHdrAttr *ptr) const {
    return decode_hts_idinfo_(id, BCF_HL_FLT, ptr);
}

int32_t bcfio::BcfHeader::k_fmt(const char *id) const {
    if (!id)
        return -1;

    BcfHdrAttr fmt {};

    int32_t status { 0 };

    if ((status = get_format_attr(id, &fmt)) < 0)
        return status;

    return static_cast<int32_t>(fmt.number);
}

int bcfio::BcfHeader::subset_samples(const char *filename) {
    return htslib::bcf_hdr_set_samples(hdr_, filename, 1);
}

const std::unique_ptr<std::string[]> bcfio::BcfHeader::sample_names() const {

    std::unique_ptr<std::string[]> samp_names = 
        std::make_unique<std::string[]>(n_samples()); 

    for (int i = 0; i < n_samples(); i++)
        samp_names[i] = std::string(*(hdr_->samples + i));

    return samp_names;
}

// *****************************************************************************
// class BcfFloatRecord
// *****************************************************************************

bcfio::BcfFloatRecord::~BcfFloatRecord() {
    if (rec_) htslib::bcf_destroy(rec_);
    if (dst_) free(dst_);
    rec_ = nullptr;
    dst_ = nullptr;
}

std::optional<float> bcfio::BcfFloatRecord::get(const size_t row_idx,
        const size_t col_idx) const {
    if ((row_idx * col_idx + col_idx) >= ndst_) return std::nullopt;

    return *(dst_ + row_idx * col_idx + col_idx);
}

int bcfio::BcfFloatRecord::load_data(bcfio::BcfHeader *hdr, const char *id) {
    int status { 0 };
    col_num_ = row_num_ = 0;

    status = htslib::bcf_get_format_values(hdr->hts_hdr(), 
            rec_, 
            id, 
            (void**)(&dst_),
            &ndst_, 
            BCF_HT_REAL);

    if (status < 0)
        return status;

    int32_t k { 0 };
    if ((k = hdr->k_fmt(id)) < 0) {
        col_num_ = row_num_ = 0;
        return k;
    }

    col_num_ = static_cast<size_t>(k);
    row_num_ = hdr->n_samples();

    return 0;
}


// *****************************************************************************
// class BcfRead
// *****************************************************************************
bcfio::ReadBcf::ReadBcf(const char *bcfname)
    : fname_(bcfname),
    fid_(htslib::hts_open(bcfname, "r")),
    hdr_(fid_) {};


// TODO: subset samples by those in sample_fname file
bcfio::ReadBcf::ReadBcf(const char *bcfname, const char *sample_fname)
    : fname_(bcfname),
    fid_(htslib::hts_open(bcfname, "r")),
    hdr_(fid_) {

    int status { 0 };
    // Subset samples with those found in the file sample_fname 
    if (!sample_fname || *sample_fname == '\0')
        fprintf(stdout, "No file with sample names detected, retreiving"
                " records for all samples.\n");
    else
        status = hdr_.subset_samples(sample_fname);

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
    if (fid_) htslib::hts_close(fid_);
}


// title: load next record
int bcfio::ReadBcf::next_record(bcfio::BcfFloatRecord *ptr, const char *id) {
    int status = htslib::bcf_read(fid_, hdr_.hts_hdr(), ptr->cur_rec());
    if (status != 0)
        return status;

    // Unpacking options defined in htslib/vcf.h line 419
    if (htslib::bcf_unpack(ptr->cur_rec(), BCF_UN_ALL) < 0)
        return -1;

    return ptr->load_data(&hdr_, id);
}
