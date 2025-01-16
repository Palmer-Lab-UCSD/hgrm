//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
// Input argument
//    filename: vcf with haplotpye
//
//
// Acknowledgment
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Sonnnet, the AI assistant from Anthropic
// (Jan 2025), with minor recommendations incorporated.

#include "HaplotypeVcfParser.h"

HaplotypeVcfParser::HaplotypeVcfParser(char* filename)
    : fname_(filename), file_stream_(filename), 
        n_cols_(0), n_samples_(0), n_meta_lines_(0), 
        meta_(nullptr), vcf_standard_colnames_(nullptr), sample_names_(nullptr) {

    if (file_stream_.bad())
        throw std::runtime_error("File Access error");
    else if (file_stream_.eof())
        throw std::runtime_error("File is empty");

    set_vcf_params_();
};


HaplotypeVcfParser::HaplotypeVcfParser(std::string filename)
    : fname_(filename), file_stream_(filename), 
        n_cols_(0), n_samples_(0), n_meta_lines_(0), 
        meta_(nullptr), vcf_standard_colnames_(nullptr), sample_names_(nullptr) {

    if (file_stream_.bad())
        throw std::runtime_error("File Access error");
    else if (file_stream_.eof())
        throw std::runtime_error("File is empty");

    set_vcf_params_();
};


HaplotypeVcfParser::~HaplotypeVcfParser() { file_stream_.close(); }


void HaplotypeVcfParser::pos_(int n) { 
    file_stream_.clear();
    file_stream_.seekg(n);
}


void HaplotypeVcfParser::load_meta_() {
    // ensure beginning of file
    pos_(std::ios_base::beg);

    std::string tmp {};

    // count meta lines
    n_meta_lines_ = 0;
    while ( std::getline(file_stream_, tmp) ) {

        // break at VCF header
        if ( tmp[0] == META_PREFIX && tmp[1] != META_PREFIX)
            break;

        n_meta_lines_++;
    }

    if (file_stream_.bad())
        throw std::runtime_error("File read error");

    // instantiate meta array
    meta_ = std::make_unique<std::string[]>(n_meta_lines_);


    // load meta lines
    pos_(std::ios_base::beg);
    for (int i = 0; i < n_meta_lines_; i++) {

        if (!std::getline(file_stream_, tmp) && !file_stream_.eof())
            throw std::runtime_error("File read error");

        meta_[i]=tmp;
    }
}


void HaplotypeVcfParser::load_header_() {
    pos_(std::ios_base::beg);

    // skip meta data lines
    std::string s;
    while ( std::getline(file_stream_, s) ) {

        if ( s[0] == META_PREFIX && s[1] == META_PREFIX)
            continue;

        break;
    }

    if (file_stream_.bad())
        throw("File read error");


    // if there is no header
    if (s[0] != META_PREFIX) {
        n_cols_ = 0;
        n_samples_ = 0;
        vcf_standard_colnames_ = nullptr;
        sample_names_ = nullptr;

        return;
    }

    if (std::isblank(s[0]))
        throw std::runtime_error("First element of VCF line must not be blank.");

    // analyze the mandatory vcf columns, not I require format to be
    // column number 9, then assume that column 10 on are samples
    n_cols_ = 0;
    vcf_standard_colnames_ = std::make_unique<std::string[]>(NUM_VCF_FIELDS);
    size_t i { 0 };
    for (;i < s.size() && n_cols_ < NUM_VCF_FIELDS
            && load_next_field_to_buffer_(s, i); i++) {

        vcf_standard_colnames_[n_cols_] = buffer_.data();

        if (vcf_standard_colnames_[n_cols_] != VCF_FIELD_NAMES[n_cols_])
            throw std::runtime_error("File does not adhere to VCF standard.");

        n_cols_++;
    }

    if (n_cols_ != NUM_VCF_FIELDS)
        throw std::runtime_error("File does not adhere to VCF standard.");
    
    size_t curr_file_position_ { static_cast<size_t>(file_stream_.tellg()) };
    size_t j { i };
    
    for (;j < s.size(); j++){
        load_next_field_to_buffer_(s, j);
        n_cols_++;
    }

    n_samples_ = n_cols_ - NUM_VCF_FIELDS;

    sample_names_ = std::make_unique<std::string[]>(n_samples_);
    size_t samp_count { 0 };
    pos_(curr_file_position_);

    for (; i < s.size() && load_next_field_to_buffer_(s, i); i++) {
        sample_names_[samp_count] = buffer_.data();
        samp_count++;
    }

    if (samp_count != n_samples_)
        throw std::runtime_error("Number of samples read do not agree");

}


bool HaplotypeVcfParser::load_next_field_to_buffer_(const std::string& s, size_t& i) {

    buffer_idx_ = 0;

    for (;i < s.size();i++) {

        if (i > 0 && std::isspace(s[i]) && !std::isspace(s[i-1])) {
            buffer_[buffer_idx_] = '\0';
            break;
        }

        // more than one space in a row
        if (std::isspace(s[i]))
            continue;

        //if (!std::isspace(s[i]))
        //    break;

        buffer_[buffer_idx_++] = s[i];
    }
     
    if (buffer_idx_ == 0)
        return false;

    return true;
}


void HaplotypeVcfParser::set_vcf_params_()
{
    load_meta_();
    load_header_();
}


bool HaplotypeVcfParser::load_record(HaplotypeDataRecord& record) {
    if (file_stream_.eof())
        return false;

    std::string s {};
    if (!std::getline(file_stream_, s) && !file_stream_.eof())
        throw("File read error");


    record.parse_vcf_line(s, n_cols_);

    return true;
}
