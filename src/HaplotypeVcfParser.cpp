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
    : fname_(filename),
        file_stream_(fopen(filename, "r")) {

    // check whether file stream open correctly 
    if (!file_stream_)
        throw std::runtime_error("File Access error");
    else if (feof(file_stream_))
        throw std::runtime_error("File is empty");

    // get number of characters in data record for line buffer size
    size_t nchar { get_line_num_char_() };


    // make buffer 10% larger then the number of characters read.
    line_buffer_size_ = static_cast<size_t>(nchar * 1.1);
    line_buffer_ = static_cast<char*>(malloc(line_buffer_size_));
    line_buffer_[0] = '\0';

    set_params_();
};


// HaplotypeVcfParser::HaplotypeVcfParser(std::string filename)
//     : fname_(filename),
//         file_stream_(filename) {
// 
//     if (file_stream_.bad())
//         throw std::runtime_error("File Access error");
//     else if (file_stream_.eof())
//         throw std::runtime_error("File is empty");
// 
//     // get number of characters in data record for line buffer size
//     size_t nchar { get_line_num_char_() };
// 
// 
//     // make buffer 10% larger then the number of characters read.
//     line_buffer_size_ = static_cast<size_t>(nchar * 1.1);
//     line_buffer_ = new char[line_buffer_size_];
//     line_buffer_[0] = '\0';
// 
//     pos_(std::ios_base::beg);
//     set_params_();
// };


HaplotypeVcfParser::~HaplotypeVcfParser() { 
    if(file_stream_)
        fclose(file_stream_);
    free(line_buffer_);
}


size_t HaplotypeVcfParser::get_line_num_char_() {

    bool line_record_found { false };

    size_t char_count { 0 };
    size_t max_char_count { 0 };
    int c { 0 };
    int prev { 'a' };
    while ((c = getc(file_stream_)) != EOF) {

        char_count++;
        
        if (prev == '\n' && c != META_PREFIX)
            line_record_found = true;

        if (line_record_found && c == '\n')
            break;

        if (c == '\n' && char_count > max_char_count)
            max_char_count = char_count;

        prev = c;
    }

    return char_count;
}


size_t HaplotypeVcfParser::n_samples() const { return n_samples_; }


size_t HaplotypeVcfParser::k_founders() const { return k_founders_; }


void HaplotypeVcfParser::pos_(long n) { 
    if (n == 0)
        rewind(file_stream_);

    int c = 0;
    if ((c = fseek(file_stream_, n, SEEK_SET)) != 0)
        throw std::runtime_error("Failed to relocate file stream to position.");
}


void HaplotypeVcfParser::set_params_() {
    pos_(0);

    // skip meta data lines
    ssize_t n { 0 };
    while ((n = getline(&line_buffer_,
                    &line_buffer_size_,
                    file_stream_)) != EOF || n < 0) {

        if ( line_buffer_[0] == META_PREFIX && line_buffer_[1] == META_PREFIX)
            continue;

        break;
    }

    if (ferror(file_stream_) || n < 0 && !feof(file_stream_))
        throw std::runtime_error("File error");

    // if there is no header
    if (line_buffer_[0] != META_PREFIX)
        return;


    if (std::isspace(line_buffer_[0]))
        throw std::runtime_error("First element of VCF line must not be blank.");


    // Get column number and sample number
    StringRecord line_parser_ { SPACE_DELIM, line_buffer_ };
    StringRecord field_parser_ { MEASUREMENT_DELIM };
    StringRecord hap_parser_ { HAP_DELIM };

    n_cols_ = 0;
    n_samples_ = 0;
    for (; line_parser_.next_field(); n_cols_++) {

        std::cout << line_parser_.data() <<std::endl;

        if (n_cols_ < NUM_VCF_FIELDS 
                && std::strcmp(line_parser_.data(), VCF_FIELD_NAMES[n_cols_]) != 0)
            throw std::runtime_error("File doesn't follow vcf header specification");

        if (n_cols_ >= NUM_VCF_FIELDS)
            n_samples_++;

    }

    // get k founders from record
    if ((n = getline(&line_buffer_, 
                        &line_buffer_size_,
                        file_stream_)) == -1 && !feof(file_stream_))
        throw std::runtime_error("File read error");

    std::cout << line_buffer_ << std::endl;

    line_parser_.update_str(line_buffer_);
    size_t hap_idx { 0 };
    for (int i = 0; line_parser_.next_field(); i++) {

        std::cout << line_parser_.data() << std::endl;

        if (i == NUM_VCF_FIELDS-1) {
            field_parser_.update_str(line_parser_.data());

            for (;field_parser_.next_field(); hap_idx++) {
                std::cout << field_parser_.data() << std::endl;
                if (std::strcmp(field_parser_.data(), HAP_CODE))
                    break;
            }

        } else if(i == NUM_VCF_FIELDS) {
            field_parser_.update_str(line_parser_.data()); 

            for(int i = 0; field_parser_.next_field() && i < hap_idx; hap_idx++)
                ;

            hap_parser_.update_str(field_parser_.data());
            for (;hap_parser_.next_field(); k_founders_++)
                ;

            break;
        }

    }

    if (k_founders_ == 0)
        throw std::runtime_error("Parse error");

}


bool HaplotypeVcfParser::load_record(HaplotypeDataRecord& record) {

    ssize_t n { 0 };

    if ((n = getline(&line_buffer_, &line_buffer_size_, file_stream_)) == -1) {
        if (feof(file_stream_))
            return false;
        throw std::runtime_error("File read error");
    }

    record.parse_vcf_line(line_buffer_);

    return true;
}


