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

const static size_t DEFAULT_BUFFER_SIZE { 100000 };




HaplotypeVcfParser::HaplotypeVcfParser(char* filename)
    : fname_(filename),
        file_io_(BufferedRead(filename, DEFAULT_BUFFER_SIZE)) {

    // get number of characters in data record for line buffer size
    size_t nchar { get_line_num_char_() };

    if (nchar == 0)
        throw std::runtime_error("No data to read");

    // make buffer 10% larger then the number of characters read.
    line_buffer_size_ = static_cast<size_t>(nchar * 1.1);
    line_buffer_.reset(line_buffer_size_);

    set_params_();
    pos_(fpos_record_one_);
};


HaplotypeVcfParser::HaplotypeVcfParser(char* filename, size_t buff_size)
    : fname_(filename),
        file_io_(BufferedRead(filename, buff_size)) {

    // get number of characters in data record for line buffer size
    size_t nchar { get_line_num_char_() };

    if (nchar == 0)
        throw std::runtime_error("No data to read");

    std::cout << "Num char: " << nchar << std::endl;
    // make buffer 10% larger then the number of characters read.
    line_buffer_size_ = static_cast<size_t>(nchar * 1.1);
    line_buffer_.reset(line_buffer_size_);

    set_params_();
    pos_(fpos_record_one_);
};


// HaplotypeVcfParser::HaplotypeVcfParser(std::string filename)
//     : fname_(filename),
//         fid_(filename) {
// 
//     if (fid_.bad())
//         throw std::runtime_error("File Access error");
//     else if (fid_.eof())
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


//HaplotypeVcfParser::~HaplotypeVcfParser() { 
//    if(fid_)
//        fclose(fid_);
//    // free(line_buffer_);
//}


size_t HaplotypeVcfParser::get_line_num_char_() {

    size_t char_count { 0 };
    size_t max_char_count { 0 };
    char c;

    while ((c = file_io_.get_char()) != '\0') {
        char_count++;

        if (c == '\n' && char_count > max_char_count) {
            max_char_count = char_count;
            char_count = 0;
        }
    }

    return max_char_count;
}


size_t HaplotypeVcfParser::n_samples() const { return n_samples_; }


size_t HaplotypeVcfParser::k_founders() const { return k_founders_; }


void HaplotypeVcfParser::pos_(size_t n) { 
    file_io_.reset();
    if (n != 0)
        file_io_.seek(n);
}


void HaplotypeVcfParser::set_params_() {
    pos_(0);

    // skip meta data lines

    size_t n { 0 };
    size_t num_bytes { 0 };
    while ((n = file_io_.get_line(line_buffer_)) > 0) {

        // the +1 is because I don't write newline characters to buffer
        num_bytes += line_buffer_.size()+1;

        if (line_buffer_(0) == META_PREFIX && line_buffer_(1) == META_PREFIX)
            continue;

        break;
    }

    // Store file position of first record
    fpos_record_one_ = sizeof(line_buffer_(0)) * num_bytes;

    // if there is no header
    if (line_buffer_(0) != META_PREFIX)
        return;


    if (std::isspace(line_buffer_(0)))
        throw std::runtime_error("First element of VCF line must not be blank.");


    // Get column number and sample number
    StringRecord line_parser_ { SPACE_DELIM, line_buffer_.data() };
    StringRecord field_parser_ { MEASUREMENT_DELIM };
    StringRecord hap_parser_ { HAP_DELIM };

    n_cols_ = 0;
    n_samples_ = 0;
    for (; line_parser_.next_field(); n_cols_++) {

        if (n_cols_ < NUM_VCF_FIELDS 
                && std::strcmp(line_parser_.data(), VCF_FIELD_NAMES[n_cols_]) != 0)
            throw std::runtime_error("File doesn't follow vcf header specification");

        if (n_cols_ >= NUM_VCF_FIELDS)
            n_samples_++;

    }


    // get k founders from record
    if ((n = file_io_.get_line(line_buffer_)) == 0)
        throw std::runtime_error("End of file");

    line_parser_.update_str(line_buffer_.data());
    size_t hap_idx { 0 };
    for (int i = 0; line_parser_.next_field(); i++) {

        if (i == NUM_VCF_FIELDS-1) {
            field_parser_.update_str(line_parser_.data());

            for (;field_parser_.next_field(); hap_idx++)
                if (std::strcmp(field_parser_.data(), HAP_CODE) == 0)
                    break;

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

    size_t n { 0 };

    if ((n = file_io_.get_line(line_buffer_)) == 0)
        return false;

    record.parse_vcf_line(line_buffer_.data());

    return true;
}


