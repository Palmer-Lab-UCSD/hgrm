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
#include "HaplotypeVcfParser.h"


HaplotypeDataRecord::HaplotypeDataRecord() 
    : chrom(""), pos(-1), id(""),
    ref('\0'), alt('\0'), qual(""), filter(""),
    info(""), format(""), samples_(nullptr) {};


HaplotypeDataRecord::HaplotypeDataRecord(const std::string& vcf_line, size_t n_cols)
    : chrom(""), pos(-1), id(""),
    ref('\0'), alt('\0'), qual(""), filter(""),
    info(""), format(""), samples_(nullptr) { parse_vcf_line(vcf_line, n_cols); };

// Copy constructor
// Assignment 

void HaplotypeDataRecord::parse_vcf_line(const& vcf_line, size_t n_cols) {

    if (n_cols <= 0)
        throw("No samples detected");

    if (std::isblank(vcf_line[i]))
        throw std::runtime_error("No line can begin with spaces");
    

    // A field of a vcf line record is a single string separated from others
    // by white space.
    // A sample field is a field with the data for a single sample.  Sample
    // fields have numerous : delimited records
    // The counts of the k founders in any one sample field is a comma delimited
    // element of a sample field record.
    int field_idx = 0;              // column index
    int i { 0 };                    // character index of vcf line
    int measurement_idx { 0 };      // index for a measurement in one sample
    int hap_idx { 0 };              // index with hap counts
    int k_founders { 0 };           // number of founders in data set
    int founder_idx { 0 };          // founder index of haplotype dose
    bool hap_found { false };       // determine whether hap dose is in dataset
    int sample_idx { 0 };           // sample index

    buffer_idx_ = 0;                // field buffer index
    for (; i < vcf_line.size(); i++) {

        // A new field is found and stored
        if (i > 0 && std::isblank(vcf_line[i]) && !std::isblank(vcf_line[i-1])) {
            field_idx++;
            buffer_[buffer_idx_] = '\0';
            buffer_idx_ = 0;
            
            if (field_idx == 1)
                chrom = buffer_;
            else if (field_idx == 2)
                pos = std::atoi(buffer_);
            else if (field_idx == 3)
                id = buffer_;
            else if (field_idx == 4)
                ref = buffer_[0];
            else if (field_idx == 5)
                alt = buffer_[0];
            else if (field_idx == 6)
                qual = buffer_;
            else if (field_idx == 7)
                filter = buffer_;
            else if (field_idx == 8)
                info = buffer_;
            else if (field_idx == 9) {
                format = buffer_;

                // verify in the format field that haplotype dose (HD)
                // is included in the data.  Moreover find the index (hap_idx)
                // for which haplotype count data is found in a sample field
                // record
                hap_found = false ;
                for (int i = 0; i < format.size()-1; i++) {
                    if (format[i] == MEASUREMENT_DELIM) {
                        hap_idx++;
                    } else if (format[i] == HAP_CODE[0] 
                                    && format[i+1] == HAP_CODE[1]) {
                        break;
                        hap_found = true;
                    }
                }

                if (!hap_found)
                    throw std::runtime_error("Haplotype counts are not specified");

            } else if ( field_idx > 9 && !samples_ ) {
                // when samples_ is the nullpter it will evaluate to false, 
                // consequently, this block is reserved for building the sample_
                // data matrix.

                // figure out how many founders and instantiate sample matrix
                measurement_idx = 0;
                k_founders = 1;
                sample_idx = 0;

                for (int j = 0; buffer_[j] != '\0'; j++) {
                    if (buffer_[j] == MEASUREMENT_DELIM)
                        measurement_idx++;

                    if (measurement_idx == hap_idx && 
                            buffer_[j] == HAP_DELIM)
                        k_founders++;
                }

                // instantiate sample matrix
                samples_ { k_founders, n_cols - NUM_VCF_FIELDS };
            }

            if (field_idx > 9 && samples_) {

                if (sample_idx < 0 || sample_idx >= n_samps)
                    throw std::runtime_error("Incorrect indexing");

                // find haplotype data
                measurement_idx = 0;
                int j { 0 };            // index of field string
                for (; buffer_[j] != '\0' && j < buffer_size_; j++) {

                    if (buffer_[j] == MEASUREMENT_DELIM)
                        measurement_idx += 1;

                    // found haplotype counts string index
                    if (measurement_idx == hap_idx) {
                        j++
                        break;
                    }
                }

                // decompose haplotype counts to respective founders
                founder_idx = 0;
                hap_buffer_idx_ = 0;
                for (;buffer_[j] != '\0' && j < buffer_idx_; j++) {

                    if (buffer_[j] == HAP_DELIM) {
                        hap_buffer_[hap_buffer_idx_] = '\0';
                        samples_(founder_idx, sample_idx) = std::atod(hap_buffer_);
                        hap_buffer_idx_ = 0;
                        continue;
                    }
                    if (std::isblank(buffer_[j]))
                        continue;

                    hap_buffer_[hap_buffer_idx_++] = buffer_[j];
                }

                sample_idx++;
            }

            continue;
        }

        if (i > 0 && std::isblank(vcf_line[i])
                        && std::isblank(vcf_line[i-1]))
            continue;
        else if (i == 0 && std::isblank(vcf_line[i]))
            throw std::runtime_error("First element of VCF line must not be blank.")

        if (buffer_idx_ >= buffer_size_)
            throw std::runtime_error("buffer overflow, please report");

        buffer_[buffer_idx_++] = vcf_line[i];
    }

}


double HaplotypeDataRecord::operator()(size_t i, size_t j) const {
    return samples_(i, j);
}

