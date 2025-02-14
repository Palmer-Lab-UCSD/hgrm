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


// Default constructor
HaplotypeDataRecord::HaplotypeDataRecord(size_t n_samples, size_t k_founders) 
    : n_samples_(n_samples),
        k_founders_(k_founders),
        samples_(n_samples_ > 0 && k_founders_ > 0 
                    ? std::make_unique<Matrix>(n_samples_, k_founders_) : nullptr) {

        if (n_samples_ <= 0 || k_founders_ <= 0)
            throw std::runtime_error("Data must have more than zero samples and founders");
         
    };


// access elements
const std::string& HaplotypeDataRecord::chrom() const { return chrom_; };
const long HaplotypeDataRecord::pos() const { return pos_; };
const std::string& HaplotypeDataRecord::id() const { return id_; };
const char HaplotypeDataRecord::ref() const { return ref_; };
const char HaplotypeDataRecord::alt() const { return alt_; };
const std::string& HaplotypeDataRecord::qual() const { return qual_; };
const std::string& HaplotypeDataRecord::filter() const { return filter_; };
const std::string& HaplotypeDataRecord::info() const { return info_; };
const std::string& HaplotypeDataRecord::format() const { return format_; };


void HaplotypeDataRecord::parse_vcf_line(const char* vcf_line) {


    if (std::isspace(vcf_line[0]))
        throw std::runtime_error("No line can begin with spaces");
    

    // A field of a vcf line record is a single string separated from others
    // by white space.
    // A sample field is a field with the data for a single sample.  Sample
    // fields have numerous : delimited records
    // The counts of the k founders in any one sample field is a comma delimited
    // element of a sample field record.
    size_t hap_idx { 0 };              // index with hap counts
    bool hap_found { false };       // determine whether hap dose is in dataset
    size_t sample_idx { 0 };           // sample index
    size_t founder_idx { 0 };

    line_parse_.update_str(vcf_line);

    for (int field_idx = 1; line_parse_.next_field(); field_idx++) {

        if (field_idx == 1)
            chrom_ = line_parse_.data();
        else if (field_idx == 2)
            pos_ = std::atoi(line_parse_.data());
        else if (field_idx == 3)
            id_ = line_parse_.data();
        else if (field_idx == 4)
            ref_ = line_parse_.data()[0];
        else if (field_idx == 5)
            alt_ = line_parse_.data()[0];
        else if (field_idx == 6)
            qual_ = line_parse_.data();
        else if (field_idx == 7)
            filter_ = line_parse_.data();
        else if (field_idx == 8)
            info_ = line_parse_.data();
        else if (field_idx == 9) {
            format_ = line_parse_.data();

            // verify in the format field that haplotype dose (HD)
            // is included in the data.  Find the index (hap_idx)
            // for which haplotype count data is found in a sample field
            // record
            field_parse_.update_str(line_parse_.data());

            hap_found = false;
            for (size_t i = 0; field_parse_.next_field(); i++) {
                if (std::strcmp(field_parse_.data(), HAP_CODE) == 0 ) {
                    hap_found = true;
                    hap_idx = i;
                    break;
                }
                hap_idx++;
            }

            if (!hap_found)
                throw std::runtime_error("Haplotype counts are not specified");

        } else if (field_idx > 9 && samples_) {

            if (sample_idx < 0 || sample_idx >= n_samples_)
                throw std::out_of_range("Index is out of matrix range.");

            field_parse_.update_str(line_parse_.data());
            for (size_t j = 0; field_parse_.next_field(); j++)
                if (j == hap_idx)
                    break;

            hap_parse_.update_str(field_parse_.data());

            // decompose haplotype counts to respective founders
            for (founder_idx = 0; hap_parse_.next_field(); founder_idx++)
                (*samples_)(sample_idx, founder_idx) = std::atof(hap_parse_.data());


            if (founder_idx != k_founders_)
                throw std::runtime_error("Number of founders found for sample is incorrect");
            
            sample_idx++;
        }
    }

    if (sample_idx != n_samples_)
        throw std::runtime_error("Number of samples found is not equal to that expected.");

}


const double& HaplotypeDataRecord::operator()(size_t i, size_t j) const {
    return (*samples_)(i, j);
}

std::array<size_t, 2> HaplotypeDataRecord::dims() const {
    return (*samples_).dims();
}
