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
HaplotypeDataRecord::HaplotypeDataRecord()
    : chrom_(""), pos_(-1), id_(""),
    ref_('\0'), alt_('\0'), qual_(""), filter_(""),
    info_(""), format_(""), samples_(nullptr) {};


// Constructor
HaplotypeDataRecord::HaplotypeDataRecord(const std::string& vcf_line, size_t n_cols)
    : chrom_(""), pos_(-1), id_(""),
    ref_('\0'), alt_('\0'), qual_(""), filter_(""),
    info_(""), format_(""), samples_(nullptr) { parse_vcf_line(vcf_line, n_cols); };


// Copy constructor
// remember the `this` object is being initialized, it doesn't already exist
HaplotypeDataRecord::HaplotypeDataRecord(const HaplotypeDataRecord& other)
    : chrom_(other.chrom_), pos_(other.pos_), id_(other.id_),
    ref_(other.ref_), alt_(other.alt_), qual_(other.qual_), filter_(other.filter_),
    info_(other.info_), format_(other.format_), 
    samples_(std::make_unique<Matrix>(*(other.samples_))) {};


// Copy assignment operator
HaplotypeDataRecord& HaplotypeDataRecord::operator=(const HaplotypeDataRecord& other) {
    if (this == &other)
        return *this;

    chrom_ = other.chrom_;
    pos_ = other.pos_;
    id_ = other.id_;
    ref_ = other.ref_;
    alt_ = other.alt_;
    qual_ = other.qual_;
    filter_ = other.filter_;
    info_ = other.info_;
    format_ = other.format_;

    // I assume that make_unique will call the copy constructor
    // of Matrix, which will request the memory accordingly.  Moreover
    // I assume that the copy assignment operator of unique_ptr
    // releases the original memory.
    samples_ = std::make_unique<Matrix>(*(other.samples_));

    return *this;
}


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


void HaplotypeDataRecord::parse_vcf_line(const std::string& vcf_line, size_t n_cols) {

    if (n_cols == 0)
        throw("No samples detected");

    if (std::isblank(vcf_line[0]))
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
    size_t k_founders { 0 };           // number of founders in data set
    int founder_idx { 0 };          // founder index of haplotype dose
    bool hap_found { false };       // determine whether hap dose is in dataset
    int sample_idx { 0 };           // sample index
    size_t n_samples { 0 };

    buffer_idx_ = 0;                // field buffer index
    for (; i < vcf_line.size(); i++) {

        // A new field is found and stored
        if (i > 0 && std::isspace(vcf_line[i]) && !std::isspace(vcf_line[i-1])) {
            field_idx++;
            if (buffer_idx_ >= buffer_size_)
                throw std::runtime_error("Buffer overflow, field exceeds buffer len");

            buffer_[buffer_idx_] = '\0';
            buffer_idx_ = 0;

            // Note that buffer_ is being used as a safer and size informed
            // C-array.  To convert a c-style string to std::string we need to
            // use array's member function data.
            
            if (field_idx == 1)
                chrom_ = buffer_.data();
            else if (field_idx == 2)
                pos_ = std::atoi(buffer_.data());
            else if (field_idx == 3)
                id_ = buffer_.data();
            else if (field_idx == 4)
                ref_ = buffer_[0];
            else if (field_idx == 5)
                alt_ = buffer_[0];
            else if (field_idx == 6)
                qual_ = buffer_.data();
            else if (field_idx == 7)
                filter_ = buffer_.data();
            else if (field_idx == 8)
                info_ = buffer_.data();
            else if (field_idx == 9) {
                format_ = buffer_.data();

                //std::cout << vcf_line << std::endl;
                // verify in the format field that haplotype dose (HD)
                // is included in the data.  Moreover find the index (hap_idx)
                // for which haplotype count data is found in a sample field
                // record
                hap_found = false ;
                for (int i = 0; i < format_.size()-1; i++) {

                    if (format_[i] == MEASUREMENT_DELIM) {
                        hap_idx++;
                    } else if (format_[i] == HAP_CODE[0] 
                                    && format_[i+1] == HAP_CODE[1]) {
                        hap_found = true;
                        break;
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

                for (int j = 0; buffer_[j] != '\0'; j++) {
                    if (j == buffer_size_)
                        throw std::runtime_error("Buffer at capacity error");

                    if (buffer_[j] == MEASUREMENT_DELIM)
                        measurement_idx++;

                    if (measurement_idx == hap_idx && buffer_[j] == HAP_DELIM)
                        k_founders++;
                }


                // instantiate sample matrix
                n_samples = n_cols - NUM_VCF_FIELDS;
                samples_ = std::make_unique<Matrix>(k_founders, n_samples);
            }


            if (field_idx > 9 && samples_) {
                std::cout << sample_idx << " of " << n_samples << std::endl;

                if (sample_idx < 0 || sample_idx >= n_samples)
                    throw std::runtime_error("Incorrect indexing");

                // find haplotype data
                measurement_idx = 0;
                int j { 0 };            // index of field string
                for (; buffer_[j] != '\0' && j < buffer_size_; j++) {

                    if (buffer_[j] == MEASUREMENT_DELIM)
                        measurement_idx += 1;

                    // found haplotype counts string index
                    if (measurement_idx == hap_idx) {
                        j++;
                        break;
                    }
                }

                // decompose haplotype counts to respective founders
                founder_idx = 0;
                hap_buffer_idx_ = 0;
                for (;buffer_[j] != '\0' && j < buffer_size_; j++) {

                    if (buffer_[j] == HAP_DELIM) {
                        hap_buffer_[hap_buffer_idx_] = '\0';
                        (*samples_)(founder_idx, sample_idx) = std::atof(hap_buffer_.data());
                        hap_buffer_idx_ = 0;
                        founder_idx++;
                        continue;
                    }

                    if (std::isspace(buffer_[j]))
                        continue;

                    hap_buffer_[hap_buffer_idx_++] = buffer_[j];
                }

                // Note that the last haplotype's data terminates with \0, and
                // not the HAP_DELIM
                if (buffer_[j] == '\0') {
                    hap_buffer_[hap_buffer_idx_] = '\0';
                    (*samples_)(founder_idx, sample_idx) = std::atof(hap_buffer_.data());
                    hap_buffer_idx_ = 0;
                }

                sample_idx++;
            }

            continue;
        }

        if (i > 0 && std::isspace(vcf_line[i])
                        && std::isspace(vcf_line[i-1]))
            continue;
        else if (i == 0 && std::isspace(vcf_line[i]))
            throw std::runtime_error("First element of VCF line must not be blank.");

        if (buffer_idx_ >= buffer_size_)
            throw std::runtime_error("buffer overflow, please report");

        buffer_[buffer_idx_++] = vcf_line[i];
    }

}


double HaplotypeDataRecord::operator()(size_t i, size_t j) const {
    return (*samples_)(i, j);
}

std::array<size_t, 2> HaplotypeDataRecord::dims() const {
    return (*samples_).dims();
}
