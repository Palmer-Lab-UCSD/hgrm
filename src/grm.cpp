// MAtrix
//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-10
//
//
// Acknowledgment
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Sonnet, the AI assistant from Anthropic
// (Jan 2025), with minor recommendations incorporated.
//
//

#include <grm.h>

int grm::details::num_lines_in_file(FILE *fid, size_t *num_lines) {
    // TODO: errno, need to reset?

    size_t line_num = 0;
    size_t word_len = 0;
    int c;
    while ((c = fgetc(fid)) != EOF) {

        if (c == '\n' && word_len != 0) {
            line_num++;
            word_len = 0;
        } else if (c != '\n')
            word_len++;
    }

    if (ferror(fid))
        return fseek(fid, 0, SEEK_SET) == 0 ? -1 : -3;

    if (feof(fid) == 0)
        return fseek(fid, 0, SEEK_SET) == 0 ? -2 : -3;

    *num_lines = line_num;
    return fseek(fid, 0, SEEK_SET) == 0 ? 0 : -3;
}


// int grm::details::get_size_t(FILE *fid, size_t *val) {
// 
//     std::string s { "" };
//     while (std::getline(fid, s))
//         if (s.size() < )
// 
//     int c;
//     for (int i = 0; (c = fgetc(fid)) != EOF && i < max_bitsize_size_t; i++) {
//         if (c == '\n')
//             break;
//         s[i] = c; 
//     }
//     s[i] = '\0';
// 
//     return 0;
// }
// 

grm::Coordinates::Coordinates(const char *contig, const size_t len):
    contig(contig), len(len) {}



// default constructor
grm::Grm::Grm(const size_t nrow, const size_t mcol)
    : dims_(nrow, mcol), 
        data_(size() != 0 ? std::make_unique<float[]>(size()) : nullptr) {

    if (data_)
       std::memset(data_.get(), 0, size());
}


// copy constructor
//
grm::Grm::Grm(const grm::Grm& other) 
    : nrow_(dims.other.nrow_), mcol_(dims.other.mcol_),
        data_(std::make_unique<float[]>(other.size())) {
    std::memset(data_.get(), 0, size());
}


// TODO: check this.
grm::Grm::Grm(grm::Grm&& other) 
    : nrow_(dims_.other.nrow_), dims.mcol_(other.mcol_), 
        data_(std::move(other.data_)) {};


float grm::Grm::operator()(const size_t i, const size_t j) const {
    return data_[i*mcol_ + j];
}


float& grm::Grm::operator()(const size_t i, const size_t j) {
    return data_[i*mcol_ + j];
}


grm::STATUS grm::Grm::get(const size_t i, const size_t j, float *val) const {
    size_t idx = 0;
    grm::STATUS status = grm::STATUS::FAILED;
    if ((status = midx_to_arr_(i, j, &idx)) != grm::STATUS::SUCCESS) 
        return status;

    *val = data_[idx];

    return status;
}


grm::STATUS grm::Grm::set(const size_t i, const size_t j, const float val) {
    size_t idx = 0;
    grm::STATUS status = grm::STATUS::FAILED;
    if ((status = midx_to_arr_(i, j, &idx)) != grm::STATUS::SUCCESS) 
        return status;

    data_[idx] = val;

    return status;
}


const grm::Dims& grm::Grm::dims() const { return dims_; };


grm::STATUS grm::Grm::midx_to_arr_(const size_t i, const size_t j, size_t *idx) const {

    if (i >= nrow_ || j >= mcol_)
        return grm::STATUS::ERROR_IDX_ARR_BOUNDS;

    *idx = i*mcol_ + j;
    return grm::STATUS::SUCCESS;
}


size_t grm::Grm::size() const { return dims_.nrow_ * dims_.mcol_; };


grm::STATUS grm::Grm::write(const char *filename) const {

    std::unique_ptr<FILE> fid = make_unique<FILE>(fopen(filename, "wb"));

    size_t size_written = fwrite(&dims_, sizeof(Dims), 1, fid.get());
    if (size_written < 1)
        return grm::STATUS::ERROR_ON_WRITE;

    size_written = fwrite(data_.get(), 
            sizeof(float), 
            size(),
            fid.get());

    if (size_written < size())
        return grm::STATUS::ERROR_ON_WRITE;

    return grm::STATUS::SUCCESS;
}


grm::STATUS grm::Grm::read(const char *filename, grm::Grm *grm) {
    return grm::STATUS;
}
