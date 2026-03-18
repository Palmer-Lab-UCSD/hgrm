// Palmer Lab at UCSD
//
//
// ACKNOWLEDGMENT
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Sonnet, the AI assistant from Anthropic
// (Jan 2025), with minor recommendations incorporated.
//
//

#include <grm.h>

// int grm::details::num_lines_in_file(FILE *fid, size_t *num_lines) {
//     // TODO: errno, need to reset?
// 
//     size_t line_num = 0;
//     size_t word_len = 0;
//     int c;
//     while ((c = fgetc(fid)) != EOF) {
// 
//         if (c == '\n' && word_len != 0) {
//             line_num++;
//             word_len = 0;
//         } else if (c != '\n')
//             word_len++;
//     }
// 
//     if (ferror(fid))
//         return fseek(fid, 0, SEEK_SET) == 0 ? -1 : -3;
// 
//     if (feof(fid) == 0)
//         return fseek(fid, 0, SEEK_SET) == 0 ? -2 : -3;
// 
//     *num_lines = line_num;
//     return fseek(fid, 0, SEEK_SET) == 0 ? 0 : -3;
// }


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

////////////////////////////////////////////////////////////////////
// Coordinates class
////////////////////////////////////////////////////////////////////

grm::Coordinates::Coordinates(const char *contig, const size_t len):
    contig(contig), len(len), pos(std::make_unique<size_t>(len)) {};


grm::STATUS write(io::FileIO *fio, const Coordinates *coords) {

    if (!fio->fid)
        return grm::ERROR_FILE_NOT_OPEN;

    // write contig name to file
    size_t nwritten = 0;
    size_t nchar = coords->contig.size();
    nwritten = fwrite(&nchar, sizeof(nchar), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    nwritten = fwrite(coords->contig.c_str(), 
            sizeof(coords->name[0]), 
            nchar,
            fio->fid);
    if (nwritten != nchar)
        return grm::ERROR_ON_WRITE;

    // write positions
    size_t npos = coords->len;
    nwritten = fwrite(&npos, sizeof(npos), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    nwritten = fwrite(coords->pos.get(), 
            sizeof(coords->pos[0]),
            npos,
            fio->fid);
    if (nwritten != pos)
        return grm::ERROR_ON_WRITE;

    return grm::SUCCESS;
}


grm::STATUS read(io::FileIO* fio, Coordinates* coords) {
}


////////////////////////////////////////////////////////////////////
// GRM class
////////////////////////////////////////////////////////////////////
//
// Recall that the GRM is a symmetric matrix, therefore we only need
// to store the upper triagonal and diagonal element values.  
// Consequently, the size of the array storing the data is n*(n +1)/2.
//
grm::Grm::Grm(const size_t n_samples)
    : n_samples_(n_samples)
        data_(size() != 0 ? std::make_unique<float[]>(size()) : nullptr) {

    if (data_)
       std::memset(data_.get(), 0, size());
}

// The number of upper diagonal + diagonal elements of the GRM
size_t grm::Grm::size() const { return n_samples * (n_samples + 1) / 2; };


grm::STATUS grm::Grm::midx_to_arr_(const size_t i, const size_t j, size_t *idx) const {

    if (i >= n_samples_ || j >= n_samples_)
        return grm::ERROR_IDX_ARR_BOUNDS;

    // remember that by symmetry, the matrix is equal to its transpose
    if (i <= j)
        *idx = MATRIX_IDX_TO_ARRAY(i, j, n_samples_);
    else 
        *idx = MATRIX_IDX_TO_ARRAY(j, i, n_samples_);

    return grm::SUCCESS;
}


float grm::Grm::operator()(const size_t i, const size_t j) const {
    return data_[MATRIX_IDX_TO_ARRAY(i, j, n_samples_)];
}


float& grm::Grm::operator()(const size_t i, const size_t j) {
    return data_[MATRIX_IDX_TO_ARRAY(i, j, n_samples_)];
}


grm::STATUS grm::Grm::get(const size_t i, const size_t j, float *val) const {
    size_t idx = 0;
    grm::STATUS status = grm::FAILED;
    if ((status = midx_to_arr_(i, j, &idx)) != grm::SUCCESS) 
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


grm::STATUS grm::Grm::write(io.FileIO *fio, const Hdr *hdr) const {

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


grm::Grm grm::Grm::read(io.FileIO *fio) {
    return grm::STATUS;
}
