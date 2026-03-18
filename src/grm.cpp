// Palmer Lab at UCSD
//
//
// ACKNOWLEDGMENT
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Opus 4.6, the AI assistant from Anthropic
// with minor recommendations incorporated.
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
// COORDINATES CLASS
////////////////////////////////////////////////////////////////////

grm::Coordinates::Coordinates(Coordinates&& other)
    : len(0), contig(""), pos(nullptr) {
    len = other.len;
    contig = other.contig;
    pos = std::move(other.pos);

    other.pos=nullptr;
    other.len = 0;
    other.contig = "";
}

grm::Coordinates& grm::Coordinates::operator=(Coordinates&& other) {
    if (this == &other)
       return *this; 

    len = other.len;
    contig = other.contig;
    pos = std::move(other.pos);

    other.pos=nullptr;
    other.len = 0;
    other.contig = "";
}

// remember that Coordinates* should be uninstantiated
grm::STATUS write(io::FileIO* fio, const Coordinates* coords) {

    if (!fio)
        return grm::ERROR_NULLPTR_ARG;

    if (!fio->fid)
        return grm::ERROR_NULLPTR_ARG;

    if (!coords)
        return grm::ERROR_NULLPTR_ARG;

    size_t nwritten = 0;

    // write contig name to file
    size_t nchar = coords->contig.size();
    nwritten = fwrite(&nchar, sizeof(nchar), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    nwritten = fwrite(coords->contig.c_str(), 
            sizeof(char), 
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
            sizeof(size_t),
            npos,
            fio->fid);
    if (nwritten != npos)
        return grm::ERROR_ON_WRITE;

    return grm::SUCCESS;
}


grm::STATUS read(io::FileIO* fio, Coordinates* coords) {
    if (!fio)
        return grm::ERROR_NULLPTR_ARG;

    if (!fio->fid)
        return grm::ERROR_NULLPTR_ARG;

    if (!coords)
        return grm::ERROR_NULLPTR_ARG;

    // I create a temporary Coordinates class, because I don't want
    // the input coords instance to partial update upon an error
    grm::Coordinates tmpc {};

    // will store the number of bytes read at each step
    size_t nread = 0;

    // read contig name
    size_t size_contig_name = 0;
    nread = fread(&size_contig_name, sizeof(size_t), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;

    std::unique_ptr<char[]> buffer = std::make_unique<char[]>(size_contig_name + 1);
    std::memset(buffer.get(), '\0', size_contig_name + 1);

    nread = fread(buffer, sizeof(char), size_contig_name, fio->fid);
    if (nread != size_contig_name)
        return grm::ERROR_ON_READ;

    tmpc.contig = std::string(buffer);

    // read in positions
    size_t npos = 0;
    nread = fread(&npos, sizeof(size_t), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;

    tmpc.len = npos;

    tmpc.pos = std::make_unique<size_t[]>(npos);
    nread = fread(tmpc.pos.get(), sizeof(size_t), npos, fio->fid);    
    if (nread != npos)
        return grm::ERROR_ON_READ;

    *coords = std::move(tmpc);

    return grm::SUCCESS;
}

////////////////////////////////////////////////////////////////////
// SAMPLES CLASS
////////////////////////////////////////////////////////////////////

grm::Samples::Samples(grm::Samples&& other) 
    : len(other.len), names(nullptr) {
        names = std::move(other.names);
        other.len = 0;
        other.names = nullptr;
}


grp::Samples& grm::Samples::operator=(grm::Samples&& other) {
    if (this == &other)
        return *this;
    
    len = other.len;
    names = std::move(other.names);
    
    other.len = 0;
    other.names = nullptr;
    
    return *this;
}


grm::STATUS write(io::FileIO* fio, const grm::Samples* samples) {
    if (!fio)
        return grm::ERROR_NULLPTR_ARG;

    if (!fio->fid)
        return grm::ERROR_NULLPTR_ARG;

    if (!samples)
        return grm::ERROR_NULLPTR_ARG;

    size_t nwritten = 0;
    size_t nsamps = samples->len;
    nwritten = fwrite(&nsamps, sizeof(size_t), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    // When it comes time to read the data, I need to make a character
    // buffer to temporarily place the read string.  To make this
    // buffer, I need to know the length of string with the greatest
    // number of characters.  Here I find that number and store in
    // the binary file.
    
    size_t nchar_max = 0;
    size_t tmp = 0;
    for (size_t n = 0; n < nsamps; n++)
        if ((tmp = samples->name[n].size()) > nchar_max) nchar_max = tmp;

    if (nchar_max == 0)
        return grm::ERROR_INVALID_ARG;

    nwritten = fwrite(&nchar_max, sizeof(size_t), 1, fio->fid);
    if (nwritten != 1)
        return grm::ERROR_ON_WRITE;

    // Write each string to file;
    size_t nchar = 0;
    for (size_t n = 0; n < nsamps; n++) {
        nchar = samples->names[n].size(); 

        nwritten = fwrite(&nchar, sizeof(size_t), 1, fio->fid);
        if (nwritten != 1)
            return grm::ERROR_ON_WRITE;

        nwritten = fwrite(samples->names[n].c_str(),
                sizeof(char), 
                nchar,
                fio->fid);
        if (nwritten != nchar)
            return grm::ERROR_ON_WRITE;
    }

    return grm::SUCCESS;
}

grm::STATUS read(io::FileIO* fio, grm::Samples* samples) {
    if (!fio)
        return grm::ERROR_NULLPTR_ARG;

    if (!fio->fid)
        return grm::ERROR_NULLPTR_ARG;

    // I create a temporary Sample class, because I don't want the 
    // input samples instance to partial update upon an error
    grm::Samples tmp_samps();

    size_t nread;
    size_t n_samples = 0;

    nread = fread(&n_samples, sizeof(size_t), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;
    
    tmp_samps->len = n_samples;

    // Get the number of characters of the longest string
    size_t nchar_max = 0;
    nread = fread(&nchar_max, sizeof(size_t), 1, fio->fid);
    if (nread != 1)
        return grm::ERROR_ON_READ;


    std::unique_ptr<char[]> buffer = std::make_unique<char[]>(nchar_max + 1);
    std::memset(buffer.get(), '\0', nchar_max + 1);

    size_t nchar = 0;
    for (size_t n = 0; n < n_samples: n++) {
        nread = fread(&nchar, sizeof(size_t), 1, fio->fid);
        if (nread != 1)
            return grm::ERROR_ON_READ;
        
        nread = fread(buffer.get(), sizeof(char), nchar, fio->fid);
        if (nread != nchar)
            return grm::ERROR_ON_READ;

        tmp_samps->names[n] = std::string(buffer);

        nchar = 0;
        std::memset(buffer.get(), '\0', nchar);
    }

    *samples = std::move(tmp_samps);

    return grm::SUCCESS;
}

////////////////////////////////////////////////////////////////////
// HDR CLASS
////////////////////////////////////////////////////////////////////


grm::STATUS write(io::FileIO* fio, const Hdr* hdr) {
    if (!fio)
        return grm::ERROR_NULLPTR_ARG;

    if (!fio->fid)
        return grm::ERROR_NULLPTR_ARG;

    if (!hdr)
        return grm::ERROR_NULLPTR_ARG;

    fwrite(&hdr->version.size(), sizeof(size_t), 1, fio->fid);
    fwrite(hdr->version.c_str(), sizeof(char), hdr->version.size(), fio->fid);

    fwrite(&hdr->grm_type, sizeof(GrmType), 1, fio->fid);
    grm::STATUS status;
    if ((status = write(fio, hdr->coords)) != grm::SUCCESS)
        return status;

    if ((status = write(fio, hdr->samples)) != grm::SUCCESS)
        return status;

    return grm::SUCCESS;
}

grm::STATUS read(io::FileIO* fio, Hdr* hdr) {
    if (!fio)
        return grm::ERROR_NULLPTR_ARG;

    if (!fio->fid)
        return grm::ERROR_NULLPTR_ARG;

    if (!hdr)
        return grm::ERROR_NULLPTR_ARG;

}
////////////////////////////////////////////////////////////////////
// GRM CLASS
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
