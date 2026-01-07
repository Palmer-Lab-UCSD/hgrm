
#include <cstdio>
#include <string>


#ifndef HEADER_TEXTIO_H
#define HEADER_TEXTIO_H

namespace details {
const size_t DEFAULT_BUF_SIZE = 100;
}


// @title: Parsing text files
// @description: This class manages the lifetime of a C-style file stream
//      by RAII, line retrieval, and getting line unumber.
class TextIO {
public:

    TextIO(const char *filename, const char *mode);
    TextIO(const char *filename, const char *mode, const size_t buf_size);
    ~TextIO();

    int num_lines();
    int get_line();

private:
    std::string fname_;
    const size_t buf_size_;
    char *buf_;
    size_t buf_line_len_ = 0;

    FILE *fid_ = nullptr;

};


TextIO text_open(const char *filename, const char *mode);

#endif
