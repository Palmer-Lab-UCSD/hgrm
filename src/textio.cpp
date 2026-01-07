
#include <textio.h>

TextIO::TextIO(const char *filename):
    fname_(filename), 
    buf_size_(details::DEFAULT_BUF_SIZE), 
    buf_(new char[buf_size_]) {
    
    std::memset(buf_, '\0', buf_size_);
};


TextIO::TextIO(const char *filename, const size_t buf_size):
    fname_(filename), 
    buf_size_(buf_size), 
    buf_(new char[buf_size_]) {

    std::memset(buf_, '\0', buf_size_);
};


int TextIO::num_lines() {

    size_t line_num = 0;
    size_t word_len = 0;
    int c;
    while (get_line(fid)) {

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


int get_line(FILE *fid) {
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
