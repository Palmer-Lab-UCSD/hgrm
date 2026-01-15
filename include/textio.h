
#include <cstdio>
#include <string>
#include <memory>


#ifndef HEADER_TEXTIO_H
#define HEADER_TEXTIO_H


namespace textio {

enum STATUS { 
    SUCCESS,
    FERROR, 
    FEOF, 
    INVALID_ARG_ERROR, 
    FSEEK_ERROR, 
    FEOF_ERROR,
    END_OF_BUF_ERROR
};


// @title: Parsing text files
// @description: This class manages the lifetime of a C-style file stream
//      by RAII, line retrieval, and getting line unumber.
struct TextIO {
    TextIO(FILE *fid);
    ~TextIO();

    int bseek();

    FILE *fid;
};


std::unique_ptr<TextIO> open(const char *filename, const char *mode);


struct FileStats {
    size_t nchar = 0;
    size_t nwords = 0;
    size_t nlines = 0;
    size_t nblanklines = 0;
}

STATUS wc(TextIO *tio, FileStats *fs);


template<typename T>
struct Array {
    Array(size_t size_in): size(size_in), 
        data(size > 0 ? new T[size] : nullptr) {};

    ~Array() { if (data) delete[] data; };

    size_t size;
    T *data;
    size_t len = 0;

    //unsafe referencing
    T operator[](size_t i) { return data[i]; };
    T& operator[](size_t i) { return data[i]; };

    STATUS append(T val) {
        if (len >= size-1)
            return END_OF_BUF_ERROR;

        data[len++] = val;
        return SUCCESS;
    }

    void fill(T val) {
        std::memset(data, val, size);
        len = 0;
    }
}


STATUS getline(TextIO *tio, Array<char> linebuf);
}

#endif
