
#include <cstdio>
#include <string>
#include <memory>


#ifndef HEADER_TEXTIO_H
#define HEADER_TEXTIO_H


namespace io {

enum STATUS { 
    SUCCESS,
    FERROR, 
    FEOF, 
    INVALID_ARG_ERROR, 
    FSEEK_ERROR, 
    FEOF_ERROR,
    END_OF_BUF_ERROR
};


struct FileIO {
    FileIO(FILE *fid): fid(fid) {}; 
    ~FileIO() { if (fid) { fclose(fid); fid = nullptr; } };

    FILE *fid;
};


// @title: file object
// @description: This class manages the lifetime of a C-style file stream
//      by RAII.  To contruct an instance of the class use the "open" function
//      below.
// @param fid: an opened C-style file stream
int bseek(FileIO *fio);

// @title: open a file and instantiate a TextIO object
// @description:
// @param filename: name and path of file to open
// @param mode: a mode in the set of those in the C library function fopen
// @return a unique_ptr<TextIO> if the file stream was successfully opened
//      and TextIO instance created.  Otherwise, return a nullptr.
FileIO *open(const char *filename, const char *mode) {
    if (!mode or !filename)
        return nullptr;

    fid = fopen(filename, mode);
    if (ferror(fid))
        return nullptr;

    if (*mode == 'b')
        return 
}


// @title: File statistics
// @description: This object is returned by any function meant to calculate
//      file character statistics.
struct FileStats {
    size_t nchar = 0;
    size_t nwords = 0;
    size_t nlines = 0;
    size_t nblanklines = 0;
}


// @title: word count
// @description: Similar to the UNIX/Linux wc command line program, wc
//      calculates the number of characters, words, lines, etc. that
//      the specified file contains.
// @param tio: an instance of TextIO
// @param fs: the structure that the file statistics will be stored
// @return a STATUS code that specifies whether the function was successful
//      or failed.
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
