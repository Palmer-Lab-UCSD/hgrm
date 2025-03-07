#ifndef HEADER_UTILS_H
#define HEADER_UTILS_H


#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <string>
#include <functional>
#include <memory>
#include <stdexcept>



class CharBuffer {
public:
    CharBuffer();
    CharBuffer(size_t);
    CharBuffer(const CharBuffer&)=delete;
    CharBuffer(CharBuffer&&)=delete;

    const char& operator()(size_t) const;
    const size_t& size() const;
    const size_t& buffer_size() const;
    
    void append(char s);
    void reset();           // set buffer_idx_ to zero
    void reset(size_t buffer_size);           // set buffer_idx_ to zero

    // remember that it is the caller's responsibility to not
    // dereference raw pointer after the CharBuffer instance is 
    // destructed
    const char* data() const;

private:
    size_t buffer_size_;
    size_t buffer_idx_;
    std::unique_ptr<char[]> buffer_;
};



class StringRecord {
public:
    // input only delimiters
    StringRecord()=delete;
    StringRecord(const char);
    StringRecord(const char, const size_t);
    StringRecord(const char, const char*);

    void update_str(const char* s);

    const char* data() const;
    void reset();
    size_t size();
    bool next_field();


private:
    size_t size_ { 0 };
    size_t idx_ { 0 };
    const char* str_;
    const char delim_;
    
    CharBuffer buf_;
    std::function<bool(char)> is_delim_;
    bool char_is_delim_(char) const;

};


class BufferedRead {
public:
    BufferedRead()=delete;
    BufferedRead(const BufferedRead&)=delete;
    BufferedRead(BufferedRead&&)=delete;
    BufferedRead& operator=(const BufferedRead&)=delete;

    BufferedRead(char* filename, size_t buff_size);
    ~BufferedRead();

    size_t get_line(CharBuffer& line_buf);
    // size_t get_line(std::unique_ptr<char[]> line_buf);
    char get_char();
    void seek(size_t n);
    size_t tell();
    void reset();

private:
    const char* filename_;
    const size_t buff_size_;

    FILE* fid_;
    std::unique_ptr<char[]> buffer_;

    size_t buffer_pos_ { 0 };

    size_t update_buffer_();

};
#endif
