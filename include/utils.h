#ifndef HEADER_UTILS_H
#define HEADER_UTILS_H


#include <cstdlib>
#include <cstring>
#include <string>
#include <memory>
#include <stdexcept>


const size_t DEFAULT_BUFFER_SIZE { 1000 };

class CharBuffer {
public:
    CharBuffer()=delete;
    CharBuffer(size_t);
    CharBuffer(const CharBuffer&)=delete;
    CharBuffer(CharBuffer&&)=delete;

    const char& operator()(size_t) const;
    const size_t& size() const;
    
    void append(char s);
    void reset();           // set buffer_idx_ to zero

    // remember that it is the caller's responsibility to not
    // dereference raw pointer after the CharBuffer instance is 
    // destructed
    const char* data() const;

private:
    const size_t buffer_size_;
    size_t buffer_idx_;
    std::unique_ptr<char[]> buffer_;
};



class StringRecord {
public:
    // input only delimiters
    StringRecord()=delete;
    StringRecord(const std::string);
    StringRecord(const std::string, const size_t);
    StringRecord(const std::string, const char*);

    void update_str(const char* s);

    const char* data() const;
    void reset();
    size_t size();
    bool next_field();


private:
    size_t size_ { 0 };
    size_t idx_ { 0 };
    const char* str_;
    const std::string delim_;
    
    CharBuffer buf_;
    bool is_delim_(char);
};


#endif
