
#include "../include/utils.h"

CharBuffer::CharBuffer(size_t buffer_size)
    : buffer_size_(buffer_size),
        buffer_idx_(0),
        buffer_(buffer_size > 0 ? std::make_unique<char[]>(buffer_size) : nullptr) {

        if (buffer_size_ <= 0)
            throw std::runtime_error("CharBuffer needs to have length > 0");
}


const char& CharBuffer::operator()(size_t idx) const {
    if (idx >= size() || idx >= buffer_idx_)
        throw std::out_of_range("index too large for buffer.");

    if (!buffer_)
        throw std::runtime_error("error");
    return buffer_[idx];
}


const size_t& CharBuffer::size() const { return buffer_size_; }


void CharBuffer::append(char s) {
    if (buffer_idx_ > buffer_size_-2)
        throw std::out_of_range("CharBuffer full");

    buffer_[buffer_idx_++] = s;
    buffer_[buffer_idx_] = '\0';
}


void CharBuffer::reset() {
    buffer_idx_ = 0;
    buffer_[buffer_idx_] = '\0';
}


const char* CharBuffer::data() const {
    return buffer_.get();
}


StringRecord::StringRecord(const std::string delim)
    : str_(nullptr),
        delim_(delim),
        buf_(DEFAULT_BUFFER_SIZE) {};


StringRecord::StringRecord(const std::string delim, const size_t buf_size)
    : str_(nullptr),
        delim_(delim),
        buf_(buf_size) {};


StringRecord::StringRecord(const std::string delim, const char* str)
    : str_(str),
        delim_(delim),
        buf_(DEFAULT_BUFFER_SIZE) {
        
    for (;str_[size_] != '\0'; size_++)
        ;
        
};


// size including null character
void StringRecord::update_str(const char* str) {
    reset();
    str_ = str;

    for (;str_[size_] != '\0'; size_++)
        ;

    return ;
}


void StringRecord::reset() { 
    idx_ = 0;
    size_ = 0;
    str_ = nullptr;
    buf_.reset();
};


size_t StringRecord::size() {
    if (!str_)
        return 0;

    return size_;
}


const char* StringRecord::data() const {
    return buf_.data();
}


bool StringRecord::next_field() {
    if (!str_ || idx_ == size()) 
        return false;

    // eliminate preceeding white space
    for(; idx_ < size(); idx_++)
        if (!is_delim_(str_[idx_]))
            break;

    buf_.reset();
    for (; idx_ < size(); idx_++) {
        if (!is_delim_(str_[idx_]))
            buf_.append(str_[idx_]);

        if (is_delim_(str_[idx_]) && !is_delim_(str_[idx_-1])) {
            idx_++;
            break;
        }
    }

    return true;
}


bool StringRecord::is_delim_(char s) {
    for (size_t i=0; i < delim_.size(); i++) {

        if (s == delim_[i])
            return true;
    }

    return false;
}

