#include "../include/utils.h"

const static size_t DEFAULT_BUFFER_SIZE { 1000 };




CharBuffer::CharBuffer()
    : buffer_size_(0),
        buffer_idx_(0),
        buffer_(nullptr) {}


CharBuffer::CharBuffer(size_t buffer_size)
    : buffer_size_(buffer_size),
        buffer_idx_(0),
        buffer_(buffer_size > 0 ? std::make_unique<char[]>(buffer_size+1) : nullptr) {

        if (buffer_ == nullptr)
            throw std::runtime_error("CharBuffer needs to have length > 0");

        buffer_[buffer_idx_] = '\0';
        buffer_[buffer_size_] = '\0';
}


const char& CharBuffer::operator()(size_t idx) const {
    if (idx >= size() || idx >= buffer_idx_)
        throw std::out_of_range("index too large for buffer.");

    if (!buffer_)
        throw std::runtime_error("error");
    return buffer_[idx];
}


const size_t& CharBuffer::buffer_size() const { return buffer_size_; }
const size_t& CharBuffer::size() const { return buffer_idx_; }


void CharBuffer::append(char s) {
    if (buffer_ == nullptr)
        throw std::runtime_error("No data stored in buffer.");

    if (buffer_idx_ >= buffer_size_)
        throw std::out_of_range("CharBuffer full");

    buffer_[buffer_idx_++] = s;
    buffer_[buffer_idx_] = '\0';
}


void CharBuffer::reset(size_t buffer_size) {
    buffer_size_ = buffer_size;
    buffer_ = std::make_unique<char[]>(buffer_size_);
    buffer_[0] = '\0';
    buffer_[buffer_size_] = '\0';
    buffer_idx_ = 0;
}


void CharBuffer::reset() {
    if (buffer_ == nullptr)
        throw std::runtime_error("No data stored in buffer.");

    buffer_idx_ = 0;
    buffer_[buffer_idx_] = '\0';
}


const char* CharBuffer::data() const {
    if (buffer_ == nullptr)
        throw std::runtime_error("No data stored in buffer.");

    return buffer_.get();
}



StringRecord::StringRecord(const char delim)
    : str_(nullptr),
        delim_(delim),
        buf_(DEFAULT_BUFFER_SIZE) { 

    if (std::isspace(delim))
        is_delim_ = [](char c){ return std::isspace(c); };
    else
        is_delim_ = [this](char c){ return char_is_delim_(c);  };
};

StringRecord::StringRecord(const char delim, const size_t buf_size)
    : str_(nullptr),
        delim_(delim),
        buf_(buf_size) { 

    if (std::isspace(delim))
        is_delim_ = [](char c){ return std::isspace(c); };
    else
        is_delim_ = [this](char c){ return char_is_delim_(c); }; 
};

StringRecord::StringRecord(const char delim, const char* str)
    : str_(str),
        delim_(delim),
        buf_(DEFAULT_BUFFER_SIZE) {

    if (std::isspace(delim))
        is_delim_ = [](char c){ return std::isspace(c); };
    else
        is_delim_ = [this](char c){ return char_is_delim_(c); };

    for (;str_[size_] != '\0'; size_++)
        ;
};


bool StringRecord::char_is_delim_(char c) const {
    if (c == delim_)
        return true;
    return false;
}

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


// bool StringRecord::is_delim_(char s) {
//     return s == delim_;
// }


// buff_size is the count of characters to load
BufferedRead::BufferedRead(char* filename, size_t buff_size)
    : filename_(filename), 
        buff_size_(buff_size),
        fid_(std::fopen(filename, "r")),
        buffer_(buff_size > 0 ? std::make_unique<char[]>(buff_size+1) : nullptr) {

    if (!fid_)
        throw std::runtime_error("File Access error");

    if (std::feof(fid_))
        throw std::runtime_error("File is empty");

    if (buffer_ == nullptr)
        throw std::runtime_error("Buffer wasn't properly set");

    buffer_pos_ = 0;
    buffer_[0] = '\0';
}


BufferedRead::~BufferedRead() {
    if (fid_)
        fclose(fid_);
}

// Note:
// Return count == 0 indicates that we have reached the end of file.
//
size_t BufferedRead::get_line(CharBuffer& line_buffer) {

    line_buffer.reset();

    size_t count { 0 };
    size_t buff_length { 1 }; 

    if (buffer_[buffer_pos_] == '\0'
            || buffer_pos_ == buff_size_)
        buff_length = update_buffer_();


    if (buff_length == 0)
        return count;

    while (buffer_[buffer_pos_] != '\n') {

        line_buffer.append(buffer_[buffer_pos_++]);
        count++;

        if (buffer_[buffer_pos_] == '\0' || buffer_pos_ == buff_size_)
            buff_length = update_buffer_();

        if (buff_length == 0)
            break;
    }

    if (buffer_[buffer_pos_] == '\n')
        buffer_pos_++;

    return count;
}


char BufferedRead::get_char() {
    size_t buff_length { 1 };


    if (buffer_[buffer_pos_] == '\0')
        buff_length = update_buffer_();

    // when all entries of a file are read, detected by
    // buff_length == 0, return the null character
    if (buff_length == 0)
        return '\0';

    return buffer_[buffer_pos_++];    
}


size_t BufferedRead::update_buffer_() {

    if (std::ferror(fid_))
        throw std::runtime_error("File read error");

    // note, if end of file, then we set n = 0;
    size_t n { 0 };
    if (!std::feof(fid_)) 
        n = fread(buffer_.get(), sizeof(buffer_[0]), buff_size_, fid_);

    buffer_pos_ = 0;
    buffer_[n] = '\0';

    return n;
}


void BufferedRead::seek(size_t n) {

    int c = 0;
    if ((c = std::fseek(fid_, n, SEEK_SET)) != 0)
        throw std::runtime_error("Failed to relocate file stream to position.");

}


size_t BufferedRead::tell() {
    return std::ftell(fid_);
}


void BufferedRead::reset() {
    buffer_pos_ = 0;
    buffer_[buffer_pos_] = '\0';
    seek(0);
}


