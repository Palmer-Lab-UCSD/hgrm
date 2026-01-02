
#include <logger.h>


Logger::Logger(): 
    t_(time(nullptr)),
    time_point_(localtime(&t_)) { empty_bufs_(); };


void Logger::empty_bufs_() {
    std::memset(time_buf_, '\0', time_buf_len_);
    std::memset(str_buf_, '\0', str_buf_len_);
}


int Logger::load_time_buf_() {
    // get time and format time string
    t_ = time(nullptr);
    time_point_ = localtime(&t_);

    // strftime returns the number of characters written to buffer,
    // a 0 returned indicates an error has occured.
    return strftime(time_buf_, time_buf_len_,"%FT%H:%M:%S", time_point_);
}


void Logger::vprintf_(FILE *stream, const char *log_type, const char *format, 
        va_list arg_ptr) {

    if ((status_ = load_time_buf_()) <= 0)
        strncpy(str_buf_, 
                "logger time buffer failure, please notify maintainer.",
                max_str_);
    else if ((status_ = vsnprintf(str_buf_, max_str_, format, arg_ptr)) <= 0)
        strncpy(str_buf_,
                "logger msg failure, please notify maintainer.", 
                max_str_);
    else if (status_ >= max_str_) {
        status_ = -status_;
        strncpy(str_buf_,
                "logger msg too long, please shorten msg.", 
                max_str_);
    }

    if (status_ <= 0)
        fprintf(stderr, "%s\t%s\t%s\n", time_buf_, err_str_, str_buf_);
    else
        fprintf(stream, "%s\t%s\t%s\n", time_buf_, log_type, str_buf_);

    empty_bufs_();
}


int Logger::info(const char *format, ...) { 
    va_list arg_ptr;
    va_start(arg_ptr, format);
    vprintf_(stdout, info_str_, format, arg_ptr);
    va_end(arg_ptr);
    return status_;
}

int Logger::warn(const char *format, ...) { 
    va_list arg_ptr;
    va_start(arg_ptr, format);
    vprintf_(stdout, warn_str_, format, arg_ptr);
    va_end(arg_ptr);
    return status_;
}

int Logger::error(const char *format, ...) {
    va_list arg_ptr;
    va_start(arg_ptr, format);
    vprintf_(stderr, err_str_, format, arg_ptr);
    va_end(arg_ptr);
    return status_;
}

// int Logger::info(const char *msg) { 
//     return info("%s", msg);
// }
// 
// int Logger::warn(const char *msg) { 
//     return warn("%s", msg);
// }
// 
// int Logger::error(const char *msg) {
//     return error("%s", msg);
// }

