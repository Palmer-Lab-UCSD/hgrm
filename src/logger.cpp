
#include <logger.h>


Logger::Logger(): 
    t_(time(nullptr)),
    time_point_(localtime(&t_)) {

    std::memset(time_buf_, '\0', time_buf_len_); 
    std::memset(str_buf_, '\0', str_buf_len_);
};

int Logger::print_(FILE *stream, 
        const char *log_type, 
        const char *format, 
        const char *msg) {

    // get time and format time string
    t_ = time(nullptr);
    time_point_ = localtime(&t_);


    // strftime returns the number of characters written to buffer,
    // a 0 returned indicates an error has occured.
    time_len_ = strftime(time_buf_, time_buf_len_,"%FT%H:%M:%S", time_point_);

    if(time_len_ == 0) {
        std::memset(time_buf_, '\0', time_buf_len_);

        fprintf(stderr, "%s\t%s\t%s\n", time_buf_, 
            err_str_, "logger time buf failure, please notify maintainer.");
        return -1;
    }

    // construct logging message
    // TODO: truncation of msg notification when msg exceeds buffer
    // Recall that snprintf returns int less than 0 if an error occurs
    msg_len_ = snprintf(str_buf_, max_str_, format, msg);
    if (msg_len_ < 0) {
        fprintf(stderr, "%s\t%s\t%s\n", time_buf_, 
            err_str_, "logger msg failure, please notify maintainer.");
        return -1;
    }

    fprintf(stream, "%s\t%s\t%s\n", time_buf_, log_type, str_buf_);

    return 0;
}

int Logger::info(const char *format, const char *msg) { 
    return print_(stdout, info_str_, format, msg); 
}

int Logger::warn(const char *format, const char *msg) { 
    return print_(stdout, warn_str_, format, msg); 
}

int Logger::error(const char *format, const char *msg) {
    return print_(stderr, err_str_, format, msg);
}

int Logger::info(const char *msg) { 
    return print_(stdout, info_str_, "%s", msg); 
}

int Logger::warn(const char *msg) { 
    return print_(stdout, warn_str_, "%s", msg); 
}

int Logger::error(const char *msg) {
    return print_(stderr, err_str_, "%s", msg);
}

