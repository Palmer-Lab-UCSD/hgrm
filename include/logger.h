
#ifndef HEADER_LOGGER_H
#define HEADER_LOGGER_H

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>

class Logger {
public:
    Logger();

    // TODO: right now only accepts a single msg string, I should make
    // this arbitrary message elements using va_list, this makes the
    // interface match that of sprintf
    int info(const char *format, const char *msg);
    int warn(const char *format, const char *msg);
    int error(const char *format, const char *msg);
    
    int info(const char *msg);
    int warn(const char *msg);
    int error(const char *msg);

private:
    time_t t_;
    tm *time_point_;

    int msg_len_ { 0 };
    size_t time_len_ { 0 };

    static constexpr size_t time_buf_len_ { 30 };
    static constexpr size_t str_buf_len_ { 500 };
    static constexpr size_t max_str_ { 450 };

    char time_buf_[time_buf_len_];
    char str_buf_[str_buf_len_];

    static constexpr char err_str_[] = { "ERROR" };
    static constexpr char warn_str_[] = { "WARN" };
    static constexpr char info_str_[] = { "INFO" };
    
    int print_(FILE *stream, const char *log_type, 
            const char *format, const char *msg);
};

#endif
