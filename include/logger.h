
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <memory>

class Logger {
public:
    Logger();
    ~Logger();

    int info(const char *format, const char *msg);
    int warn(const char *format, const char *msg);
    int error(const char *format, const char *msg);
    
    int info(const char *msg);
    int warn(const char *msg);
    int error(const char *msg);

private:
    time_t t_;
    tm *time_point_;
    char *time_buf_;
    char *str_buf_;

    int msg_len_ { 0 };
    size_t time_len_ { 0 };

    static const size_t time_buf_len_ { 30 };
    static const size_t str_buf_len_ { 500 };
    static const size_t max_str_ { 450 };

    static constexpr char err_str_[] = { "ERROR" };
    static constexpr char warn_str_[] = { "WARN" };
    static constexpr char info_str_[] = { "INFO" };
    
    int print_(FILE *stream, const char *log_type, 
            const char *format, const char *msg);
    void empty_time_buf_();
};
