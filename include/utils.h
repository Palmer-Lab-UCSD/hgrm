#ifndef HEADER_UTILS_H
#define HEADER_UTILS_H

#include <cstdint>

namespace utils {

// Claude recommends that I define serialization explicitly so
// that the bit packed field order is not ABI dependent
struct Version {
    uint8_t major;
    uint16_t minor;
    uint16_t micro;

    // pack as 8, 12, 12 for 32 bits in total
    uint32_t pack() const {
        return (static_cast<uint32_t>(major) << 24
            | static_cast<uint32_t>(minor & 0x0FFF) << 12 
            | static_cast<uint32_t>(micro & 0x0FFF));
    }

    static Version unpack(uint32_t vnum) {
        return Version {
            static_cast<uint8_t>(vnum >> 24),
            static_cast<uint16_t>(vnum >> 12 & 0x0FFF),
            static_cast<uint16_t>(vnum & 0x0FFF)
        };
    }
};

// template<typename T>
// struct Array {
//     Array(size_t size_in): size(size_in), 
//         data(size > 0 ? new T[size] : nullptr) {};
// 
//     ~Array() { if (data) delete[] data; };
// 
//     size_t size;
//     T *data;
//     size_t len = 0;
// 
//     //unsafe referencing
//     T operator[](size_t i) { return data[i]; };
//     T& operator[](size_t i) { return data[i]; };
// 
//     STATUS append(T val) {
//         if (len >= size-1)
//             return END_OF_BUF_ERROR;
// 
//         data[len++] = val;
//         return SUCCESS;
//     }
// 
//     void fill(T val) {
//         std::memset(data, val, size);
//         len = 0;
//     }
// };

}

#endif
