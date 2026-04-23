#ifndef PTOLEMY_TYPES_H
#define PTOLEMY_TYPES_H

#include <cstdint>
#include <complex>
#include <cstring>

// Fortran type aliases for verbatim transliteration
using real8 = double;
using real4 = float;
using integer = int;
using integer2 = int16_t;
using integer4 = int32_t;
using integer8 = int64_t;
using logical = int;       // Fortran LOGICAL is typically 4 bytes
using logical1 = int8_t;   // LOGICAL*1
using complex16 = std::complex<double>;

// Fortran logical constants
constexpr int TRUE_F  = 1;
constexpr int FALSE_F = 0;

// Helper to convert Fortran logical to C++ bool
inline bool ftobool(int val) { return val != 0; }
inline bool ftobool(int8_t val) { return val != 0; }

// Fortran CHARACTER*8 helper
struct char8 {
    char data[8];
    char8() { std::memset(data, ' ', 8); }
    char8(const char* s) {
        std::memset(data, ' ', 8);
        if (s) {
            size_t len = std::strlen(s);
            if (len > 8) len = 8;
            std::memcpy(data, s, len);
        }
    }
    bool operator<=(const char8& other) const { return std::memcmp(data, other.data, 8) <= 0; }
    bool operator>=(const char8& other) const { return std::memcmp(data, other.data, 8) >= 0; }
    bool operator==(const char8& other) const { return std::memcmp(data, other.data, 8) == 0; }
    bool operator!=(const char8& other) const { return std::memcmp(data, other.data, 8) != 0; }
    bool operator<(const char8& other) const { return std::memcmp(data, other.data, 8) < 0; }
    bool operator>(const char8& other) const { return std::memcmp(data, other.data, 8) > 0; }
};

#endif // PTOLEMY_TYPES_H
