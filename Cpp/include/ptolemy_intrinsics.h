#ifndef PTOLEMY_INTRINSICS_H
#define PTOLEMY_INTRINSICS_H

#include <cmath>
#include <cstdlib>
#include <algorithm>

// Fortran intrinsic function wrappers
// These map Fortran intrinsic names to C++ equivalents

// Double precision math
inline double DABS(double x)           { return std::fabs(x); }
inline double DSQRT(double x)          { return std::sqrt(x); }
inline double DLOG(double x)           { return std::log(x); }
inline double DLOG10(double x)         { return std::log10(x); }
inline double DEXP(double x)           { return std::exp(x); }
inline double DSIN(double x)           { return std::sin(x); }
inline double DCOS(double x)           { return std::cos(x); }
inline double DTAN(double x)           { return std::tan(x); }
inline double DASIN(double x)          { return std::asin(x); }
inline double DACOS(double x)          { return std::acos(x); }
inline double DATAN(double x)          { return std::atan(x); }
inline double DATAN2(double y, double x) { return std::atan2(y, x); }
inline double DMAX1(double a, double b)  { return std::max(a, b); }
inline double DMIN1(double a, double b)  { return std::min(a, b); }
inline double DMAX1(double a, double b, double c) { return std::max({a, b, c}); }
inline double DMIN1(double a, double b, double c) { return std::min({a, b, c}); }
inline double DSIGN(double a, double b)  { return b >= 0.0 ? std::fabs(a) : -std::fabs(a); }
inline double DMOD(double a, double b)   { return std::fmod(a, b); }
inline double DFLOAT(int i)            { return static_cast<double>(i); }
inline double DBLE(float x)            { return static_cast<double>(x); }
inline double DBLE(int x)              { return static_cast<double>(x); }
inline double DBLE(double x)           { return x; }
inline double DIMAG(std::complex<double> z) { return z.imag(); }
inline double DREAL(std::complex<double> z) { return z.real(); }
inline double DNINT(double x)          { return std::round(x); }
inline double ANINT(double x)          { return std::round(x); }

// Single precision math
inline float AMAX1(float a, float b)   { return std::max(a, b); }
inline float AMIN1(float a, float b)   { return std::min(a, b); }
inline float ALOG(float x)             { return std::log(x); }
inline float ALOG10(float x)           { return std::log10(x); }

// Generic math (Fortran generic names)
inline double ABS(double x)            { return std::fabs(x); }
inline int    ABS(int x)               { return std::abs(x); }
inline float  ABS(float x)             { return std::fabs(x); }
inline double SQRT(double x)           { return std::sqrt(x); }
inline double LOG(double x)            { return std::log(x); }
inline double EXP(double x)            { return std::exp(x); }
inline double SIN(double x)            { return std::sin(x); }
inline double COS(double x)            { return std::cos(x); }
inline double ATAN2(double y, double x) { return std::atan2(y, x); }

// Integer functions
inline int IABS(int x)                 { return std::abs(x); }
inline int ISIGN(int a, int b)         { return b >= 0 ? std::abs(a) : -std::abs(a); }
inline int MOD(int a, int b)           { return a % b; }
inline int MAX0(int a, int b)          { return std::max(a, b); }
inline int MIN0(int a, int b)          { return std::min(a, b); }
inline int MAX0(int a, int b, int c)   { return std::max({a, b, c}); }
inline int MIN0(int a, int b, int c)   { return std::min({a, b, c}); }
inline int MAX0(int a, int b, int c, int d) { return std::max({a, b, c, d}); }
inline int MIN0(int a, int b, int c, int d) { return std::min({a, b, c, d}); }
inline int MAX(int a, int b)           { return std::max(a, b); }
inline int MIN(int a, int b)           { return std::min(a, b); }
inline double MAX(double a, double b)  { return std::max(a, b); }
inline double MIN(double a, double b)  { return std::min(a, b); }
inline int IDINT(double x)             { return static_cast<int>(x); }
inline int INT(double x)               { return static_cast<int>(x); }
inline int INT(float x)                { return static_cast<int>(x); }
inline int NINT(double x)              { return static_cast<int>(std::round(x)); }
inline double FLOAT(int i)             { return static_cast<double>(i); }
inline double SNGL(double x)           { return x; }  // no-op in our context

// Type conversion
inline int IFIX(float x)               { return static_cast<int>(x); }
inline int IDNINT(double x)            { return static_cast<int>(std::round(x)); }

// Bit manipulation (Fortran intrinsics)
inline bool BTEST(int val, int bit)     { return (val >> bit) & 1; }
inline int IBSET(int val, int bit)      { return val | (1 << bit); }
inline int IBCLR(int val, int bit)      { return val & ~(1 << bit); }
inline int ISHFT(int val, int shift)    { return shift >= 0 ? (val << shift) : (val >> (-shift)); }
inline int IAND(int a, int b)           { return a & b; }
inline int IOR(int a, int b)            { return a | b; }
inline int IEOR(int a, int b)           { return a ^ b; }
inline int NOT(int a)                   { return ~a; }

// Character functions
inline int LEN(const char* s)           { return static_cast<int>(std::strlen(s)); }
inline int ICHAR(char c)               { return static_cast<int>(static_cast<unsigned char>(c)); }
inline char CHAR_F(int i)              { return static_cast<char>(i); }

// MAX/MIN with mixed int/double (Fortran allows this)
inline double MAX1(float a, float b)   { return std::max(a, b); }
inline double MIN1(float a, float b)   { return std::min(a, b); }

#endif // PTOLEMY_INTRINSICS_H
