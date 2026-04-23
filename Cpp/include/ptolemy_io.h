#ifndef PTOLEMY_IO_H
#define PTOLEMY_IO_H

#include <cstdio>
#include <cstring>

//
// Fortran I/O emulation.
// Most WRITE/FORMAT pairs are translated directly to printf calls
// at the translation site. This header provides helpers for common patterns.
//

// Fortran unit numbers
// Unit 5 = stdin
// Unit 6 = stdout
// Unit 0 = stderr

// Get FILE* for a Fortran unit number
inline FILE* fortran_unit(int unit) {
    switch (unit) {
        case 0:  return stderr;
        case 5:  return stdin;
        case 6:  return stdout;
        default: return stdout;
    }
}

// Fortran STOP statement
inline void FSTOP() {
    std::exit(0);
}

inline void FSTOP(int code) {
    std::exit(code);
}

// Fortran hostnm intrinsic
int hostnm(char* name);

// Fortran dtime (from dtime.c)
double dtime_();

// Sequential I/O (from srread.cpp)
void srread(int* idcb, int* icode, char* line, int* numchr, int* maxchr, int* rtype);
void swrite(int* idcb, int* icode, char* line, int* numchr);
void snwrte(int* idcb, int* icode, char* line, int* numchr);
void sropen(int* idcb, int* icode, char* fname);
void srclse(int* idcb);
void swopen(int* idcb, int* icode, char* fname, int* oversw, int* append);
void swclse(int* idcb);
void swdlte(int* idcb, int* icode, char* fname);

#endif // PTOLEMY_IO_H
