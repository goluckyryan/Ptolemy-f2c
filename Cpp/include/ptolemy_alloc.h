#ifndef PTOLEMY_ALLOC_H
#define PTOLEMY_ALLOC_H

#include "ptolemy_commons.h"

//
// Type-punned accessors for the ALLOC memory pool.
// Matches Fortran EQUIVALENCE:
//   EQUIVALENCE ( facfr4, ALLOC(1), ILLOC(1), ILLOC2(1),
//                 ALLOC4(1), chlloc(1), dalloc(1) )
//
// Fortran ALLOC(1) = FACFR4 = byte offset 0 of /ALLOCS/
// C++ index: ALLOC(i) uses [i-1] for 1-based Fortran indexing
//

// Access as REAL*8 (double) — index i is 1-based
inline double& ALLOC(int i) {
    return reinterpret_cast<double*>(&ALLOCS.FACFR4)[i-1];
}

// Access as REAL*4 (float) — index i is 1-based
inline float& ALLOC4(int i) {
    return reinterpret_cast<float*>(&ALLOCS.FACFR4)[i-1];
}

// Access as INTEGER (int) — index i is 1-based
inline int& ILLOC(int i) {
    return reinterpret_cast<int*>(&ALLOCS.FACFR4)[i-1];
}

// Access as INTEGER*2 (short) — index i is 1-based
inline int16_t& ILLOC2(int i) {
    return reinterpret_cast<int16_t*>(&ALLOCS.FACFR4)[i-1];
}

// Access as CHARACTER*8 — index i is 1-based
inline char8& DALLOC(int i) {
    return reinterpret_cast<char8*>(&ALLOCS.FACFR4)[i-1];
}

// Access as CHARACTER*1 — index i is 1-based
inline char& CHLLOC(int i) {
    return reinterpret_cast<char*>(&ALLOCS.FACFR4)[i-1];
}

// Fortran array-passing helper: returns a pointer such that ptr[1] = ALLOC(i)
// Use when passing ALLOC(L) to a Fortran subroutine that uses 1-based array indexing
// Fortran: CALL FOO(ALLOC(L)) → FOO's VRAY(1) = ALLOC(L)
// C++:     FOO(ALLOC_base(L)) → FOO's VRAY[1] = ALLOC(L)
inline double* ALLOC_base(int i) {
    return &ALLOC(i) - 1;  // so that ptr[1] = ALLOC(i)
}

inline float* ALLOC4_base(int i) {
    return &ALLOC4(i) - 1;  // so that ptr[1] = ALLOC4(i)
}

inline int* ILLOC_base(int i) {
    return &ILLOC(i) - 1;  // so that ptr[1] = ILLOC(i)
}

#endif // PTOLEMY_ALLOC_H
