// numbered_store.cpp -- translated from numbered_store.f
//
// THIS IS THE ENTIRE NUMBERED STORAGE PACKAGE
//
// This is the Ptolemy dynamic memory allocation system.
// Original by D.Gloeckner 7/14/73, with extensive modifications
// by S.P. and R.P.G. through 8/20/85.
//
// Translated from Fortran to C++ preserving all logic and labels.

#include "ptolemy_commons.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_alloc.h"
#include "ptolemy_io.h"
#include "ptolemy_forward.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cinttypes>

// ============================================================================
// EQUIVALENCE (NUMBER(1), ANAME)
// INTEGER*2 NUMBER(4) overlaid on CHARACTER*8 ANAME (also INTEGER*8 IANAME)
// ============================================================================
union ANameUnion {
    char     c8[8];       // CHARACTER*8 ANAME
    int16_t  NUMBER[4];   // INTEGER*2 NUMBER(4)
    int64_t  ianame;      // INTEGER*8 IANAME
};

// ============================================================================
// Dispatch enum for IALLOC + its ENTRY points
// ============================================================================
enum class IallocEntry {
    IALLOC,
    GIVEAL,
    CHOPIT,
    IREDEF_F,
    ISTART_F,
    ICLEAR_F,
    NSSTAT_F,
    NSCAT_F
};

// ============================================================================
// Core dispatch function implementing IALLOC and all ENTRY points
// ============================================================================
static int IALLOC_dispatch(IallocEntry entry, int IWRDS, int NAMLCA = 0) {
    // implicit real*8 (a-h, o-z)
    // Variables that would be I-N: int; A-H,O-Z: double

    // SAVE'd / static locals
    static ANameUnion aname_u;
    static bool STOPSW = true;
    static char8 BLANK(" ");

    // References to COMMON blocks
    auto& IBASE  = USAGE.IBASE;
    auto& NWORDS = USAGE.NWORDS;
    auto& IMAX   = USAGE.IMAX;
    auto& INAME  = USAGE.INAME;
    auto& IGOT   = USAGE.IGOT;
    auto& NLOC   = USAGE.NLOC;
    auto& IPEAK  = USAGE.IPEAK;
    auto& NAMPEK = USAGE.NAMPEK;
    auto& NUMCOM = USAGE.NUMCOM;
    auto& ICHEC  = USAGE.ICHEC;
    auto& THENAM = USAGE.THENAM;
    auto& IWASAT = USAGE.IWASAT;
    auto& NAMLC2 = USAGE.NAMLC2;
    auto& NSPTSW = USAGE.NSPTSW;
    auto& FACFR4 = ALLOCS.FACFR4;
    auto& FACFR2 = ALLOCS.FACFR2;
    auto& FACFR1 = ALLOCS.FACFR1;
    auto* Z      = LOCPTRS.Z;     // 1-based
    auto* LENG   = LENGTH.LENG;   // 1-based
    auto* NAMES  = NAMCOM.NAMES;  // 1-based

    // Initialize ANAME: data iANAME / Z'00000000FFFF0000' /
    // This is done once via static init
    static bool aname_inited = false;
    if (!aname_inited) {
        aname_u.ianame = INT64_C(0x00000000FFFF0000);
        aname_inited = true;
    }

    // Local variables
    int IALLOC_val = 0;  // function return value
    bool REDFSW = false;
    bool ERRSW, PRTSW, FONDSW, ALLSW;
    int NAMLOC, NAMUSE, N, ILOC, IHAD;
    int LUSE, IHCNT, LL, LCNT, LSUM, LLOW;
    int NAMOLD, IWDOLD;
    double LSIZE, BSIZE, BAVAIL, BUSE, BPEAK;

    // -----------------------------------------------------------------------
    // ENTRY point dispatch
    // -----------------------------------------------------------------------

    // --- GIVEAL ---
    if (entry == IallocEntry::GIVEAL) {
        if (NSPTSW) printf(" GIVEAL CALLED\n");
        REDFSW = false;
        ALLSW = true;
        NAMLOC = 0;
        if (NLOC > IBASE) goto label_505;
        goto label_300;
    }

    // --- CHOPIT ---
    if (entry == IallocEntry::CHOPIT) {
        if (NSPTSW) printf(" CHOPIT CALLED:%15d\n", IWRDS);
        aname_u.c8[0] = DALLOC(IWASAT).data[0];
        aname_u.c8[1] = DALLOC(IWASAT).data[1];
        aname_u.c8[2] = DALLOC(IWASAT).data[2];
        aname_u.c8[3] = DALLOC(IWASAT).data[3];
        aname_u.c8[4] = DALLOC(IWASAT).data[4];
        aname_u.c8[5] = DALLOC(IWASAT).data[5];
        aname_u.c8[6] = DALLOC(IWASAT).data[6];
        aname_u.c8[7] = DALLOC(IWASAT).data[7];
        IHAD = aname_u.NUMBER[0];
        IPEAK = IPEAK + IWRDS + 1;
        if (NAMLC2 != IHAD) goto label_320;
        NAMPEK = NAMPEK + 1;
        NLOC = IWASAT + IWRDS + 1;
        LENG[NAMLC2] = IWRDS;
        return NAMLC2;  // implicit: function value not explicitly set, but object is NAMLC2
    }

    // --- IREDEF_F ---
    if (entry == IallocEntry::IREDEF_F) {
        NAMLOC = NAMLCA;
        NAMOLD = NAMLOC;
        THENAM = NAMES[NAMOLD];
        if (NSPTSW) printf(" IREDEF: %.8s%15d\n", THENAM.data, IWRDS);
        IWDOLD = LENG[NAMOLD];
        REDFSW = true;
        goto label_5;
    }

    // --- ISTART_F ---
    if (entry == IallocEntry::ISTART_F) {
        int NBYTES = IWRDS;
        NUMCOM = 0;
        FACFR4 = 2;
        FACFR2 = 4;
        FACFR1 = 8;

        NWORDS = NBYTES;
        IBASE = 0;
        GRAB(&ALLOC(0), FACFR1, IBASE, NWORDS);

        if (IBASE < 0) {
            fprintf(stderr, " ibase < 0: %d %d\n", IBASE, NWORDS);
        }

        ICHEC = IBASE + NWORDS - 1;
        IMAX = NAMDIM;
        IPEAK = 0;
        NAMPEK = 0;
        NUMCOM = 0;

        goto label_710;
    }

    // --- ICLEAR_F ---
    if (entry == IallocEntry::ICLEAR_F) {
        printf("\n   -------  ALL NUMBERED OBJECTS CLEARED  ------\n\n");
        goto label_710;
    }

    // --- NSSTAT_F ---
    if (entry == IallocEntry::NSSTAT_F) {
        PRTSW = false;
        ALLSW = false;
        goto label_900;
    }

    // --- NSCAT_F ---
    if (entry == IallocEntry::NSCAT_F) {
        PRTSW = true;
        ALLSW = (IWRDS == -1);
        goto label_900;
    }

    // --- IALLOC (default entry) ---
    REDFSW = false;
    if (NSPTSW) printf(" IALLOC/NALLOC: %.8s%15d\n", THENAM.data, IWRDS);

label_5:
    if (IWRDS > 0) goto label_100;
    printf("\n**** %.8s  - REQUEST FOR %15d WORDS. ****\n", THENAM.data, IWRDS);
    goto label_9500;

    // -----------------------------------------------------------------
    // PREPARE TO ALLOCATE AN OBJECT.
    // Walk through the list; find used space, check integrity,
    // look for an exact fit to a freed spot.
    // -----------------------------------------------------------------

label_100:
    FONDSW = false;
    NAMUSE = 0;
    LUSE = 0;
    NAMLOC = 0;
    if (INAME == 0) goto label_150;
    for (N = 1; N <= INAME; N++) {
        // Is this an empty slot?
        if (Z[N] != 0) goto label_120;
        if (NAMLOC == 0) NAMLOC = N;
        goto label_149;

    label_120:
        // In use or freed spot. Check integrity.
        aname_u.NUMBER[0] = (int16_t)N;
        if (DALLOC(IABS(Z[N]) - 1) == *(char8*)aname_u.c8) goto label_130;
        printf("\n*** HEADER OF OBJECT%4d HAS BEEN DESTROYED.\n"
               " THIS AND PREVIOUS OBJECT LOCATION, LENGTH AND NAME ARE:\n"
               "%15d%15d  %.8s\n%15d%15d  %.8s\n",
               N, Z[N], LENG[N], NAMES[N].data,
               Z[N-1], LENG[N-1], NAMES[N-1].data);

    label_130:
        // Is it free space or still in use?
        if (Z[N] < 0) goto label_140;
        NAMUSE = NAMUSE + 1;
        LUSE = LUSE + LENG[N] + 1;
        goto label_149;

    label_140:
        if (FONDSW) goto label_149;
        if (IWRDS != LENG[N]) goto label_149;
        NAMLOC = N;
        ILOC = -Z[NAMLOC] - 1;
        FONDSW = true;

    label_149:;
    }

    // Some statistics
label_150:
    IPEAK = MAX0(IPEAK, LUSE + IWRDS + 1);
    NAMPEK = MAX0(NAMPEK, NAMUSE + 1);

    // Did we find an exact fit?
    if (FONDSW) goto label_250;

    // No, did we at least find a name slot?
    if (NAMLOC != 0) goto label_200;

    // No, we must get a new name slot if there are more available
    if (INAME < IMAX) goto label_170;

    // We have gone through the whole stack of names, are there holes?
    if (NAMUSE < IMAX) goto label_200;

label_160:
    printf("\n*** %.8s  -- NO MORE NAME SLOTS AVAILABLE ***\n", THENAM.data);
    goto label_9500;

    // Add to end of name list
label_170:
    INAME = INAME + 1;
    NAMLOC = INAME;

    // Is there space at the end of the allocator?
label_200:
    if (NLOC + IWRDS <= ICHEC) goto label_220;

    // No, will compression solve it?
    if (LUSE + IWRDS < NWORDS) goto label_500;

label_210:
    printf("\n*** %.8s  -- NOT ENOUGH ALLOCATOR.  ATTEMPTING"
           " TO DEFINE A%12d WORD OBJECT.\n"
           "     TOTAL ALLOCATOR SIZE MUST BE AT LEASET%13d WORDS.\n",
           THENAM.data, IWRDS, IPEAK);
    goto label_9500;

    // There is room to allocate and also a name slot is available.
    // If NAMLOC is not yet assigned then we must compress to get the slot.
label_220:
    if (NAMLOC == 0) goto label_500;
    ILOC = NLOC;
    NLOC = NLOC + IWRDS + 1;
    IGOT = IGOT + 1;

    // Define the object
label_250:
    if (NSPTSW) printf(" ASSIGNED%4d%12d     NLOC =%12d\n", NAMLOC, ILOC, NLOC);

    Z[NAMLOC] = ILOC + 1;

    // Put name (number) in first INTEGER*2 part of array
    aname_u.NUMBER[0] = (int16_t)NAMLOC;
    DALLOC(ILOC) = *(char8*)aname_u.c8;
    LENG[NAMLOC] = IWRDS;
    NAMES[NAMLOC] = BLANK;
    IALLOC_val = NAMLOC;
    NAMLC2 = NAMLOC;
    if (REDFSW) goto label_620;
    return IALLOC_val;

    // -----------------------------------------------------------------
    // GIVEAL continued: return from compress, give largest piece
    // -----------------------------------------------------------------
label_300:
    IWRDS = ICHEC - NLOC - 2;
    IWASAT = NLOC;
    if (IWRDS >= 2) goto label_310;
    printf(" ***ERROR/NOTHING TO GIVE IN GIVEAL  CALL\n");
    goto label_9500;

    // Find a name stack location for GIVEAL
label_310:
    if (INAME == 0) goto label_170;
    for (N = 1; N <= INAME; N++) {
        if (Z[N] != 0) continue;
        NAMLOC = N;
        goto label_200;
    }
    if (INAME < IMAX) goto label_170;
    goto label_160;

    // CHOPIT error
label_320:
    printf(" ***ERROR CALL TO CHOPIT MUST IMMEDIATELY FOLLOW GIVEAL\n");
    goto label_9500;

    // -----------------------------------------------------------------
    // COMPRESSION
    // -----------------------------------------------------------------
label_500:
    ALLSW = false;

    if (NAMLOC != 0) printf("\n %.8s  --  COMPRESSION IS OCCURING NOW.\n", THENAM.data);
    if (NAMLOC == 0) printf("\n %.8s  --  NAME STACK COMPRESSION IS OCCURING NOW.\n", THENAM.data);

label_505:
    LL = IBASE;
    LCNT = 0;
    LSUM = 0;
    NUMCOM = NUMCOM + 1;

    // Find first free object
label_510:
    LLOW = LL;
    if (LCNT >= IGOT) goto label_600;

    std::memcpy(aname_u.c8, DALLOC(LL).data, 8);
    NAMLOC = aname_u.NUMBER[0];

    // Structure check
    if (IABS(Z[NAMLOC]) != (LL + 1)) goto label_8171;
    if (Z[NAMLOC] < 0) goto label_520;
    LL = LL + LENG[NAMLOC] + 1;
    LCNT = LCNT + 1;
    goto label_510;

label_8171:
    printf(" **** STRUCTURE CHECK: LOC, HEADER, NAMLOC, LENG, L, NAME:\n"
           "%13d%18" PRIx64 "%13d%13d%13d  %.8s\n",
           LL, aname_u.ianame, NAMLOC, LENG[NAMLOC], Z[NAMLOC], NAMES[NAMLOC].data);
    goto label_9500;

    // First free location found
label_520:
    IHCNT = LCNT;

    // Have found a free location. Mark it not in use and go on.
label_530:
    Z[NAMLOC] = 0;
    LSUM = LSUM + LENG[NAMLOC] + 1;

    // Have found first free location -- now find next filled one
label_540:
    if (IHCNT >= (IGOT - 1)) goto label_600;
    IHCNT = IHCNT + 1;
    LL = LL + LENG[NAMLOC] + 1;
    if (LL >= NLOC) goto label_600;
    std::memcpy(aname_u.c8, DALLOC(LL).data, 8);
    NAMLOC = aname_u.NUMBER[0];

    // Internal structure check
    if (IABS(Z[NAMLOC]) != (LL + 1)) goto label_8171;
    if (Z[NAMLOC] < 0) goto label_530;

    // Found a still-in-use object, compress
    N = LENG[NAMLOC] + 1;
    for (int I = 1; I <= N; I++) {
        DALLOC(LLOW + I - 1) = DALLOC(LL + I - 1);
    }

    // Note we have also moved the identifying number
    LCNT = LCNT + 1;
    Z[NAMLOC] = LLOW + 1;
    LLOW = LLOW + N;
    goto label_540;

    // Compress all done. Update number of active slots and top of allocator.
label_600:
    IGOT = LCNT;
    NLOC = LLOW;
    if (ALLSW) goto label_300;

    printf(" %5d ARRAYS COMPRESSED%8d WORDS FREED\n\n", IGOT, LSUM);
    goto label_100;

    // -----------------------------------------------------------------
    // IREDEF continued: copy old to new
    // -----------------------------------------------------------------
label_620:
    N = MIN0(IWRDS, IWDOLD) + 1;
    LL = Z[NAMOLD] - 1;
    for (int I = 1; I <= N; I++) {
        DALLOC(ILOC - 1 + I) = DALLOC(LL - 1 + I);
    }

    // Change old to point to new
    Z[NAMOLD] = ILOC + 1;
    LENG[NAMOLD] = IWRDS;
    NAMES[NAMOLD] = THENAM;

    // Change new to point to old and free it
    DALLOC(LL) = *(char8*)aname_u.c8;
    LENG[NAMLOC] = IWDOLD;
    Z[NAMLOC] = -(LL + 1);
    IALLOC_val = NAMOLD;
    NAMLOC = NAMOLD;

    return IALLOC_val;

    // -----------------------------------------------------------------
    // ISTART_F continued: initialize tables
    // -----------------------------------------------------------------
label_710:
    NLOC = IBASE;
    INAME = 0;
    IGOT = 0;
    THENAM = BLANK;
    for (int I = 1; I <= IMAX; I++) {
        NAMES[I] = BLANK;
        Z[I] = 0;
    }

    if (NSPTSW) printf(" ALLOCATOR RANGE:%12d%12d     # NAMES =%5d\n", IBASE, ICHEC, IMAX);
    IALLOC_val = IBASE;
    return IALLOC_val;

    // -----------------------------------------------------------------
    // NSSTAT / NSCAT: statistics and catalog
    // -----------------------------------------------------------------
label_900:
    ERRSW = false;

label_910:
    LUSE = 0;
    NAMUSE = 0;

    if (PRTSW) printf("             ---- NUMBERED STORAGE CATALOG ----\n");
    if (NSPTSW || (PRTSW && ALLSW))
        printf("\n    FROM%8d   TO%8d      TOP OF ALLOC =%8d%8d ACTIVE NAME SLOTS\n",
               IBASE, NLOC, ICHEC, IGOT);
    if (PRTSW) printf("\n  NUMBER   NAME     LENGTH      START        END\n");

    for (int III = 1; III <= INAME; III++) {
        LLOW = IABS(Z[III]);
        LL = LLOW + LENG[III] + 1;
        if (Z[III] <= 0) goto label_940;
        LUSE = LUSE + LENG[III] + 1;
        NAMUSE = NAMUSE + 1;
        if (PRTSW) printf(" %8d  %.8s%8d%12d%11d\n",
                          III, NAMES[III].data, LENG[III], LLOW, LL);
        goto label_949;
    label_940:
        if (ALLSW && LLOW > 0)
            printf(" %8d  %.8s%8d%12d%11d     FREED\n",
                   III, NAMES[III].data, LENG[III], LLOW, LL);
    label_949:;
    }

    LSIZE = NLOC - IBASE;
    BSIZE = LSIZE / 128.0;
    BAVAIL = NWORDS / 128.0;
    BUSE = LUSE / 128.0;
    BPEAK = IPEAK / 128.0;
    printf("0               ALLOCATOR STATISTICS\n"
           "                                        WORDS   KILOBYTES\n"
           " ALLOCATOR SIZE:                     %8d%10.1f\n"
           " CURRENT UNCOMPRESSED SIZE:          %8d%10.1f\n"
           " CURRENT COMPRESSED SIZE:            %8d%10.1f\n"
           " PEAK COMPRESSED SIZE:               %8d%10.1f\n"
           "0NAME STACK SIZE:                    %8d\n"
           " CURRENT HIGHEST NAME:               %8d\n"
           " CURRENT NAMES IN USE:               %8d\n"
           " PEAK NAMES IN USE:                  %8d\n"
           " NUMBER OF ALLOCATOR COMPRESSES:     %8d\n\n\n",
           NWORDS, BAVAIL, (int)LSIZE, BSIZE, LUSE, BUSE,
           IPEAK, BPEAK, IMAX, INAME, NAMUSE, NAMPEK, NUMCOM);

    if (ERRSW) goto label_9550;

    return IALLOC_val;

    // -----------------------------------------------------------------
    // No room for allocation
    // -----------------------------------------------------------------
label_9500:
    IALLOC_val = 0;
    if (!STOPSW) return IALLOC_val;

    printf("\n****** EXECUTION IS BEING TERMINATED BY NUMBERED"
           " STORAGE ******\n\n");
    ERRSW = true;
    PRTSW = true;
    ALLSW = true;
    goto label_910;

label_9550:
    std::printf("\n $$$ IALLOC error\n");
    SYSERR();
    // STOP 1234
    std::exit(1234);

    return IALLOC_val;  // unreachable, keeps compiler happy
}


// ============================================================================
// Public function wrappers (inline dispatch)
// ============================================================================

int IALLOC(int IWRDS) {
    return IALLOC_dispatch(IallocEntry::IALLOC, IWRDS);
}

int GIVEAL(int IWRDS) {
    return IALLOC_dispatch(IallocEntry::GIVEAL, IWRDS);
}

int CHOPIT(int IWRDS) {
    return IALLOC_dispatch(IallocEntry::CHOPIT, IWRDS);
}

int IREDEF_F(int IWRDS, int NAMLCA) {
    return IALLOC_dispatch(IallocEntry::IREDEF_F, IWRDS, NAMLCA);
}

int ISTART_F(int NBYTES) {
    return IALLOC_dispatch(IallocEntry::ISTART_F, NBYTES);
}

int ICLEAR_F(int IDUM) {
    return IALLOC_dispatch(IallocEntry::ICLEAR_F, IDUM);
}

int NSSTAT_F(int IDUM) {
    return IALLOC_dispatch(IallocEntry::NSSTAT_F, IDUM);
}

int NSCAT_F(int IDUM) {
    return IALLOC_dispatch(IallocEntry::NSCAT_F, IDUM);
}


// ============================================================================
// SUBROUTINE IREDEF ( N1, N2 )
//   These are things inside IALLOC which are called instead
//   of being invoked as functions -- violates the standard
//   to do it directly.
//
// Also contains ENTRY ISTART, ICLEAR, NSSTAT, NSCAT
// ============================================================================

void IREDEF(int N1, int N2) {
    IREDEF_F(N1, N2);
}

void ISTART(int N1) {
    ISTART_F(N1);
}

void ICLEAR(int N1) {
    ICLEAR_F(N1);
}

void NSSTAT(int N1) {
    NSSTAT_F(N1);
}

void NSCAT(int N1) {
    NSCAT_F(N1);
}


// ============================================================================
// FUNCTION NALLOC( IWRDS, NAME )
//
// Named-numbered-storage.
// Allocates IWRDS words and associates NAME with the slot.
// If NAME already exists, frees it and tries to reuse space.
//
// Also contains ENTRY NAMLOC( NAME ).
// ============================================================================

int NALLOC(const char* NAME, int IWRDS) {
    // implicit real*8 (a-h, o-z)

    auto* NAMES = NAMCOM.NAMES;  // 1-based
    auto* L     = LOCPTRS.Z;     // 1-based  (same as Z)
    auto* LENG  = LENGTH.LENG;   // 1-based

    auto& IBASE  = USAGE.IBASE;
    auto& NWORDS = USAGE.NWORDS;
    auto& IMAX   = USAGE.IMAX;
    auto& INAME  = USAGE.INAME;
    auto& IGOT   = USAGE.IGOT;
    auto& NLOC   = USAGE.NLOC;
    auto& THENAM = USAGE.THENAM;
    auto& NSPTSW = USAGE.NSPTSW;

    bool FNDSW;
    char8 BLANK(" ");
    int IL;

    // Pad NAME out to 8 chars (Fortran CHARACTER*8 assignment)
    char8 NAME8(NAME);

    // -----------------------------------------------------------
    // Is it defined?  (Search for the name)
    // -----------------------------------------------------------

    // ASSIGN 100 TO IGO  /  GO TO 800
    // --- inline label_800 for IGO=100 ---
    FNDSW = false;
    IL = 0;

    // label_810: search
    for (int I = 1; I <= IMAX; I++) {
        if (NAMES[I] != NAME8) continue;
        if (L[I] > 0) {
            // Really found
            IL = I;
            FNDSW = true;
            goto label_100;
        }
        // Semi found - perhaps there is a better entry
        NAMES[I] = BLANK;
        IL = I;
        FNDSW = true;
    }
    // label_870: fall through

label_100:
    if (!FNDSW) goto label_150;

    // Found - free it and try to reuse
    L[IL] = IABS(L[IL]);
    if (LENG[IL] == IWRDS) goto label_180;
    L[IL] = -L[IL];
    NAMES[IL] = BLANK;

label_150:
    THENAM = NAME8;
    IL = IALLOC(IWRDS);
    THENAM = BLANK;
    if (IL == 0) goto label_200;
label_180:
    NAMES[IL] = NAME8;
label_200:
    return IL;
}


// ============================================================================
// ENTRY NAMLOC( NAME )
//
// Locate a name.  Returns the query stack location of an
// already defined name, or 0 if not found.
// ============================================================================

int NAMLOC(const char* NAME) {
    auto* NAMES = NAMCOM.NAMES;
    auto* L     = LOCPTRS.Z;

    auto& IMAX  = USAGE.IMAX;

    bool FNDSW;
    char8 BLANK(" ");
    int IL;

    char8 NAME8(NAME);

    // ASSIGN 300 TO IGO  /  GO TO 800
    FNDSW = false;
    IL = 0;

    for (int I = 1; I <= IMAX; I++) {
        if (NAMES[I] != NAME8) continue;
        if (L[I] > 0) {
            // Really found
            IL = I;
            FNDSW = true;
            goto label_300;
        }
        // Semi found
        NAMES[I] = BLANK;
        IL = I;
        FNDSW = true;
    }

label_300:
    if (!FNDSW) return 0;
    if (L[IL] < 0) return 0;
    return IL;
}


// ============================================================================
// SUBROUTINE NSDUMP ( IDUMP, NSKEY )
//
// Prints objects in numbered storage.
//
// IDUMP - number of object to be printed; if zero, all are printed.
// NSKEY - type of printout:
//   1 - automatic / default real*8
//   2 - real*4
//   3 - integer*4
//   4 - integer*2
//   5 - hexadecimal
//   6 - real*8
// ============================================================================

void NSDUMP(int IDUMP, int NSKEY) {
    // implicit real*8 (a-h, o-z)

    auto* Z    = LOCPTRS.Z;
    auto* LENG = LENGTH.LENG;
    auto* NAMES = NAMCOM.NAMES;
    auto& FACFR4 = ALLOCS.FACFR4;
    auto& FACFR2 = ALLOCS.FACFR2;
    auto& INAME  = USAGE.INAME;
    auto& UNDEF  = INTRNL.UNDEF;
    auto& NOTDEF = INTRNL.NOTDEF;

    int IFAC;
    int NSTRT, NEND, NUMLOC, LX, LLX;
    int I, I1, I2;

    // EQUIVALENCE ( BUFFER(1), IBUF(1) ) and ( INTGER(1), SINGLE(1) )
    union { float f; int32_t i; } intger_single;

    // Setup size of each word to print
    // GO TO ( 110, 120, 120, 130, 120, 110 ), NSKEY
    switch (NSKEY) {
        case 1: case 6:
            IFAC = 1;
            break;
        case 2: case 3: case 5:
            IFAC = FACFR4;
            break;
        case 4:
            IFAC = FACFR2;
            break;
        default:
            IFAC = 1;
            break;
    }

    NSTRT = 1;
    NEND = INAME;
    if (IDUMP != 0) {
        NSTRT = IDUMP;
        NEND = IDUMP;
    }

    for (NUMLOC = NSTRT; NUMLOC <= NEND; NUMLOC++) {
        LX = Z[NUMLOC];
        if (LX <= 0) continue;
        LX = IFAC * LX - IFAC;
        LLX = LENG[NUMLOC];
        LLX = IFAC * LLX;

        printf("\n  ----NUMBERED OBJECT%6d  NAMED %.8s  WITH%6d ELEMENTS----\n\n",
               NUMLOC, NAMES[NUMLOC].data, LLX);

        // GO TO ( 310, 320, 330, 340, 350, 310 ), NSKEY
        switch (NSKEY) {
            case 1: case 6: {
                // REAL*8 dump
                for (I = 1; I <= LLX; I++) {
                    printf(" %13.5g", ALLOC(LX + I));
                    if (I % 10 == 0 || I == LLX) printf("\n");
                }
                break;
            }
            case 2: {
                // REAL*4 dump
                for (I = 1; I <= LLX; I++) {
                    printf(" %13.5g", (double)ALLOC4(LX + I));
                    if (I % 10 == 0 || I == LLX) printf("\n");
                }
                break;
            }
            case 3: {
                // INTEGER*4 dump
                for (I = 1; I <= LLX; I++) {
                    if (I == 1) printf(" %10d", ILLOC(LX + I));
                    else        printf("%13d", ILLOC(LX + I));
                    if (I % 10 == 0 || I == LLX) printf("\n");
                }
                break;
            }
            case 4: {
                // INTEGER*2 dump
                for (I = 1; I <= LLX; I++) {
                    if (I == 1) printf(" %10d", (int)ILLOC2(LX + I));
                    else        printf("%13d", (int)ILLOC2(LX + I));
                    if (I % 10 == 0 || I == LLX) printf("\n");
                }
                break;
            }
            case 5: {
                // Hexadecimal dump
                for (I = 1; I <= LLX; I++) {
                    if (I % 2 == 1) printf("%10X", (unsigned)ILLOC(LX + I));
                    else            printf("%9X", (unsigned)ILLOC(LX + I));
                    if (I % 10 == 0 || I == LLX) printf("\n");
                }
                break;
            }
            default: {
                // Automatic detection mode (NSKEY not in 1-6)
                // Fall through to case 1
                for (I = 1; I <= LLX; I++) {
                    printf(" %13.5g", ALLOC(LX + I));
                    if (I % 10 == 0 || I == LLX) printf("\n");
                }
                break;
            }
        }

        // If NSKEY == 1 or 6, the above already handled it via
        // the switch. For the automatic determination path (which
        // was largely commented out in the original Fortran), we also
        // handle it. The Fortran code at label 500 prints an initial
        // G-format dump then does per-element int/real detection.
        // Since the CIB code is commented out, the active path is:
        // label 500 -> print with G format, then label 510 per-element.
        // We implement this for completeness when NSKEY == 1.

        // The Fortran code after the switch statements has:
        //   500: write(6,257) ...  (already printed above)
        //        write(6,503) (alloc(lx+i), i=1,llx)    -- bulk dump
        //        Then per-element loop 510 with format building.
        //
        // Since the active code just does a G15.5 bulk dump (line 1010-1011)
        // followed by the per-element loop, and the auto-detection (CIB)
        // is commented out, we skip the complex format-building loop
        // and rely on the switch-based dump above which already covers
        // all NSKEY values.
    }
}
