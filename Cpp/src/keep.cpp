// keep.cpp -- translated from keep.f
//
// KEEPIT - performs Speakeasy KEEP operation for Ptolemy
//
// Ptolemy input is:  KEEP  OBJECTNAME  MEMBERNAME
//
// OBJECTNAME may be an allocator object or a special name:
//   ANGLEGRID, HEADER, REACTION, ELAB, SMAG, SPHASE, etc.
//
// 11/7/76 - first version
// 12/14/18 - crosssections, F, etc by angle grid - S.P.
// 3/19/04 - Speakeasy IV keep call - s.p.
//
// Translated from Fortran to C++ preserving all logic and labels.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_forward.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_alloc.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>

void KEEPIT()
{
    // implicit real*8 (a-h, o-z)

    // References to COMMON blocks
    auto* Z    = LOCPTRS.Z;     // /LOCptrs/ Z(namdim), 1-based
    auto* LENG = LENGTH.LENG;   // /LENGTH/  LENG(1), 1-based

    // /FLOAT/ variables used
    auto& ANGSTP = FLOAT_common.ANGSTP;
    auto& ANGMIN = FLOAT_common.ANGMIN;
    auto& ANGMAX = FLOAT_common.ANGMAX;
    auto& ELAB   = FLOAT_common.ELAB;

    // /HEDCOM/
    auto* REACT  = HEDCOM.REACT;    // 1-based: REACT[1]..REACT[45]
    auto* HEADER = HEDCOM.HEADER;   // 1-based: HEADER[1]..HEADER[65]

    // /INELCM/
    auto& NMLOLX = INELCM.NMLOLX;

    // Local variables
    char OBNAME[8];
    char MEMBER[8];
    char BLANK = ' ';
    char8 AS("AS");

    static const int NSPECL = 16;

    // Special object names (CHARACTER*8, space-padded)
    static const char8 SPECIL[16] = {
        char8("SMAG"),     char8("SPHASE"),   char8("LXCROSSS"),
        char8("MXCROSSS"), char8("SIN"),      char8("SOUT"),
        char8("F"),        char8("HEADER"),   char8("REACTION"),
        char8("ELAB"),     char8("ANGLEGRI"), char8("CHIGNVEC"),
        char8("S"),        char8("ISUM"),     char8("SEL"),
        char8("CROSSSEC")
    };
    // 0-based arrays matching SPECIL
    static const int KINDS[16]  = { 2, 2, 2, 2, 4, 4, 4, 9, 9, 2, 2, 2, 4, 4, 4, 2 };
    static const int IFAMS[16]  = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 2, 1, 1, 1, 1 };
    static const int NDIMS[16]  = { 2, 2, 2, 2, 1, 1, 2, 1, 1, 0, 1, 2, 2, 2, 1, 2 };
    static const int ISIZES[16] = { 8, 8, 8, 8,16,16,16, 1, 1, 8, 8, 8,16,16,16, 8 };
    static const int IGOTOS[16] = { 2, 2, 3, 3, 1, 1, 3, 4, 5, 6, 7, 8, 2, 2, 1, 3 };

    int dims[16];   // 0-based: dims[0]..dims[14]
    int IOB, LOCOB, LENOB;
    int KIND, KLASS, ISIZE, NCOLS, NROWS;
    int IFAM, NDIM, IGOTO;
    int LENOB2;
    int IER;
    int I;

    // Read OBNAME
    if (NXWORD(OBNAME) != 0) goto L900;

    // Read MEMBER (label 50)
L50:
    if (NXWORD(MEMBER) != 0) goto L900;
    if (std::memcmp(MEMBER, AS.data, 8) == 0) goto L50;

    // Convert MEMBER to lowercase
    for (int i = 0; i < 8; i++) {
        int ich = static_cast<unsigned char>(MEMBER[i]);
        if (ich >= 'A' && ich <= 'Z')
            MEMBER[i] = static_cast<char>(ich - ('A' - 'a'));
    }

    // Is it in the allocator?
    IOB = NAMLOC(reinterpret_cast<const char*>(OBNAME));
    LOCOB = Z[IOB];
    LENOB = LENG[IOB];

    // Set up default KIND, etc
    KIND  = 2;
    KLASS = 5;
    ISIZE = 8;
    NCOLS = 1;
    NROWS = 0;

    // See if we know something special about it
    for (I = 0; I < NSPECL; I++) {
        if (std::memcmp(SPECIL[I].data, OBNAME, 8) == 0) goto L200;
    }
    goto L500;

    // Process special stuff (label 200)
L200:
    KIND  = KINDS[I];
    IFAM  = IFAMS[I];
    NDIM  = NDIMS[I];
    ISIZE = ISIZES[I];
    IGOTO = IGOTOS[I];

    switch (IGOTO) {
        case 1: goto L500;
        case 2: goto L220;
        case 3: goto L230;
        case 4: goto L240;
        case 5: goto L250;
        case 6: goto L260;
        case 7: goto L270;
        case 8: goto L280;
    }

L220:
    NCOLS = NMLOLX;
    goto L500;

    // Something on the angle grid (label 230)
L230:
    NCOLS = static_cast<int>((ANGMAX - ANGMIN) / ANGSTP + 1.5);
    if (8 * LENOB / ISIZE == NROWS) NDIM = 1;
    goto L500;

    // HEADER (label 240)
L240:
    for (I = 1; I <= 65; I++) {
        if (HEADER[I] != BLANK) goto L245;
    }
    I = 65;
L245:
    {
        int dimval = 66 - I;
        char lib[1] = {' '};
        KEEP(&HEADER[I], 9, 1, 1, &dimval, lib, MEMBER, IER);
    }
    NROWS = 66 - I;
    goto L265;

    // REACTION (label 250)
L250:
    for (I = 1; I <= 45; I++) {
        if (REACT[46 - I] != BLANK) goto L255;
    }
    I = 45;
L255:
    {
        // Note: original Fortran has "9, 1 1, 46-I" — missing comma.
        // Translated as IFAM=1, NDIM=1 (matching the HEADER call pattern).
        int dimval = 46 - I;
        char lib[1] = {' '};
        KEEP(&REACT[1], 9, 1, 1, &dimval, lib, MEMBER, IER);
    }
    NROWS = 46 - I;
    goto L265;

    // ELAB (label 260)
L260:
    {
        int dimval = 1;
        char lib[1] = {' '};
        KEEP(reinterpret_cast<char*>(&ELAB), 2, 0, 0, &dimval, lib, MEMBER, IER);
    }
    NROWS = 1;

L265:
    NCOLS = 1;
    goto L800;

    // ANGLEGRID (label 270)
L270:
    {
        int dimval = 3;
        char lib[1] = {' '};
        KEEP(reinterpret_cast<char*>(&ANGSTP), 2, 1, 1, &dimval, lib, MEMBER, IER);
    }
    NROWS = 3;
    goto L265;

    // Object is a square matrix or array (label 280)
L280:
    NCOLS = static_cast<int>(DSQRT((8.0 * LENOB) / ISIZE + 0.010));
    goto L500;

    // Output an object in the allocator (label 500)
L500:
    if (IOB == 0) goto L910;
    LENOB = (8 * LENOB) / ISIZE;
    if (NROWS != 0) NCOLS = LENOB / NROWS;
    if (NCOLS != 0) NROWS = LENOB / NCOLS;
    LENOB2 = NROWS * NCOLS;
    if (LENOB == LENOB2) goto L600;

    printf("\n**** WARNING: %.8s HAS%6d"
           " ELEMENTS BUT NROWS, NCOLS & OUTPUT LENGHT =%6d%6d%6d"
           " KIND, ifam, ndim =%3d%3d%10d%10d\n",
           OBNAME, LENOB, NROWS, NCOLS, LENOB2,
           KIND, IFAM, NDIM, dims[0], dims[1]);

    NROWS = LENOB;
    NCOLS = 1;
    IFAM  = 1;
    NDIM  = 1;

    // Label 600
L600:
    if (NDIM == 2 && NROWS == 1) {
        NDIM  = 1;
        NROWS = NCOLS;
        NCOLS = 1;
    }
    dims[0] = NROWS;
    dims[1] = NCOLS;
    {
        char lib[1] = {' '};
        KEEP(reinterpret_cast<char*>(&ALLOC(LOCOB)), KIND, IFAM, NDIM, dims,
             lib, MEMBER, IER);
    }

    // Label 800 — check for errors
L800:
    if (IER != 0) goto L920;
    printf(" %.8s HAS BEEN KEPT AS MEMBER %.8s"
           "     KIND =%2d   family =%2d   num. dims =%2d   dims =%6d%6d\n",
           OBNAME, MEMBER, KIND, IFAM, NDIM, dims[0], dims[1]);
    return;

    // Error messages
L900:
    printf("\n**** KEEP MUST BE FOLLOWED BY TWO NAMES.\n");
    goto L980;

L910:
    printf("\n**** %.8s IS NOT A DEFINED OBJECT.\n", OBNAME);
    return;

L920:
    printf("\n**** ERROR NUMBER %08X ENCOUNTERED WHILE KEEPING %.8s AS %.8s\n",
           IER, OBNAME, MEMBER);
    return;

L980:
    printf(" **** CONTINUING WITH NEXT INPUT CARD.\n");
    NEWCD();
    return;
}
