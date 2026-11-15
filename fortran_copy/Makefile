FC      = gfortran
CC      = gcc
EXEC    = ptolemy.x

# -m32, -fdefault-real-8 and -O0 are all required; see compile.sh for why.
FFLAGS  = -m32 -g -O0 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -fallow-argument-mismatch -fno-range-check -fallow-invalid-boz
CFLAGS  = -m32 -std=gnu89 -Wno-implicit-function-declaration

FSRC    = av18.f fortlib.f gfortran_stuff.f keep.f keepsub.f keptsub.f \
          linkule.f linkulesfitters.f masstable.f numbered_store.f phiffer.f \
          ptolemy-main.f source.f

CSRC    = dtime.c srread.c

FOBJ    = $(FSRC:.f=.o)
COBJ    = $(CSRC:.c=.o)

all: $(EXEC)

$(EXEC): $(FOBJ) $(COBJ)
	$(FC) -m32 -g -o $@ $^

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o $(EXEC)

.PHONY: all clean
