SRCS := $(wildcard *.f)
OBJS = $(SRCS:.f=.o)
INCS = $(wildcard *.inc)

LIBS=-lvffpack

# For Pentium IV
FC=ifort
FFLAGS=-xP -tpp7 -ftz -fpe0 -O3 -ip -parallel
LIBS=-lvffpack

# Intel Compiler for any x86 machine
# FC=ifort
# FFLAGS=-ftz -fpe0 -O3 -ip -parallel

# Debugging for Pentium IV
# FC=ifort
# FFLAGS=-arch pn4 -xN -tpp7 -tune pn4 -ftz -fpe0 -g

# g77 for Pentium IV
# FC=g77
# FFLAGS=-march=pentium IV -O3 -mfpmath=sse -fomit-frame-pointer -ffast-math

# g77 for AMD Athlon MP
# FC=g77
# FFLAGS=-march=athlon-mp -O3 -mfpmath=sse -fomit-frame-pointer -ffast-math

# g77 debugging for AMD Athlong MP
# FC=g77
# FFLAGS=-march=athlon-mp -g

iim2d: $(SRCS) $(INCS) makefile
	$(FC) -o $@ $(FFLAGS) $(SRCS) $(LIBS)
	strip $@

clean:
	rm -rf *.o
