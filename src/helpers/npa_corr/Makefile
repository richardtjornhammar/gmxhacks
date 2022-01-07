# Generated automatically from Makefile.in by configure.
#
# This is a Gromacs 3.0 template makefile for your own utility programs.
#
# Copy this file to whatever directory you are using for your own
# software and add more targets like the template one below.
#
# If you are using gmake it is relatively straightforward to add
# an include based on environment variables (like previous Gromacs versions)
# to select compiler flags and stuff automatically, but below it is static:
#

# Variables set by the configuration script:
LIBS         = -lgmx -lmd -lgmxana -lnsl -lfftw3f -lm   -lSM -lICE -lX11   
LDFLAGS      = -L/home/richardt/prg/gmx_mpqc/lib   
CFLAGS	     = -O3 -fomit-frame-pointer -finline-functions -Wall -Wno-unused -funroll-all-loops   -I/home/richardt/prg/gmx_mpqc/include -I/home/richardt/prg/gmx_mpqc/include/gromacs	
CC           = cc
LD           = $(CC)
GMXLIB       =/home/richardt/prg/gromacs-4.0.7/src/gmxlib

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...
OBJ          = npa_corr_large.o 

#g_dipole:       tools.o template.o
npa_autocor:    $(OBJ)
		$(LD) $(LDFLAGS) -o $@ $(GMXLIB)/xvgr.o $(OBJ) $(LIBS)
