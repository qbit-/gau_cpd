#!/bin/sh
PGI=/opt/pgi
export PGI
PATH=/opt/pgi/linux86-64/13.3/bin:$PATH
export PATH
MANPATH=$MANPATH:/opt/pgi/linux86-64/13.3/man
export MANPATH
gdvroot=/home/shiva/software/gdv-h21
GAUSS_MEMDEF=67108864
GAUSS_SCRDIR=/tmp
export gdvroot GAUSS_MEMDEF GAUSS_SCRDIR
. $gdvroot/gdv/bsd/gdv.profile
export MKL_NUM_THREADS=4

alias mk=" make  GAU_DIR=/home/shiva/software/gdv-h21/gdv BLAS='' UTIL_NAME='{util,archlib,mdarch,util,archlib,mdarch,bsd/libf77blas-corei764sse3,bsd/libatlas-corei764sse3}.a' MACHTY=nehalem-64 GAULIBU=util.a I8FLAG=-i8 R8FLAG=-r8 MMODEL='-mcmodel=medium' OPTOI= I8CPP1=-DI64 I8CPP2=-DP64 I8CPP3=-DPACK64 I8CPP4=-DUSE_I2 NJSEC=-DDEFJSEC=512 X86TYPE=-DX86_TYPE=S GAUDIM=2500 FCN='pgf77 -Bstatic_pgi' FC='-mp=nonuma -tp nehalem-64 -mcmodel=medium' FC2='-lpthread -lm -lc'"
