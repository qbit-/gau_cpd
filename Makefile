SHELL=/bin/csh
GAU_DIR = $(gdvroot)/gdv
GAU_DIRL = $(GAU_DIR)
GAU_DIRA = $(GAU_DIR)
NUTIL=$(GAU_DIR)/util.a

MKLROOT=/opt/intel/composerxe-2013.3.174/mkl
INTELROOT=/opt/intel/composerxe-2013.3.174/compiler
INTELPATH = $(INTELROOT)/lib/intel64
MKLPATH = $(MKLROOT)/lib/intel64
MKLINCLUDE = -I$(MKLROOT)/include -I$(INTELROOT)/include
MKLLIB = -L$(MKLPATH) -L$(INTELPATH) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm 
PROFFLAG = #-Mprof=func
#OPTFLAG  = -O3 -mp 
OPTFLAG   = -O0 -g 
#FC0      = pgf77 -i8 -r8 -mp -O3
FC0      = pgf77 -i8 -r8 -O0 -g
#FC0 = pgf90 -Wl,-z,muldefs

.SUFFIXES:
.SUFFIXES: .lo .o .F

.F.o:
	rm -f $*.o
	$(MAKE) -f $(GAU_DIR)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' OPTFLAG='$(OPTFLAG)' $*.o

.F.lo:
	$(MAKE) -f $(GAU_DIR)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' OPTFLAG='$(OPTFLAG)' $*.lo

all: l325.exe
	gdv run/testri.gjf
	tail -40 run/testri.log

# = link 325 =

MAIN325 = ml325.o

OBJ325  = aclearf.o atquadwrt.o cpd2int.o convquad.o ctrmemest1.o ctrmemest_nd.o ctrmemest2.o frm2eints.o frm2eri.o frmemt.o frmemt_nd.o frmemtri.o frmquad.o frmextquad.o frmria.o frmrib.o frmspovinv.o \
	  frmz.o updmatf.o
TEMPOBJ = asubf.o cpdfock.o cpdfkmem.o cpdexmem.o ctrcou.o ctrexc.o frmspovinv_blas.o normfro.o outcsv.o outoctfmt.o reconstr.o symmetric.o toeplitz.o \
	  readfmt.o  localize.o moquadwrt.o 

OBJGAU     = geninv.o genin1.o
TEMPOBJGAU = 

all:     l325.exe l330.exe

l325.exe: $(MAIN325) $(OBJ325) $(TEMPOBJ) $(OBJGAU) $(TEMPOBJGAU)
	$(FC0) $(OPTFLAG) $(PROFFLAG) $(MKLINCLUDE) -o l325.exe $(MAIN325) $(OBJ325) $(TEMPOBJ) $(OBJGAU) $(TEMPOBJGAU) $(NUTIL) \
	$(MKLLIB)
	chmod o-rx l325.exe

ck325:
	cat $(MAIN325:.o=.F) $(OBJ325:.o=.F) $(OBJUT:.o=.F) \
	$(NUTIL:.o=.F) >x.x
	checkf x.x x
	rm -f x.x

# = link 330 =

MAIN330 = ml330.o

OBJ330  = dump_all.o unpack_t2.o

l330.exe: $(MAIN330) $(OBJ330) $(OBJ325) $(TEMPOBJ)
	$(FC0) $(OPTFLAG) $(PROFFLAG) $(MKLINCLUDE) -o l330.exe $(MAIN330) $(OBJ330) $(OBJ325) $(TEMPOBJ) $(OBJGAU) $(TEMPOBJGAU) $(NUTIL) \
	$(MKLLIB)
	chmod o-rx l330.exe

