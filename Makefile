SHELL=/bin/csh
GAU_DIR = $(gdvroot)/gdv
GAU_DIRL = $(GAU_DIR)
GAU_DIRA = $(GAU_DIR)
NUTIL=$(GAU_DIR)/util.a

MKLPATH = /opt/intel/composerxe-2015.0.090/mkl/lib/intel64
MKLINCLUDE = /opt/intel/composerxe-2015.0.090/mkl/include
PROFFLAG =
FC0 = pgf77
#FC2 = -Wl"-M /dev/null -D DUPENTRY=NOTE -D FORCE=OFF -f indef"

#LAPext = /usr/lib64/libreflapack.so
LAPext = -L$(MKLPATH) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 

.SUFFIXES:
.SUFFIXES: .lo .o .F

.F.o:
	rm -f $*.o
	$(MAKE) -f $(GAU_DIR)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.o

.F.lo:
	$(MAKE) -f $(GAU_DIR)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.lo

all: l325.exe
	gdv run/test.gjf && tail -40 run/test.log

bdrys.o:
	gau-get bdrys utilam
	make -f $(GAU_DIR)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.o
	rm -f bdrys.F

# = link 325 =

MAIN325 = ml325.o

OBJ325 =  atquadwrt.o cpd2int.o ctrmemest.o gtinct.o frm2eints.o frmemt.o  frmquad.o frmspovinv.o

l325.exe: $(MAIN325) $(OBJ325) 
	$(FC0) -g -o l325.exe $(MAIN325) $(OBJ325) $(NUTIL) \
	$(LAPext) $(FC1) $(BLAS)
	chmod o-rx l325.exe

ck325:
	cat $(MAIN325:.o=.F) $(OBJ325:.o=.F) $(OBJUT:.o=.F) \
	$(NUTIL:.o=.F) >x.x
	checkf x.x x
	rm -f x.x


