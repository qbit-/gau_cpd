SHELL=/bin/csh
GAU_DIR = $(gdvroot)/gdv
GAU_DIRL = $(GAU_DIR)
GAU_DIRA = $(GAU_DIR)
UTIL_NAME = util.a

FLC = flc
#BLAS = -lcxml
#BLAS = -lblas
#BLAS = -lesslp2
BLAS = -lblas
IBM_FC = xlf_r -q64 -qsmp=noauto -qextname -qintsize=8 -qrealsize=8
SGI_FC = f77 -w -i8 -r8 -r8const -mips4 -64 -mp -r10000 -align64	\
-trapuv -Wl,-Xlocal,tstampc_ -Wl,-Xlocal,savedt_
SGI_FC2 = -lfastm
SUN_FC = f95 -fast -xtypemap=real:64,double:64,integer:64	\
-xtarget=native -xarch=v9 -xcache=generic
ALP_FC = f90 -O5 -transform_loops -omp -automatic -i8 -r8 -align	\
dcommons -tune host -trapuv -assume noaccuracy_sensitive -math_library	\
fast -reentrancy threaded
LIN_FC = pgf77 -mp -Mnostdlib
PGILIBS = $(PGI)/linux86/lib/libpgthread.a $(PGI)/linux86/lib/libpgc.a	\
$(PGI)/linux86/lib/libpgftnrtl.a $(PGI)/linux86/lib/libpgc.a		\
$(PGI)/linux86/lib/libpgmp.a $(PGI)/linux86/lib/libpgc.a
LIN_FC2 = /usr/local/lib/libf77blas.a /usr/local/lib/libatlas.a $(PGILIBS)
CRY_FC = f90 
CRY_FC2 = -Wl"-M /dev/null -D DUPENTRY=NOTE -D FORCE=OFF -f indef"
#CRY_FC2= -Wl,-z,muldefs
HP_FC = f77 +U77
HP_FC2 = -lveclib -lcps -lpthread -lm
PROFFLAG = 
#NUTILM = ../nutil/nutil-mark
#NUTIL = ../nutil/*.o
#NUTILL = ../nutill/*.lo ../nutill/*.o
#NUTILLM = ../nutill

FC0 = pgf90 -g $(FC) $(PROFFLAG)
FC2 = $(CRY_FC2)
FC1 = $(GAU_DIRL)/$(UTIL_NAME) $(FC2)

LAPext = ./lapack_LINUX.a

.SUFFIXES:
.SUFFIXES: .lo .o .F

.F.o:
	rm -f $*.o
	$(MAKE) -f $(GAU_DIRL)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.o

.F.lo:
	$(MAKE) -f $(GAU_DIRL)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.lo

all: l320.exe

#bdrys.o:
#	gau-get bdrys utilam
#	make -f $(GAU_DIRL)/bsd/gdv.make MAKE='$(MAKE)' \
#	PROFFLAG='$(PROFFLAG)' $*.o
#	rm -f bdrys.F

OBJUT = 
OBJPARUT =
OBJLIN = $(GAU_DIRL)/mdutil-linda.o $(GAU_DIRL)/drum.lo	\
	 $(GAU_DIRL)/caldsu.lo $(GAU_DIRL)/coulsu.lo	\
	 $(GAU_DIRL)/prnxtv.lo $(GAU_DIRL)/lindev.lo	\
	 $(GAU_DIRL)/chain.lo $(GAU_DIRL)/linint.lo	\
	 $(GAU_DIRL)/glini1.lo $(GAU_DIRL)/prlin3.lo	\
	 $(GAU_DIRL)/glinlo.lo $(GAU_DIRL)/prrfsu.lo	\
	 $(GAU_DIRL)/glinre.lo $(GAU_DIRL)/prsmsu.lo

$(COBJS):
	$(MAKE) -f $(GAU_DIRL)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.o


# = link 320 =

MAIN320 = ml320.o

OBJ320 = aobstf.o  binom.o   bldngg.o  bldpgg.o  cnvtst.o  detmat.o  \
	 detmt2.o  dmblck.o  erictr.o  evalhmt.o evals2.o  evalsmt.o \
	 evalyg.o  evlygg.o  fixmat.o  fixmt2.o  fockdg.o  focksel.o \
	 formfk.o  formgg.o  formng.o  formpg.o  formyg.o  gaufmt.o  \
	 gauleg.o  hfdmdg.o  hffock.o  iguess.o  indgrd.o  mltdg.o   \
	 phfcyc.o  phfdis.o  phfdm.o   phfdm1.o  phfdmx.o  phfdrv.o  \
	 prtdmt.o  prtocc.o  rtfact.o  solvci.o  sptblk.o  sptmat.o  \
	 trpgrd.o  trprct.o  wigarr.o  wignerd.o clgord.o

l320.exe: $(MAIN320) $(OBJ320) $(NUTILM)
	$(FC0) -g -o l320.exe $(MAIN320) $(OBJ320) $(NUTIL) \
	$(LAPext) $(FC1) $(BLAS)
	chmod o-rx l320.exe

ck320:
	cat $(MAIN320:.o=.F) $(OBJ320:.o=.F) $(OBJUT:.o=.F) \
	$(NUTIL:.o=.F) >x.x
	checkf x.x x
	rm -f x.x



