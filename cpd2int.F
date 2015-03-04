

      subroutine cpd2int (v, ngot)

      implicit none


C +------------------------------------------------------------+
C |  CPD-2INT   ---   the Scuseria group, 02.2015              |
C |                                                            |
C |                                                   v0.1 -   |
C |                                                            |
C |    ( link 320, part of the Gaussian system of programs )   |
C |                                                            |
C |                                                            |
C |  This link performs canonical decomposition of the 2       |
C |  electron integrals as described in                        |
C |                                                            |
C |    http://scuseria.rice.edu/                               |
C |                                                            |
C |                                                            |
C |  IMPORTANT NOTE:                                           |
C |                                                            |
C |    We do not guarantee that this software is bug-free.     |
C |    The use of this code is sole responsibility of the      |
C |    user.                                                   |
C |                                                            |
C +------------------------------------------------------------+

C +------------------------------------------------------------+
C |                                                            |
C |  Description of IOps                                       |
C |  ===================                                       |
C |                                                            |
C |                                                            |
C |  iop(6) - printing level                                   |
C |     = 0, [default] print results at convergence            |
C |     = 1, print energy at every iteration                   |
C |     = 2, print CI matrices at every iteration              |
C |     = 3, print initial and converged density matrices      |
C |     = 4, print density matrix at every iteration           |
C |     = 5, print (gauge angle) matrices at every iteration   |
C |     ...                                                    |
C |     = 10, print two-electron integrals                     |
C |                                                            |
C +------------------------------------------------------------+


C +------------------------------------------------------------+
C |  cpd2int   ---   RSS, 02.2015                              |
C |                                                            |
C |  Driver routine for CPD decomposition of 2 electron        |
C |  integrals                                                 |
C |                                                            |
C +------------------------------------------------------------+

C     include some Gaussian files
      
#include "commonmol.inc"
#include "commonb.inc"
#include "repall.inc"

C     input / output variables

      real*8      v(*)
      integer     ngot

C     iop common block

      integer     iop, maxiop
      parameter   ( maxiop = 200 )
      common      /iop/ iop(maxiop)

C     other common blocks

      integer     in, iout, ipunch
      integer     ibf, isym2e

      common      /io/ in, iout, ipunch
      common      /ibf/ ibf(30)

      equivalence (isym2e, ibf(30))

C     general array

      integer     lengen
      parameter   ( lengen = 55 )
      real*8      dgen(lengen)


C     prism variables

      integer     ifmm, ipflag, fmflag, fmflg1, nfxflg, ihmeth, iseall
      integer     lseall, momega, nomega
      logical     allowp(50), fmm
      parameter   ( momega = 5, nomega = 6 )
      real*8      omega(momega,nomega)


C     general options variables

      integer     iprint
      integer     iradan, iranwt, irangd, iradanu, iranwtu, irangdu
      
C     memory allocation variables

      integer     jend, lenv, mdv, junk
      real*8      xx

C     other variables
      
      integer     iopcl, nbas6d
      integer     ngridmx, ntheta, nphi, mxradu, mxangu, ngridu, iprune
      integer     iopwrt, leniopwrt, lenscrwrt

      parameter  ( leniopwrt = 5 )
      dimension   iopwrt(leniopwrt)
      
C     CPD2E RWF files

      integer     irwfac, irwspo, irwz, irwe

      parameter   ( irwfac = 825, irwspo = 826, irwz = 827, irwe = 828 )

C     Gaussian RWF files
C       - irwgen - general array (see Link 1 for details)
C       - irwibf - file with IBF array
C       - irweig - MO orbital energies
C       - irwc?  - MO (x-spin) coefficients
C       - irwp?  - density matrix (x-spin)
C       - irwmol - /mol/ common
      
      integer     irwgen, irwibf, irweig
      integer     irwca, irwcb, irwx, irwpa, irwpb, irws, irwh
      integer     irwmol
      
      parameter   ( irwgen = 501, irwibf = 508, irweig = 522 )
      parameter   ( irwca = 524, irwcb = 526, irwpa = 528, irwpb = 530 )
      parameter   ( irws = 514, irwx = 685, irwh = 515 )
      parameter   ( irwmol = 997 )
      
C       - neq, neqshl, neqsh2 - files useful for symmetry

      integer     neq, neqshl, neqsh2

      parameter   ( neq = 580, neqshl = 565, neqsh2 = 726 )


C     symmetry related variables

      integer     jsym2e
      integer     nop1, nop2, nopuse
      integer     lenneq, lennes, lenne2
      integer     ineq, ineqsh, ineqs2


C     functions called

      real*8      scfacc, decacc
      integer     intowp, itqry
      external    atquadwrt, mofunx, denfnx 
      

C     format statements

 1001 format (' CPD2E was invoked with the following options:')
 1010 format (' iprint = ', I2)
 1011 format ('   grid selection:')      
 1012 format (' iradan = ', I6, ' iranwt = ', I5, ' irangd = ', I5)
 1020 format ('   iop')
 1021 format (' exact integrals will be used ')
 1022 format (' density fitted integrals will be used ')
 1040 format (' integrals will be replicated using symmetry')
      write (iout, *) ' '
      write (iout, *) ' +-------------------------------------------+'
      write (iout, *) ' |                                           |'
      write (iout, *) ' |  LINK 320 - Canonical decomposition of 2e |'
      write (iout, *) ' |     integrals                             |'
      write (iout, *) ' |                                           |'
      write (iout, *) ' |                                           |'
      write (iout, *) ' |         the Scuseria group, v0.1, 02.2015 |'
      write (iout, *) ' +-------------------------------------------+'
      write (iout, *) ' '

      call drum (v, ngot)
      jend = 1
      lenv = ngot

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Options for CPD2E            %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      write (iout, *)
      write (iout, 1001)
      write (iout, *)


C     Set printing level.

      iprint = iop(6)

      if ( iprint .lt. 0 .or. iprint .gt. 10 ) iprint = 0
      write (iout, 1010) iprint
      
C     Set up grid options 
C     =======================

C     Read grid options from ilsw
      
      call ilsw (2, 13, iradan)      
      call ilsw (2, 35, iranwt)
      call ilsw (2, 36, irangd)
      
      iradanu = iop(60)
      if ( iradanu .ne. 0 ) iradan = iradanu

      iranwtu = iop(61)
      if ( iranwtu .ne. 0 ) iranwt = iranwtu

      irangdu = iop(62)
      if ( irangdu .ne. 0 ) irangd = irangdu

      write (iout, 1011)
      write (iout, 1012) iradan, iranwt, irangd

C     Use of symmetry
C     ===============

C     We currently set jsym2e = 2 in the code (users cannot change that).

C       jsym2e = 2 uses integral symmetry by replicating integrals, but
C       does not force the density or Fock matrices to be symmetric.

      jsym2e = 0

C     Load IBF array.

      call fileio (2, -irwibf, intowp (30), ibf, 0)

      if ( isym2e .eq. 0 ) jsym2e = 0


      write (iout, 1040)

      if ( jsym2e .eq. 2 ) write (iout, 1040)

      write (iout, *)

C     Read /mol/

      call rwmol(2,997)

C     Get 6d basis

      call getnb6(nbas6d)

C     Read iopcl from ILSW file.

C       iopcl = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)
      
      call ilsw (2, 1, iopcl)

      if ( iopcl .gt. 6 .or. iopcl .lt. 1 ) then
         write (iout, *) 'iopcl', iopcl
C        call gauerr ('Incorrect iopcl in cpd2e.')
              
      endif

      
C     Get number of symmetry operations.
C       - nop1,    order of the concise abelian subgroup
C       - nop2,    order of the full abelian subgroup
C       - nopuse,  actual oder to be used

C     %%%%%%%%%%%%%%%%%%%%%%%
C     %  Memory allocation  %
C     %%%%%%%%%%%%%%%%%%%%%%%

C     Initialize CPD2E RWF files.

C      call conddf (irwx, szp)

C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Retrieve symmetry information  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      call getnop (nop1, nop2)
      call pccki (0, junk, neqall, mxatso, nopall, natoms, nopall)

      nopuse = nopall

C     Memory allocation for symmetry related quantities.

      lenneq = 0
      lennes = 0
      lenne2 = 0

      if ( isym2e .eq. 1 ) then
        lenneq = intowp (nop1*nbasis)
        lennes = max (itqry (neqshl), 0)
        lenne2 = max (itqry (neqsh2), 0)
      endif

      ineq   = jend
      ineqsh = ineq   + lenneq
      ineqs2 = ineqsh + lennes
      jend   = ineqs2 + lenne2

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'cpd2e')


C     Fill arrays with symmetry related quantities.

      if ( lenneq .ne. 0 ) then
        call fileio (2, -neq, lenneq, v(ineq), 0)
      endif

      if ( lennes .ne. 0 ) then
        call fileio (2, -neqshl, lennes, v(ineqsh), 0)
      endif

      if ( lenne2 .ne. 0 ) then
        call fileio (2, -neqsh2, lenne2, v(ineqs2), 0)
      endif


C     Compute the collocation matrix X(mu,G) and save it to the RWF irwfac
      
      call DecPrn(iout,iprint,iradan,0,iprune,ntheta,nphi,
     &  mxradu,mxangu,xx) 
      ngridmx = mxradu * mxangu * natoms

      call conddf (irwfac, ngridmx*nbas6d)

C     Prepare the arguments passed/received from the writing routine

      lenscrwrt  = lenv / ngridmx
      
C     iopwrt integer array:
C        1 - RWF to store X(mu,G)
C        2 - total number of grid points (incremented by working routine)
      
      iopwrt(1) = irwfac
      iopwrt(2) = 0
      
C      call NIn3Su(iout,iprint,iopcl,ParOK,ParDig,NPMax,ReadB,           
C     &  UseBfn,UseMO,UseOcc,IUseDn,IDerBf,IDerRh,IRadAn,IRanWt,IRanGd,        
C     &  IRType,NMO,NAtoms,NBas6D,IAn,IAtTyp,AtmChg,C,AtNetC,CMO,OccMO,
C     &  PA,PB,HvRMax,RQMax,BfnFnc,IntBfn,WPBfn,LScrBf,MOFunc,IntMO,WPMO,
C     &  LScrMO,DenFnc,IntDen,WPDen,LScrDn,V,MDV)

      call NIn3Su(iout,iprint,iopcl,.false.,.true.,1,.true.,
     &  .true.,.false.,.false.,0,0,0,iradan,iranwt,irangd,
     &  2,junk,natoms,nbas6d,ian,iattyp,atmchg,c,xx,xx,xx,
     &  xx,xx,.false.,xx,atquadwrt,iopwrt,xx,lenscrwrt,MOFunX,junk,
     &  xx,junk,DenFnX,junk,xx,junk,v,lenv)

C     Read back the output from the writing routine
      ngridu = iopwrt(2)
      
      write (iout, *) 'nin3su completed, npoints = ', ngridu
      call chainX (0)

      return
      end


