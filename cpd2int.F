

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
#include "commonb2.inc"
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
      integer     iopcl

C     general options variables

      integer     iprint
      integer     iradan, iranwt, irangd, iradanu, iranwtu, irangdu
      
C     memory allocation variables

      integer     jend, lenv
      
C     other variables

      integer     nbsq
      integer     ngrid, nbas6d   
      real*8      tolinv, accdes

C     CPD2E RWF files

      integer     irwfac, irwfact, irwspo, irwz, irwe

      parameter   ( irwfac = 824, irwfact = 825, irwspo = 826 )
      parameter   ( irwe = 827, irwz = 828)

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
      
C    - neq, neqshl, neqsh2 - files useful for symmetry in PRISM

      integer     neq, neqshl, neqsh2
      parameter   ( neq = 580, neqshl = 565, neqsh2 = 726 )

C     symmetry related variables

      integer     jsym2e

C     functions called

      real*8      decacc
      integer     intowp

C     debug temporary stuff
      
      integer     dims(4)
      logical     dbg

C------------RI testing-------------------------------------------------

      integer     irwria, irwrib, ir12
      integer     jovoff, jov, j2eri, j2erit
      integer     ndbf, ndbf6d, nbas6d
      integer     j2erec
      real*8      XX
      data        XX/0.0d0/

      integer     jtmp, next, irw2efull, ii
      integer     janb, nshelb, maxtyb
C     format statements

 1001 format (' CPD2E was invoked with the following options:')
 1010 format (' iprint = ', I2)
 1011 format ('   grid selection:')      
 1012 format (' iradan = ', I6, ' iranwt = ', I5, ' irangd = ', I5)
 1020 format (' accuracy = ',D8.3, ' pseudoinv. tolerance = ', D8.3)
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

C----------------------------------------------------------------------+
C       Process options for CPD2E                                      |
C----------------------------------------------------------------------+


      write (iout, *)
      write (iout, 1001)
      write (iout, *)


C     Set printing level.

      iprint = iop(6)

      if ( iprint .lt. 0 .or. iprint .gt. 10 ) iprint = 0
      write (iout, 1010) iprint
      
C     Set up grid options 
C     Read grid parameters from ilsw, but substitute 
C     with user supplied if available
      
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

C     Desired accuracy and pseudoinversion tolerance for spatial overlap

      accdes = decacc (0, 0.0d0, .false.)
      tolinv = 1.0e-6
      write (iout, 1020) accdes, tolinv
      
C     Use of symmetry
C     We currently set jsym2e = 2 in the code (users cannot change that).

C       jsym2e = 2 uses integral symmetry by replicating integrals, but
C       does not force the density or Fock matrices to be symmetric.

      jsym2e = 0

C     Restore /IBF/

      call fileio (2, -irwibf, intowp (30), ibf, 0)

      if ( isym2e .eq. 0 ) jsym2e = 0
      write (iout, 1040)

      if ( jsym2e .eq. 2 ) write (iout, 1040)

C     Read iopcl from ILSW file.
C       iopcl = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)
      
      call ilsw (2, 1, iopcl)

      write (iout, *) 'got iopcl = ', iopcl
C     call gauerr ('Incorrect iopcl in cpd2e.')
                    
C     Read /mol/

      call rwmol(2,997)

C     Get 6d basis used with all quadratures

      call getnb6(nbas6d)

      write (iout, *)
      
C-----------------------------------------------------------------------------+
C     Compute the collocation matrix X(mu,G) and save it to the RWF irwfac    |
C-----------------------------------------------------------------------------+

      call frmquad (iout, iprint, iradan, iranwt, irangd, iopcl,
     & ngrid, irwfac, irwfact, v(jend), lenv - jend)

C-----------------------------------------------------------------------------+      
C     Compute spatial overlap using collocation from irwfac                   |
C-----------------------------------------------------------------------------+      

      call frmspovinv (iout, iprint, nbas6d, ngrid, tolinv,
     &     irwfac, irwspo, v(jend), lenv - jend)
      
C-----------------------------------------------------------------------------+      
C     Project integrals onto the quadrature space                             |
C-----------------------------------------------------------------------------+

      call frmemt(iout, iprint, nbas6d, ngrid, jsym2e, neq,
     &     neqshl, neqsh2, accdes, irwfac, irwe, v(jend), lenv - jend)

C-----------------------------------------------------------------------------+
C     Contract projected integrals with the spatial overlap                   |
C-----------------------------------------------------------------------------+

      call frmz(iout, iprint, ngrid, irwspo, irwe, irwz, v(jend), lenv
     & - jend) 

      dbg = .true. 
      if (dbg) then
         call fileio(2,-irwfac,nbasis*ngrid,v(jend),0)
         dims(1) = nbasis
         dims(2) = ngrid
         call outoctfmt(iout, 'x.mat', 'x', 2, dims, v(jend))

         call fileio(2,-irwspo,ngrid*ngrid,v(jend),0)
         dims(1) = ngrid
         dims(2) = ngrid
         call outoctfmt(iout, 's.mat', 's', 2, dims, v(jend))
         
         call fileio(2,-irwe,ngrid*ngrid,v(jend),0)
         dims(1) = ngrid
         dims(2) = ngrid
         call outoctfmt(iout, 'e.mat', 'e', 2, dims, v(jend))
         
         call fileio(2,-irwz,ngrid*ngrid,v(jend),0)
         call outoctfmt(iout, 'z.mat', 'z', 2, dims, v(jend))         
      endif

C-----------------------------------------------------------------------------------

C     Load /B2/ with the density basis.
      Call RdWrB(-10,0)
      call getnb6(nbas6d)
      
C     Figure out NDBF and NDBF6D.
      call GetMBN(I5DB2,I7FB2,NShelB,ShelTB,ShelCB,NDBF)
      call GetMBN(1,1,NShelB,ShelTB,ShelCB,NDBF6D)

      ir12 = 1

      call frm2eri (iout, iprint, jsym2e, neq, neqshl, 
     &     neqsh2, ir12, .true., .true., tolinv, accdes, jovoff,
     &     next, irwria, irwrib, v(jend), lenv - jend)

      jov     = jend
      j2eri   = jov    + jovoff
      j2erit  = jend   + next
      jtmp    = j2erit + nbas6d*nbas6d*ndbf
      j2erec  = jtmp   + nbas6d*nbas6d*ndbf
      jend    = j2erec + nbas6d*nbas6d*nbas6d*nbas6d
       
      call TrCor1(1,1,ndbf,nbas6d*nbas6d,3,V(j2eri),XX,XX,V(j2erit),XX,
     &     XX,.False.,XX,XX)
      call dgemm ('n', 'n', nbas6d*nbas6d, ndbf, ndbf, 1.0d0, v(j2eri),
     &     nbas6d*nbas6d, v(jov), ndbf, 0.0d0, v(jtmp), nbas6d*nbas6d )
      call dgemm ('n', 'n', nbas6d*nbas6d, nbas6d*nbas6d, ndbf, 1.0d0,
     &     v(jtmp), nbas6d*nbas6d, v(j2erit), ndbf, 0.0d0,
     &     v(j2erec), nbas6d*nbas6d)

      call frm2eints (iout, iprint, jsym2e, neq, neqshl, 
     &     neqsh2, accdes, .true., next, irw2efull, v(jend),lenv- jend)
      
      call daxpy (nbas6d*nbas6d*nbas6d*nbas6d, -1.0d0, v(j2erec), 1,
     &     v(jend), 1)

      write (iout,'10E12.3') (v(jend + ii),
     &     ii=0,nbas6d*nbas6d*nbas6d*nbas6d-1)

      call frmemtri (iout, iprint, nbasis, ngrid, jsym2e, neq, neqshl,
     &     neqsh2, ir12, tolinv, accdes, irwfac, irwe, v(jend),
     &     lenv - jend)

      call fileio(2, -irwe, ngrid*ngrid, v(jend), 0)
      dims(1) = ngrid
      dims(2) = ngrid
      call outoctfmt(iout, 'eri.mat', 'eri', 2, dims, v(jend))
      
      call chainX (0)

      return
      end


