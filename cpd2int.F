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

C     general option variables

      integer     iprint
      integer     iradan, iranwt, irangd, iradanu, iranwtu, irangdu
      integer     iuseri, ir12

      logical     ri

      integer     iradandef, iranwtdef, irangddef
      parameter   ( iradandef = 8012, iranwtdef = -1, irangddef = 0 )
      
C     memory allocation variables

      integer     jend, lenv
      
C     other variables

      integer     nbsq
      integer     ngrid, nbas6d
      integer     nrank, itol   
      real*8      accdes

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

      real*8      decacc, gfloat
      integer     intowp

C     debug temporary stuff
      
      integer     dims(4)
      integer     jj, irw2efull, irwrec
      integer     irwria, irwrib
      integer     ndbf6d
      real*8      normfro, froful, frodiff
      real*8      outbuf(4)
      logical     dbg, tests

      logical     neg
      integer     jrho
      real*8      xx
      
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

C     Defaults if nothing is set
      if ( iradan .eq. 0 ) iradan = iradandef
      if ( iranwt .eq. 0 ) iranwt = iranwtdef
      if ( irangd .eq. 0 ) irangd = irangddef
      
      write (iout, 1011)
      write (iout, 1012) iradan, iranwt, irangd

C     Desired accuracy and pseudoinversion tolerance for spatial overlap

      accdes = decacc (0, 0.0d0, .false.)
      itol = 10
      write (iout, 1020) accdes, Gfloat(10)**(-itol)
      
C     Use of symmetry
C     We currently set jsym2e = 2 in the code (users cannot change that).

C       jsym2e = 2 uses integral symmetry by replicating integrals, but
C       does not force the density or Fock matrices to be symmetric.

      jsym2e = 0

C     Use of RI (N^4 scaling algorithm) and related options

      iuseri = mod(iop(109), 10)
      if(iuseri .eq. 1) then
         ri = .true.
         write (iout, 1022)
      else
         ri = .false.
      endif
      
C     Metric to use with RI: 1 - 1/r, 2 - overlap      
      ir12 = mod (iop(109), 100000) / 10000
      if (ir12 .eq. 0) ir12 = 1
      
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
      if (nbasis .ne. nbas6d)
     &  call gauerr ('Sperical bases are not supported yet')
      
      write (iout, *)
      
C------------------------------------------------------------------------------+
C     Fixme: frmquad breaks PRISM settings, hence integral formation done first
C------------------------------------------------------------------------------+     

c$$$      if (.not. ri) then
c$$$         call frm2eints (iout, iprint, jsym2e, neq, neqshl, neqsh2,
c$$$     &     accdes, .false., jj, irw2efull, v(jend), lenv)
c$$$
c$$$      else
c$$$         call frm2eri (iout, iprint, jsym2e, neq, neqshl, 
c$$$     &        neqsh2, ir12, .false., .false., itol, accdes, jj,
c$$$     &        jj, irwria, irwrib, v(jend), lenv - jend)
c$$$
c$$$      endif

C-----------------------------------------------------------------------------+
C     Compute the collocation matrix X(mu,G) and save it to the RWF irwfac    |
C-----------------------------------------------------------------------------+

      call frmquad (iout, iprint, iradan, iranwt, irangd, iopcl,
     & ngrid, irwfac, irwfact, v(jend), lenv - jend)

C-----------------------------------------------------------------------------+      
C     Compute spatial overlap using collocation from irwfac                   |
C-----------------------------------------------------------------------------+      

      call frmspovinv (iout, iprint, nbasis, ngrid, itol,
     &     irwfac, nrank, irwspo, v(jend), lenv - jend)
      
C-----------------------------------------------------------------------------+      
C     Project integrals onto the quadrature space                             |
C-----------------------------------------------------------------------------+

      if (.not. ri) then
         call frmemt_nd(iout, iprint, nbasis, ngrid, jsym2e, neq,
     &        neqshl, neqsh2, accdes, irwfac, irw2efull, irwe,
     &        v(jend), lenv - jend)
         
      else
         call frmemtri(iout, iprint, nbasis, ngrid, jsym2e, neq, neqshl,
     &        neqsh2, ir12, itol, accdes, ndbf6d, irwfac,
     &        irwria, irwrib, irwe, v(jend), lenv - jend)
      endif
      
C-----------------------------------------------------------------------------+
C     Contract projected integrals with the spatial overlap                   |
C-----------------------------------------------------------------------------+

      call frmz(iout, iprint, ngrid, irwspo, irwe, irwz, v(jend), lenv
     & - jend) 

C----------------------------------------------------------------------------+
C     Testing and staging stuff                                              |
C----------------------------------------------------------------------------+
      
      dbg = .true. 
      if (dbg) then
         call fileio(2,-irwfac,nbasis*ngrid,v(jend),0)
         dims(1) = nbasis
         dims(2) = ngrid
         call outoctfmt(iout, 'x.mat', 1, 'x', 2, dims, v(jend))

         call fileio(2,-irwspo,ngrid*(ngrid+1)/2,v(jend),0)
         call square (v(jend), v(jend), ngrid, ngrid, 0)
         dims(1) = ngrid
         dims(2) = ngrid
         call outoctfmt(iout, 's.mat', 1, 's', 2, dims, v(jend))
         
         call fileio(2,-irwe,ngrid*ngrid,v(jend),0)
         dims(1) = ngrid
         dims(2) = ngrid
         call outoctfmt(iout, 'e.mat', 1, 'e', 2, dims, v(jend))
         
         call fileio(2,-irwz,ngrid*ngrid,v(jend),0)
         call outoctfmt(iout, 'z.mat', 1, 'z', 2, dims, v(jend))         
      endif

      tests = .true.
      if (tests) then
         call frm2eints (iout, iprint, jsym2e, neq, neqshl, neqsh2,
     &     accdes, .false.,.true., jj, irw2efull, v(jend), lenv)
         
         froful=normfro(iout, iprint, nbasis, irw2efull, v(jend), lenv)
         
         call reconstr (iout, iprint, nbasis, ngrid,
     &     irwfac, irwz, irwrec, v(jend), lenv)
         call asubf(iout,iprint,nbasis,irw2efull,irwrec, v(jend), lenv) 
         frodiff=normfro (iout, iprint, nbasis, irwrec, v(jend), lenv)

 1100    format(A,I20)
 1200    format(A,E20.12)
         
         outbuf(1) = dble(nbasis)
         outbuf(2) = dble(ngrid)
         outbuf(3) = froful
         outbuf(4) = frodiff
         
         call outcsv (iout, 'summary.csv', 2, 4, outbuf)
                  
         write (iout, 1100) ' nbasis = ', nbasis
         write (iout, 1100) ' ngrid = ', ngrid
         write (iout, 1100) ' nrank = ', nrank
         write (iout, 1200) ' full = ', froful
         write (iout, 1200) ' diff = ', frodiff
         write (iout, 1200) ' frac % = ', frodiff / froful * 100
         if (ri) write (iout, 1100) ' ndbf = ', ndbf6d
      endif

      jrho = jend
      jend = jrho + nbasis*nbasis
      call toeplitz (iout, nbasis, 0, v(jrho))
      dims(1) = nbasis
      dims(2) = nbasis
      call outoctfmt(iout, 'rho.mat', 1, 'rho', 2, dims, v(jrho))         
      
      call linear (v(jrho), v(jrho), nbasis, nbasis)
      
      call cpdfock (iout, iprint, 1, .true., nbasis, ngrid, irwfac,
     &     irwz, 1, 0, 1, v(jrho), v(jrho), v(jend),
     &     v(jend+nbasis*nbasis),
     &     v(jend+2*nbasis*nbasis), lenv-2*nbasis*nbasis)

c$$$      call ascale (nbasis*(nbasis+1)/2, -1.0d0, v(jend), v(jend))
c$$$      call square (v(jend), v(jend), nbasis, nbasis, -1)

      call square (v(jend), v(jend), nbasis, nbasis, 0)
      
      call outoctfmt(iout, 'fock.mat', 1, 'fock', 2, dims, v(jend))         
      
      call chainx(0)
      
      end


