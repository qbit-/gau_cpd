      subroutine cpd2int (v, ngot)

      implicit none


C +------------------------------------------------------------+
C |  CPD-2INT   ---   the Scuseria group, 02.2015              |
C |                                                            |
C |                                                   v0.2 -   |
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
C |  iop(33) - printing level                                  |
C |     = 0, [default] print basic information                 |
C |         from Quad and PRISM                                |
C |     = 10, print two-electron integrals                     |
C |                                                            |
C |  iop(59) - integrals decomposition type                    |
C |       = 1 decompose AOs  (default)                         |
C |       = 2 decompose MOs                                    |
C |  iop(60) - number of radial and angular points in          |
C |            spatial quadrature:                             |
C |            IIIJJJ .. III radial points, JJJ angular points |
C |            from appropriate Lebedev grid. The default is   |
C |            10 radial and 26 angular points (7th order grid)|
C |                                                            |
C |  iop(61) - selects the radial weights?                     |
C |  iop(62) - selects the radial grid type?                   |     
C |  iop(63) - selects the tolerance of the pseudoinversion    |
C |                of spatial overlap. Default: 10 (10^-10)    |     
C |                                                            |
C |  iop(109)-  RI options:                                    |
C |           0 ... Default (2).                               |
C |           1 ... If there is a density basis file already   |
C |                 set up, then do RI.                        |
C |           2 ... Ignore any density basis file do not do RI.|
C |        NN00 ... Set the tolerance for discarding singular  |
C |                 values in SVD in RI to 10^-NN              |
C |                 (default 10^-6).                           | 
C |       10000 ... Use the Coulomb operator (default).        |
C |       20000 ... Use the overlap operator.                  |
C |                                                            |
C |  iop(77) - debug printing to Octave files                  |
C |                0 - do not print                            |
C |               >0 - print all matrices (X,S,E,Z)            |
C |  iop(78) - debugging of decomposition accuracy (N^4 cost)  |
C |                0 - do not perform                          |
C |               >0 - compare decomposed integrals to full    |
C |              >10 - print summary to a CSV file             |      
C |  iop(79) - requests to read external grid from grid.dat    |
C |                0 - build grid normally                     |
C |               >0 - read external grid                      |     
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

      logical     extgrid

      integer     dectype
      logical     dotran
      
C     memory allocation variables

      integer     jend, lenv
      
C     other variables

      integer     nbsq
      integer     ngrid, nbas6d
      integer     nrank, itol, ritol   
      integer     ii, junk
      real*8      accdes

C     transformation related variables

      integer     jxa, jxb, jscrx, jscrmap
      integer     lenx, kk
      integer     iopclx, methlc
      logical     units
      
      integer     nspblk, nri
      logical     cbfn, spinor
      
C     CPD2E RWF files

      integer     irwfac, irwspo, irwz, irwe

      parameter   ( irwfac = 824, irwspo = 826 )
      parameter   ( irwe = 827, irwz = 828 )

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

C     debug printing

      logical     dbgprint
      integer     dims(4)

C     debug integral decomposition
      
      logical     tests, testsprint

      integer     irw2efull, irwrec
      integer     irwria, irwrib
      integer     ndbf6d
      real*8      normfro, froful, frodiff
      real*8      outbuf(6)
      
C     debug Fock matrices
      
      integer     jrhoa, jrhob, jfa, jfb, jscr, ntt
      integer     nmat, nmats, nmatt
      logical     initrand
      character   mtname*10
      real*8      tr1a, tr1b, tr2a, tr2b, xx, time
      
C     functions for debug fock
      real*8      scftrc, omp_get_wtime
      
C     .. prism variables ..

#include "repall.inc"      
      integer     ipflag, nsymop, nop1
      integer     lseall, momega, nomega
      parameter   ( momega = DEFMOMEGA, nomega = DEFNOMEGA )
      real*8      omega(momega,nomega)
      integer     jneq, jneqshl, jneqsh2
      integer     lenneq, lenneqshl, lenneqsh2
      logical     allowp(50), fmm
      integer     jpflsta
      integer     ifmm, iseall
      integer     fmflag, fmflg1, nfxflg, ihmeth
      integer     itqry

      logical     keepp, initf, dopurp, dopurf, dospar
C       .. PBC variables ..

      integer     npdir, jpbcsta, jpbcend
      real*8      celvec(3,3)


      
C     format statements

 1001 format (' CPD2E was invoked with the following options:')
 1010 format (' iprint = ', I2)
 1011 format ('   grid selection:')      
 1012 format (' iradan = ', I6, ' iranwt = ', I5, ' irangd = ', I5)
 1013 format (' custom grid from file "grid.dat" will be used')
 1020 format (' accuracy = ',D8.3, ' pseudoinv. tolerance = ', D8.3)
 1021 format (' exact integrals will be used ')
 1022 format (' density fitted integrals will be used ')
 1023 format (' RI metric pseudoinv. tolerance = ', D8.3)
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

      iprint = iop(33)

      if ( iprint .lt. 0 .or. iprint .gt. 10 ) iprint = 0
      write (iout, 1010) iprint

C     Set up job type

      if (iop(59) .le. 0) then
         dectype = 1    ! AOs
      elseif (iop(59) .gt. 2 ) then
         call gauerr ('Incorrect decomposition type in cpd2e')
      else
         dectype = iop(59)
      endif
      
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

C     Decide on external grid
      
      if (iop(79) .gt. 0) then
         extgrid = .true.
         write (iout, 1013)
      endif
      
C     Desired accuracy and pseudoinversion tolerance for spatial overlap

      accdes = decacc (0, 0.0d0, .false.)
      itol  = mod(iop(63),100)
      if (itol .le. 0) then
         itol  = 10
      endif
      write (iout, 1020) accdes, Gfloat(10)**(-itol)
      
C     Use of symmetry
C     We currently set jsym2e = 0 in the code (users cannot change that).

      jsym2e = 0

C     debug printing options

      if ( iop(77) .gt. 0 ) then
         dbgprint = .true.
      else
         dbgprint = .false.
      endif

C     testing decomposition options

      if ( iop(78) .gt. 0 ) then
         tests = .true.
         write(iout,'(A)') ' Testing of the decomposition is requested '
         if ( ( iop(78) / 10 ) .gt. 0 ) then
            testsprint = .true.
            write(iout,'(A)') ' will output results to CSV file'
         else
            testsprint = .false.
         endif
      else
         tests = .false.
      endif
      
C     Use of RI (N^4 scaling algorithm) and related options

      iuseri = mod(iop(109), 10)
      if(iuseri .eq. 1) then
         ri = .true.
         write (iout, 1022)
      else
         ri = .false.
      endif

C     RI inversion tolerance
      if (ri) then
         ritol = mod(iop(109), 10000) / 100
         if (ritol .le. 0) then
            ritol = 10
         endif
         write (iout, 1023) Gfloat(10)**(-ritol)
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
     &  call gauerr ('cpd2int: Spherical bases are not supported yet')
      write (iout, *)

C-----------------------------------------------------------------------------+
C     Prepare transformation used with the quadrature                         |
C-----------------------------------------------------------------------------+

      if (dectype .eq. 1) then
         dotran = .false.
      elseif (dectype .eq. 2) then
         dotran = .true.
      elseif (dectype .eq. 3) then
         dotran = .false.
      endif
      
c$$$      if (dotran) then
c$$$         call decosc(iopcl,nspblk,nri,junk,cbfn,junk,spinor)
c$$$
c$$$         if (cbfn) then
c$$$            call gauerr('cpd2int: Complex bases are not supported yet')
c$$$         endif
c$$$
c$$$         lenx = itqry(irwx)
c$$$         if (lenx .ne. nbasis*nbsuse) then
c$$$            write (iout, '(A,I4,A,I4)') ' lenx expected = ',
c$$$     &           nbasis*nbsuse, ' but got = ', lenx
c$$$            write (iout, '(A,A)') ' Make sure to run l325 before l401',
c$$$     &       ' aborting ...'
c$$$!            call gauerr ('cpd2int: localization')
c$$$         endif
c$$$
c$$$         jxa     = jend
c$$$         jxb     = jxa + nbasis*nbasis + lenx
c$$$         jscrmap = jxa + nbasis*nbasis
c$$$         jend    = jscrmap + nbasis 
c$$$         
c$$$         call tstcor(jend ,lenv,'cpd2int: localization')
c$$$
c$$$         call fileio (2, -irwx, lenx, v(jxa), 0)
c$$$         dims(1) = nbasis
c$$$         dims(2) = nbsuse
c$$$         call outoctfmt(iout, 't.mat', 1, 't', 2, dims, v(jxa))
c$$$
c$$$C     Modify iopcl as localize does not like complex
c$$$
c$$$         iopclx = 0
c$$$         units  = .false.
c$$$         methlc = 2
c$$$
c$$$        call localize(iout,5,iopclx,methlc,units,natoms,
c$$$     &        nbsuse,nbasis,v(jxa),xx,v(jscrmap),v(jend),lenv-jend)
c$$$        dims(1) = nbasis
c$$$        dims(2) = nbsuse
c$$$        call outoctfmt(iout, 'j.mat', 1, 'j', 2, dims, v(jend))
c$$$         
c$$$      endif
      
C-----------------------------------------------------------------------------+
C     Compute the collocation matrix X(mu,G) and save it to the RWF irwfac    |
C-----------------------------------------------------------------------------+

      if (extgrid) then
         call frmextquad (iout, iprint, 'grid.dat', 0, 
     &        nbas6d, ngrid, irwfac, v(jend), lenv - jend)
      else
         call frmquad (iout, iprint, iradan, iranwt, irangd, iopcl,
     &        dotran, nbsuse, v(jxa), ngrid, irwfac, v(jend),
     &        lenv - jend)
      endif

C-----------------------------------------------------------------------------+
C     Scale quadrature to coincide with PRISM scaling, possibly convert to    |
C     5d/7f                                                                   |
C-----------------------------------------------------------------------------+
      call convquad (iout, iprint, .true.,
     &     .false., .false., ngrid, irwfac, v(jend), lenv - jend)
      
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
         call frmemtri(iout, iprint, nbas6d, ngrid, jsym2e, neq, neqshl,
     &        neqsh2, ir12, ritol, accdes, ndbf6d, irwfac,
     &        irwria, irwrib, irwe, v(jend), lenv - jend)
      endif
      
C-----------------------------------------------------------------------------+
C     Contract projected integrals with the spatial overlap                   |
C-----------------------------------------------------------------------------+

      call frmz(iout, iprint, ngrid, irwspo, irwe, irwz, v(jend), lenv
     & - jend) 

C-----------------------------------------------------------------------------+
C     Print some summary                                                      |
C-----------------------------------------------------------------------------+
      if (iprint .gt. 0) then
         write (iout, '(A,I20)') ' basis size = ', nbasis
         write (iout, '(A,I20)') ' grid size = ', ngrid
         write (iout, '(A,I20)') ' overlap rank = ', nrank
         if (ri) write (iout, '(A,I20)') ' ndbf = ', ndbf6d
      endif

      
C     save number of grid points for other routines

      call ilsw (1, 80, ngrid)
      
C----------------------------------------------------------------------------+
C     debug printing                                                         |
C----------------------------------------------------------------------------+
      
      if (dbgprint) then
         call tstcor(max(ngrid*ngrid,nbasis*ngrid), lenv - jend,
     &        'cpd2int: debug print')
         
         call fileio(2,-irwfac,nbasis*ngrid,v(jend),0)
         if (.not.dotran) then
            dims(1) = nbasis
         else
            dims(1) = nbsuse
         endif
         dims(2) = ngrid
         call outoctfmt(iout, 'x.mat', 1, 1,'x', 2, dims, v(jend))

         call fileio(2,-irwspo,ngrid*(ngrid+1)/2,v(jend),0)
         call square (v(jend), v(jend), ngrid, ngrid, 0)
         dims(1) = ngrid
         dims(2) = ngrid
         call outoctfmt(iout, 's.mat', 1, 1,'s', 2, dims, v(jend))
         
         call fileio(2,-irwe,ngrid*ngrid,v(jend),0)
         call outoctfmt(iout, 'e.mat', 1, 1,'e', 2, dims, v(jend))
         
         call fileio(2,-irwz,ngrid*ngrid,v(jend),0)
         call outoctfmt(iout, 'z.mat', 1, 1,'z', 2, dims, v(jend))         

         call tstcor (nbasis**4, lenv - jend, 'cpd2int: testing')

         call frm2eints (iout, iprint, jsym2e, neq, neqshl, neqsh2,
     &        accdes, .false.,.true., junk, irw2efull, v(jend), lenv
     &    - jend)

         call fileio(2,-irw2efull,nbasis**4,v(jend),0)
         dims(1) = nbasis
         dims(2) = nbasis
         dims(3) = nbasis
         dims(4) = nbasis
         call outoctfmt(iout, 'ints.mat', 1, 1,'ints', 4, dims, v(jend))         
      endif

C----------------------------------------------------------------------------+
C     testing the decomposition ( O(N^4) cost )                              |
C----------------------------------------------------------------------------+
      
      if (tests) then
         call frm2eints (iout, iprint, jsym2e, neq, neqshl, neqsh2,
     &        accdes, .false.,.true., junk, irw2efull, v(jend),
     &   lenv - jend)

c$$$         call tstcor (nbasis**4, lenv - jend, 'cpd2int: testing')
c$$$         call fileio(2,-irw2efull,nbasis**4,v(jend),0)
c$$$         dims(1) = nbasis
c$$$         dims(2) = nbasis
c$$$         dims(3) = nbasis
c$$$         dims(4) = nbasis
c$$$         call outoctfmt(iout, 'ints.mat', 1, 1,'ints', 4, dims, v(jend))         
         
         froful=normfro(iout, iprint, nbasis, irw2efull, v(jend), lenv
     &    - jend)
         
         call reconstr (iout, iprint, nbasis, ngrid,
     &        irwfac, irwz, irwrec, v(jend), lenv - jend)
         
         call asubf(iout,iprint,nbasis,irw2efull,irwrec, v(jend), lenv
     &    - jend) 
         frodiff=normfro (iout, iprint, nbasis, irwrec, v(jend), lenv
     &    - jend)

         write (iout, '(A,I20)')    ' nbasis = ', nbasis
         write (iout, '(A,I20)')    ' nbsuse = ', nbsuse
         write (iout, '(A,I20)')    ' ngrid  = ', ngrid
         write (iout, '(A,I20)')    ' rank   = ', nrank
         if (ri) write (iout, '(A,I20)') ' ndbf = ', ndbf6d
         write (iout, '(A,E20.12)') ' full   = ', froful
         write (iout, '(A,E20.12)') ' diff   = ', frodiff
         write (iout, '(A,E20.12)') ' frac % = ', frodiff / froful * 100
         
         if (testsprint) then
            outbuf(1) = dble(nbasis)
            outbuf(2) = dble(ngrid)
            outbuf(3) = dble(nrank)
            outbuf(4) = froful
            outbuf(5) = frodiff
            outbuf(6) = frodiff / froful
            call outcsv (iout, 'l325_summary.csv', 2, 6, outbuf)
         endif
         
      endif

C---------------------------------------------------------------------------+
C     testing Fock matrices                                                 |
C---------------------------------------------------------------------------+      
      
      dbgprint = .false.
      if (dbgprint) then
         ntt   = nbasis*(nbasis+1)/2
         nbsq  = nbasis*nbasis

         nmat  = 6
         nmats = 2
         nmatt = 2
         iopcl = 1

         jscr  = jend
         jrhoa = jscr  + nbsq 
         jrhob = jrhoa + nmat*ntt
         jfa   = jrhob + nmat*ntt
         jfb   = jfa   + nmat*ntt
         jend  = jfb   + nmat*ntt
         
         dims(1) = nbasis
         dims(2) = nbasis

         initrand = .true.
         time = omp_get_wtime()
         call srand(time)
  
         do 1100 ii = 0, nmats - 1
            call symmetric (iout, .false., nbasis, 0, v(jscr))
            write(mtname, '(A,I3.3)') 'pas', ii 
            call outoctfmt (iout, 'pa.mat', 2, 1, mtname, 2, dims,
     &           v(jscr))
            call linear (v(jscr), v(jrhoa + ii*ntt), nbasis, nbasis)
 1100    continue
         do 1200 ii = 0, nmatt - 1 
            call symmetric (iout, .false., nbasis, 0, v(jscr))
            write(mtname, '(A,I3.3)') 'pat', ii 
            call outoctfmt (iout, 'pa.mat', 2, 1, mtname, 2, dims,
     &           v(jscr))
            call linear (v(jscr), v(jrhoa+ntt*nmats+ii*ntt), nbasis,
     &    nbasis)
 1200    continue
         do 1300 ii = 0, nmat - nmats - nmatt - 1 
            call symmetric (iout, .false., nbasis, 1, v(jscr))
            write(mtname, '(A,I3.3)') 'paa', ii             
            call outoctfmt (iout, 'pa.mat', 2, 1, mtname, 2, dims,
     &           v(jscr))
            call linear (v(jscr), v(jrhoa+(nmats+nmatt)*ntt+ii*ntt),
     &           nbasis, nbasis)
 1300    continue

         do 1150 ii = 0, nmats - 1
            call symmetric (iout, .false., nbasis, 0, v(jscr))
            write(mtname, '(A,I3.3)') 'pbs', ii 
            call outoctfmt (iout, 'pb.mat', 2, 1, mtname, 2, dims,
     &           v(jscr))
            call linear (v(jscr), v(jrhob + ii*ntt), nbasis, nbasis)
 1150    continue
         do 1250 ii = 0, nmatt - 1 
            call symmetric (iout, .false., nbasis, 0, v(jscr))
            write(mtname, '(A,I3.3)') 'pbt', ii             
            call outoctfmt (iout, 'pb.mat', 2, 1, mtname, 2, dims,
     &           v(jscr))
            call linear (v(jscr), v(jrhob+ntt*nmats+ii*ntt), nbasis,
     &    nbasis)
 1250    continue
         do 1350 ii = 0, nmat - nmats - nmatt - 1 
            call symmetric (iout, .false., nbasis, 1, v(jscr))
            write(mtname, '(A,I3.3)') 'pba', ii             
            call outoctfmt (iout, 'pb.mat', 2, 1, mtname, 2, dims,
     &           v(jscr))
            call linear (v(jscr), v(jrhob+(nmats+nmatt)*ntt+ii*ntt),
     &           nbasis, nbasis)
 1350    continue

C         set up PBC information
      
         call getpbc(npdir, celvec)
         jpbcsta = jend
         jpbcend = jpbcsta
         call rdipbc(1, npdir, 0, jpbcsta, jpbcend, v, lenv)
         jend = jpbcend
         
C     set up PRISM control flags

         ifmm = 0
         iseall = -1
         
         jpflsta = jend

         call setpfl (iout, iprint, ifmm, ipflag, allowp, fmm, fmflag,
     $        fmflg1, nfxflg, ihmeth, omega, iseall,lseall,jpflsta,
     $        v(jend),
     $        lenv-jend+1)

C     load symmetry information

        lenneq    = intowp (1*nbasis)
        lenneqshl = max (itqry (), 0)
        lenneqsh2 = max (itqry (), 0)

        jneq    = jend
        jneqshl = jneq     + lenneq
        jneqsh2 = jneqshl  + lenneqshl
        jend    = jneqsh2  + lenneqsh2

C     Fill arrays with symmetry related quantities.

        if ( lenneq .ne. 0 ) then
           call fileio (2, -580, lenneq, v(jneq), 0)
        endif

        if ( lenneqshl .ne. 0 ) then
           call fileio (2, -565, lenneqshl, v(jneqshl), 0)
        endif

        if ( lenneqsh2 .ne. 0 ) then
           call fileio (2, -726, lenneqsh2, v(jneqsh2), 0)
        endif

        keepp  = .true.
        initf  = .true.
        dopurp = .false.
        dopurf = .false.
        dospar = .false.
        
         call fofcou (iout,iprint,1,iopcl,0,1,nomega,
     $ omega,fmm,fmflag,fmflg1,nfxflg,ipflag,allowp,lseall,
     $ v(iseall),
     $ 1,natoms,keepp,initf,dopurp,dopurf,dospar,0,accdes,
     $ (1.0d0,0.0d0),0,nmat,nmats,nmatt,nmat,0,0,nbasis,
     $ 0,1,1,v(jneq),v(jneqshl),v(jneqsh2),rotall,neqall,0,0,v(jrhoa),
     $ v(jrhob),junk,junk,junk,junk,junk,junk,junk,0,v(jfa),v(jfb),0,
     $ junk,junk,natoms,ian,atmchg,c,iattyp,iatfrg,xx,xx,
     $ v(jpbcsta),celvec,(0.0d0,0.0d0),(1.0d0,0.0d0),xx,0,0,xx,xx,junk,
     $ junk,xx,v(jend),v(jend),lenv-jend+1)

      do 2100 ii = 0, nmats - 1
            write(mtname, '(A,I3.3)') 'fas', ii 
            call outoctfmt(iout, 'fa.mat', 2, 1, mtname, 1, ntt,
     &           v(jfa+ii*ntt))         
 2100    continue
         do 2200 ii = 0, nmatt - 1
            write(mtname, '(A,I3.3)') 'fat', ii 
            call outoctfmt(iout, 'fa.mat', 2, 1, mtname, 1, ntt,
     &           v(jfa+nmats*ntt+ii*ntt))         
 2200    continue
         do 2300 ii = 0, nmat - nmats - nmatt - 1
            write(mtname, '(A,I3.3)') 'faa', ii 
            call outoctfmt(iout, 'fa.mat', 2, 1, mtname, 1, ntt,
     &           v(jfa+(nmats+nmatt)*ntt+ii*ntt))         
 2300    continue

         do 2150 ii = 0, nmats - 1
            write(mtname, '(A,I3.3)') 'fbs', ii             
            call outoctfmt(iout, 'fb.mat', 2, 1, mtname, 1, ntt,
     &           v(jfb+ii*ntt))         
 2150    continue
         do 2250 ii = 0, nmatt - 1
            write(mtname, '(A,I3.3)') 'fbt', ii             
            call outoctfmt(iout, 'fb.mat', 2, 1, mtname, 1, ntt,
     &           v(jfb+nmats*ntt+ii*ntt))         
 2250    continue
         do 2350 ii = 0, nmat - nmats - nmatt - 1
            write(mtname, '(A,I3.3)') 'fba', ii             
            call outoctfmt(iout, 'fb.mat', 2, 1, mtname, 1, ntt,
     &           v(jfb+(nmats+nmatt)*ntt+ii*ntt))         
 2350    continue

         tr1a = scftrc(v(jrhoa), v(jfa), nbasis, nmat)
         tr1b = scftrc(v(jrhob), v(jfb), nbasis, nmat)
         
         call cpdfock (iout, iprint, iopcl, -1, nbasis, ngrid, 
     &    irwfac,
     &        irwz, nmat, nmats, nmatt, v(jrhoa), v(jrhob), v(jfa),
     &        v(jfb),
     &        v(jend), lenv-jend)

         do 3100 ii = 0, nmats - 1
            write(mtname, '(A,I3.3)') 'facs', ii             
            call outoctfmt(iout, 'fac.mat', 2, 1, mtname, 1, ntt,
     &           v(jfa+ii*ntt))         
 3100    continue
         do 3200 ii = 0, nmatt - 1
            write(mtname, '(A,I3.3)') 'fact', ii             
            call outoctfmt(iout, 'fac.mat', 2, 1, mtname, 1, ntt,
     &   v(jfa+nmats*ntt+ii*ntt))         
 3200    continue
         do 3300 ii = 0, nmat - nmats - nmatt - 1
            write(mtname, '(A,I3.3)') 'faca', ii 
            call outoctfmt(iout, 'fac.mat', 2, 1, mtname, 1, ntt,
     &           v(jfa+(nmats+nmatt)*ntt+ii*ntt))         
 3300    continue

         do 3150 ii = 0, nmats - 1
            write(mtname, '(A,I3.3)') 'fbcs', ii 
            call outoctfmt(iout, 'fbc.mat', 2, 1, mtname, 1, ntt,
     &           v(jfb+ii*ntt))         
 3150    continue
         do 3250 ii = 0, nmatt - 1
            write(mtname, '(A,I3.3)') 'fbct', ii             
            call outoctfmt(iout, 'fbc.mat', 2, 1, mtname, 1, ntt,
     &   v(jfb+nmats*ntt+ii*ntt))         
 3250    continue
         do 3350 ii = 0, nmat - nmats - nmatt - 1
            write(mtname, '(A,I3.3)') 'fbca', ii 
            call outoctfmt(iout, 'fbc.mat', 2, 1, mtname, 1, ntt,
     &           v(jfb+(nmats+nmatt)*ntt+ii*ntt))         
 3350    continue

         tr2a = scftrc(v(jrhoa), v(jfa), nbasis, nmat)
         tr2b = scftrc(v(jrhob), v(jfb), nbasis, nmat)

         write(iout, *)
         write(iout, '(A,E20.12,A,E20.12,A,E20.12)')
     &        ' full = ', tr1a, ' diff = ', (tr1a - tr2a), ' rel = ',
     &        (tr1a - tr2a) / tr1a
         write(iout, '(A,E20.12,A,E20.12,A,E20.12)')
     &        ' full = ', tr1b, ' diff = ', (tr1b - tr2b), ' rel = ',
     &        (tr1b - tr2b) / tr1b
      endif
      
      call chainx(0)
      
      end


