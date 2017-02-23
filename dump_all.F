      subroutine dump_all (v, ngot)

      implicit none


C +------------------------------------------------------------+
C |  Dump all   ---   dumps quadratures, integrals and         |
C |                   everything needed for an external        |
C |                   program to perform CPD                   |
C |                                                            |
C |                                                            |
C |    ( link 321, part of the Gaussian system of programs )   |
C |                                                            |
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
C |  iop(79) - requests to read external grid from grid.dat    |
C |                0 - build grid normally                     |
C |               >0 - read external grid                      |
C |  iop(80) - what to output                                  |
C |                0 - energies only (Default)                 |
C |                1 - dump quadrature X                       |
C |                2 - dump overlap and AO to orth AO transform|
C |                4 - dump MOs and eigenvalues                |
C |               10 - dump DF integrals and DF overlap        |
C |               20 - dump full integrals                     |
C |               40 - dump T2 amplitudes                      |
C |              100 - dump localized MOs                      |
C |         options can be added, example: 13 - DF integrals + |
C |                  quadrature + overlap                      |
C |                                                            |
C |  iop(81) - options to localize                             |
C |                0 - Boys                                    |
C |                1 - Population                              |
C |                2 - Boys + Population                       |
C |               00 - localize occ/virtuals separately        |
C |               10 - localize occ/virtuals together          |
C |                                                            |
C |  iop(99) - options for the output file                     |
C |                0 - Rewrite any previous file (default)     |
C |                1 - Append. Warning: this does not check    |
C |                    what was inside. If same variable is    |
C |                    dumped twice it is your duty to check   |
C |                    what a user program will do.            |
C |              XX0 - suffix of the dump_allXX.mat file       |
C |  iop(109)-  RI options:                                    |
C |           0 ... Default do not use RI                      |
C |           1 ... Force RI use                               |
C |           2 ... If there is a density basis file set up,   |
C |                 then use RI.                               |
C |        NN00 ... Set the tolerance for discarding singular  |
C |                 values in SVD in RI to 10^-NN              |
C |                 (default 10^-10).                          |
C |       10000 ... Use the Coulomb operator (default).        |
C |       20000 ... Use the overlap operator.                  |
C +------------------------------------------------------------+

C     include some Gaussian files

#include "commonmol.inc"
#include "commonb2.inc"

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

c     /orb/ common - various orbital dimensions

      common/orb/ispect,lspect,nrorb,noa,nva,nob,nvb,noaob,noava,noavb,
     &  nobva,nobvb,nvavb,noa2,noa3,nob2,nob3,nva2,nva3,nvb2,nvb3,
     &  novaa,novab,novbb,maxbuc,ieval,ioab,loab

      integer      ispect,lspect,nrorb,noa,nva,nob,nvb,noaob,noava,noavb
      integer      nobva,nobvb,nvavb,noa2,noa3,nob2,nob3,nva2,nva3,nvb2
      integer      nvb3,novaa,novab,novbb,maxbuc,ieval,ioab,loab
      integer      lenorb

      integer      nvava, nvbvb, novaat, novbbt, novabt

c     window information variables
      integer      iwindow(2)
      integer      lenwind
      real*8       acta(3)

c     general array

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
      logical     dopure

      integer     idump, idumpdef, idp
      logical     outquad, outaox, outmos, outint, outintdf, outamp
      integer     ilocal, methloc
      logical     outlocal

      parameter   (idumpdef = 0)

C     memory allocation variables

      integer     jend, jovoff, jscr, lenx, next, lenv, lenscr
      integer     jatmap

C     other variables

      integer     nbsq, nbsq6d
      integer     ngrid, nbas6d
      integer     ndbf
      integer     nrank, itol, ritol
      integer     ii, junk
      real*8      accdes
      real*8      xx
      data        junk/0/, xx/0.0d0/

      integer     ndims
      integer     dims(4)

      integer     ipured, ipuref


C     iopcl deciphering

      integer     nspblk, nri
      logical     cbfn, spinor

C     file handling

      integer    ios
      character*(*)   filenmdef
      parameter   ( filenmdef = 'dump_all')

      integer        fnmsuff
      character*200 filenm

      logical    rewrite

C     Needed for commonb2 to work with implicit none

      integer    nshelb, janb, maxtyb

C     RWF files

      integer     irwfac, irwspo, irwz, irwe, irwamp
      integer     irworb, irwwd

      parameter   ( irwfac = 824, irwspo = 826 )
      parameter   ( irwe = 827, irwz = 828 )
      parameter   ( irworb = 545, irwamp = 570 )
      parameter   ( irwwd = 569 )

C     Scratch RWFs

      integer     irw2efull
      integer     irwria, irwrib

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
      integer     intowp, itqry

C     format statements

 1001 format (' Dump all was invoked with the following options:')
 1010 format (' iprint = ', I2)
 1011 format ('   grid selection:')
 1012 format (' iradan = ', I6, ' iranwt = ', I5, ' irangd = ', I5)
 1013 format (' custom grid from file "grid.dat" will be used')
 1020 format (' accuracy = ',D8.3, ' pseudoinv. tolerance = ', D8.3)
 1021 format (' exact integrals will be computed ')
 1022 format (' density fitted integrals will be computed ')
 1023 format (' RI metric pseudoinv. tolerance = ', D8.3)
 1040 format (' integrals will be replicated using symmetry')
 1050 format (' from /B/ 5d/6d is: ',I2,' 7f/10f is: ',I2)
 1060 format (' from ILSW 5d/6d is: ',I2,' 7d/10f is:',I2)
 1070 format (' Orbilal sizes: NOA=',I3,' NOB=',I3)
 1080 format (' Orbital window: ',I3,' - ',I3)

      write (iout, *) ' '
      write (iout, *) ' +-------------------------------------------+'
      write (iout, *) ' |                                           |'
      write (iout, *) ' |  LINK 330 Dumping information for CPD     |'
      write (iout, *) ' |                                           |'
      write (iout, *) ' +-------------------------------------------+'
      write (iout, *) ' '

      call drum (v, ngot)
      jend = 1
      lenv = ngot

C----------------------------------------------------------------------+
C       Process options for Dump all                                   |
C----------------------------------------------------------------------+

      write (iout, *)
      write (iout, 1001)
      write (iout, *)

C     read /gen/

      call fileio (2, -irwgen, lengen, dgen, 0)

C     Set printing level.

      iprint = iop(33)

      if ( iprint .lt. 0 .or. iprint .gt. 10 ) iprint = 0
      write (iout, 1010) iprint

C     decide on the file handling

      if ( mod(iop(99), 10) .eq. 0 ) then
         rewrite = .true.
      else
         rewrite = .false.
      endif

      fnmsuff = mod(iop(99), 1000) / 10

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

C     Try to recover iradan from a previous run

      if (iradan .eq. 0) then
         call ilsw (2, 80, iradan)
      endif

C     Defaults if still nothing is set
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

C     Decide on what to output

      if (iop(80) .le. 0) then
         idump = idumpdef
      else
         idump = iop(80)
      endif

      idp = mod(idump,10)
      outmos =  ( (idp / 4) .ge. 1 )
      idp = mod(idp,4)
      outaox = ( (idp / 2) .eq. 1 )
      idp = mod(idp,2)
      outquad = (idp .eq. 1)

      idp = ( (mod(idump,100) / 10) )
      outamp = ( ( idp / 4) .ge. 1 )
      idp = mod(idp,4)
      outint = ( (idp / 2) .ge. 1 )
      idp = mod(idp,2)
      outintdf = (idp .eq. 1)

      idp = ( (mod(idump,1000) / 100) )
      outlocal = (idp .eq. 1)

C     Options to MO localization routines

      methloc = (mod(iop(81),10))
      ilocal  = (mod(iop(81),100) / 10)

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

C     Use of RI and related options

      iuseri = mod(iop(109), 10)
      if(iuseri .eq. 1) then
         ri = .true.
         write (iout, 1022)
      elseif (iuseri .eq. 2) then
         if (itqry(724) .gt. 0) then
            ri = .true.
         else
            ri = .false.
         endif
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

      call rwmol(2,irwmol)

C     Get 6d basis used with all quadratures

      call getnb6(nbas6d)
C      if (nbasis .ne. nbas6d)
C     &  call gauerr ('cpd2int: Spherical bases are not supported yet')
      write (iout, *)
      nbsq = nbasis*nbasis
      nbsq6d = nbas6d*nbas6d

      if (ri) then
C     Load /B2/ with the density basis.
         call rdwrb(-10,0)
         call getmbn(I5DB2,I7FB2,NShelB,ShelTB,ShelCB,ndbf)
      endif

C     figure out purity

      call rdwrb(2,0)

      write(iout, 1050) I5DB1, I7FB1

C     Get 6d basis used with all quadratures

      call getnb6(nbas6d)
      if (nbasis .ne. nbas6d)
C     &  call gauerr('Spherical basis is different from cart')
     & write(iout,*) 'Spherical basis is different from cart'
      write (iout, *)

C-----------------------------------------------------------------------------+
C     File initialization                                                     |
C-----------------------------------------------------------------------------+

      if (fnmsuff .ne. 0) then 
         write(filenm, '(A,i2,A)') filenmdef, fnmsuff, '.mat'
      else
         write(filenm, '(A,A)') filenmdef, '.mat'
      endif
      
      if ( rewrite ) then
C     Clean the previous file if it existed
         open(unit=771, iostat=ios, file=filenm, status='old')
         if (ios .eq. 0) then
            close(771, status='delete')
         endif
      endif

      
C-----------------------------------------------------------------------------+
C     Dump some scalars                                                       |
C-----------------------------------------------------------------------------+

      call fileio (2, -irwgen, lengen, dgen, 0)

C     SCF energy
      if ( rewrite ) then
         call outoctfmt (iout, filenm, 1, 1, 'e_scf', 1, 1, dgen(32))
      else
         call outoctfmt (iout, filenm, 2, 1, 'e_scf', 1, 1, dgen(32))
      endif

C     CC/CI energy
      call outoctfmt (iout, filenm, 2, 1, 'e_cc', 1, 1, dgen(9))
C     MP2 energy
      call outoctfmt (iout, filenm, 2, 1, 'e_mp2', 1, 1, dgen(33))


C-----------------------------------------------------------------------------+
C-----------------------------------------------------------------------------+

C     Dump quadrature

      if (outquad) then
         if (extgrid) then
C     Load X(mu,G) from grid.dat and save it to the RWF irwfac
            call frmextquad (iout, iprint, 'grid.dat', 0,
     &           nbas6d, ngrid, irwfac, v(jend), lenv - jend)
         else
C     Or compute X(mu,G) and save it to the RWF irwfac
            call frmquad (iout, 2, iradan, iranwt, irangd, iopcl,
     &           .false., nbsuse, junk, ngrid, irwfac, v(jend),
     &           lenv - jend)
         endif

         call tstcor(nbas6d*ngrid,lenv,'dump_all: quadrature')
         call fileio (2, -irwfac, nbas6d*ngrid, v(jend), 0)

         dims(1) = nbas6d
         dims(2) = ngrid
c        This is a raw 6d/10f matrix
         call outoctfmt (iout, filenm, 2, 1,'Xu', 2, dims,
     &        v(jend))

         call fileio (2, -irwfac, nbas6d*ngrid, v(jend), 0)

c        This prepares a pure collocation matrix
         if (nbasis .ne. nbas6d) then
            call convquad (iout, iprint, .false.,
     &           .true., .false., ngrid, irwfac, v(jend), lenv - jend)
            call fileio (2, -irwfac, nbasis*ngrid, v(jend), 0)
            dims(1) = nbasis
            dims(2) = ngrid
            call outoctfmt (iout, filenm, 2, 1,'X', 2, dims,
     &     v(jend))
         else
            dims(1) = nbas6d
            dims(2) = ngrid
            call outoctfmt (iout, filenm, 2, 1,'X', 2, dims,
     &           v(jend))
         endif

      endif

C-----------------------------------------------------------------------------+
C     Compute integrals                                                       |
C-----------------------------------------------------------------------------+

      if (outint) then
         call frm2eints (iout, iprint, jsym2e, neq, neqshl, neqsh2,
     &        accdes, .true.,.false., next, irw2efull, v(jend),
     &        lenv - jend)

         dims(1) = nbasis
         dims(2) = nbasis
         dims(3) = nbasis
         dims(4) = nbasis
         call outoctfmt (iout, filenm, 2, 1, 'I', 4, dims, v(jend))
      endif

      if (outintdf) then
         if (.not. ri)
     & call gauerr ('dump_all: DF integrals requested, but RI not set')

         call  frm2eri (iout, iprint, jsym2e, neq, neqshl,
     &     neqsh2, ir12, .true., .true., ritol, accdes, jovoff,
     &     next, irwria, irwrib, v, lenv)

         dims(1) = ndbf
         dims(2) = ndbf
         call outoctfmt (iout, filenm, 2, 1, 'J', 2, dims,
     &        v(jend))


C         call genaos(I5DB1, I7FB1, NShell, ShellT, ShellC, junk, AOS)

c$$$         do 150 ii = 0,ndbf-1
c$$$            call adjdf3(.true.,.true.,MaxTyp,NShell,ShellT,ShellC,
c$$$     &           v(jend + jovoff + nbsq6d*ii))
c$$$ 150     continue
c$$$
c$$$         call purcar(iout,2,.true.,.true.,.false.,
c$$$     $    .false.,.false.,.false.,.true.,.false.,.false.,nbasis,nbas6d,
c$$$     $    1,ndbf,1,0,0,I5DB1,I7FB1,v(jend+jovoff),xx,
c$$$     $        xx,xx,xx,xx,xx,xx,xx,v(jend + jovoff+nbsq6d*ndbf),
c$$$     $        lenv - jend - nbsq6d*ndbf)

c$$$         Subroutine PurCar(IOut,IPrint,ToPure,DoX,DoZ,
c$$$     $        DoCISG,DoC,DoXD,SqMat,HaveXB,DoDens,NBasis,NBas6D,
c$$$     $        NMtPBC,NMat,NMatS,NMO,NMtFkD,IPureD,IPureF,
c$$$     $        XA,XB,PZA,PZB,CIDenA,CIDenB,CMO,XDA,XDB,V,MDV)

         dims(1) = nbasis
         dims(2) = nbasis
         dims(3) = ndbf
         call outoctfmt (iout, filenm, 2, 1, 'K', 3, dims,
     &        v(jend + jovoff))

      endif

C-----------------------------------------------------------------------------+
C     Dump AO overlap and other quantities                                    |
C-----------------------------------------------------------------------------+

      call decosc(iopcl,nspblk,nri,junk,cbfn,junk,spinor)

      if (cbfn) then
         call gauerr('dump_all: Complex bases are not supported yet')
      endif

      if (outaox) then
         call tstcor(nri*nbasis*nbsuse, lenv,'dump_all: T')
         call fileio (2, -irwx, nri*nbasis*nbsuse, v(jend), 0)
         dims(1) = nbasis
         dims(2) = nbsuse
         if (nri .eq. 1) then
            call outoctfmt (iout, filenm, 2, 1, 'T', 2, dims, v(jend))
         elseif (nri .eq. 2) then
            call outoctfmt (iout, filenm, 2, 2, 'T', 2, dims, v(jend))
         endif

         call tstcor(nri*nbasis*nbasis, lenv,'dump_all: S')
         call fileio (2, -irws, nri*nbasis*(nbasis+1)/2, v(jend), 0)
         if (nri .eq. 1) then
            call square (v(jend), v(jend), nbasis, nbasis, 0)
         elseif (nri .eq. 2) then
            call zsquare (v(jend), v(jend), nbasis, nbasis, 0)
         endif
         dims(1) = nbasis
         dims(2) = nbasis
         if (nri .eq. 1) then
            call outoctfmt (iout, filenm, 2, 1, 'S', 2, dims, v(jend))
         elseif (nri .eq. 2) then
            call outoctfmt (iout, filenm, 2, 2, 'S', 2, dims, v(jend))
         endif
      endif

C--------------------------------------------------------------------------
C     Dumping of MOs
C     Here 6d/10f MOs are adjusted to 5d/7f
C--------------------------------------------------------------------------

      if (outmos) then
C     First load and save window information
C     Load /orb/
         lenorb = intowp(28)
         call fileio(2,-irworb,lenorb,ispect,0)
         write(iout,1070) noa, nob

         lenwind = intowp(2)
         call fileio(2,-irwwd,lenwind,iwindow,0)
         write(iout, 1080) iwindow(1), iwindow(2)

         acta(1) = iwindow(1)
         acta(3) = iwindow(2)
         acta(2) = iwindow(1) + noa - 1

         call outoctfmt (iout, filenm, 2, 1, 'acta', 1, 3,
     &    acta)
         if (nspblk .eq. 2) then
            acta(2) = iwindow(1) + nob - 1
            call outoctfmt (iout, filenm, 2, 1, 'actb', 1, 3,
     &           acta)
         endif

C     Now load and dump MO coefficients
         call tstcor(nri*nbasis*nbasis, lenv,'dump_all: C')
         call aclear (nri*nbasis*nbasis, v(jend))
         call fileio (2, -irwca, nri*nbasis*nbsuse, v(jend), 0)
         call trspmo(.true.,.true.,nri,nbasis,nbasis,nbsuse,
     &        v(jend),v(jend + nri*nbasis*nbasis))
         call adjmo2(iout,iprint,.true.,I5DB1,I7FB1,nbasis,
     &        nbsuse*nri,NShell,ShellT,ShellC,v(jend))
         call trspmo(.false.,.true.,nri,nbasis,nbasis,
     &        nbsuse,v(jend),v(jend + nri*nbasis*nbasis))

         dims(1) = nbasis
         dims(2) = nbsuse
         if (nri .eq. 1) then
            call outoctfmt (iout, filenm, 2, 1, 'CA', 2, dims, v(jend))
         elseif (nri .eq. 2) then
            call outoctfmt (iout, filenm, 2, 2, 'CA', 2, dims, v(jend))
         endif

         if (nspblk .eq. 2) then
            call aclear (nri*nbasis*nbasis, v(jend))
            call fileio (2, -irwcb, nri*nbasis*nbsuse, v(jend), 0)
            call trspmo(.true.,.true.,nri,nbasis,nbasis,nbsuse,
     &           v(jend),v(jend + nri*nbasis*nbasis))
            call adjmo2(iout,iprint,.true.,I5DB1,I7FB1,nbasis,
     &           nbsuse*nri,NShell,ShellT,ShellC,v(jend))
            call trspmo(.false.,.true.,nri,nbasis,nbasis,
     &           nbsuse,v(jend),v(jend + nri*nbasis*nbasis))

            dims(1) = nbasis
            dims(2) = nbsuse
            if (nri .eq. 1) then
               call outoctfmt (iout, filenm, 2, 1, 'CB', 2, dims,
     &          v(jend))
            elseif (nri .eq. 2) then
               call outoctfmt (iout, filenm, 2, 2, 'CB', 2, dims,
     &          v(jend))
            endif
         endif

         call fileio (2, -irweig, 2*nbsuse, v(jend), 0)
         call outoctfmt (iout, filenm, 2, 1, 'ee', 1, 2*nbsuse, v(jend))
      endif

C-----------------------------------------------------------------------
C     Localization (does not handle proper 5d/7f conversion for now)
C-----------------------------------------------------------------------

      if (outlocal) then
         call tstcor(nri*nbasis*nbasis + nbasis*natoms,
     &        lenv,'dump_all: loc')
         call aclear (nri*nbasis*nbasis + nbasis*natoms, v(jend))
         jatmap = jend
         jend   = jatmap + nbasis*natoms

         call fileio (2, -irwca, nri*nbasis*nbsuse, v(jend), 0)
c$$$         call trspmo(.true.,.true.,nri,nbasis,nbasis,nbsuse,
c$$$     &        v(jend),v(jend + nri*nbasis*nbasis))
c$$$         call adjmo2(iout,iprint,.true.,I5DB1,I7FB1,nbasis,
c$$$     &        nbsuse*nri,NShell,ShellT,ShellC,v(jend))
c$$$         call trspmo(.false.,.true.,nri,nbasis,nbasis,
c$$$     &        nbsuse,v(jend),v(jend + nri*nbasis*nbasis))

         call localize(iout, iprint, 0, ilocal, methloc, .false.,
     &        natoms, nae, nbe, nbsuse, nbasis, v(jend), XX, v(jatmap),
     &        v(jend+nbasis*nbasis*nri), lenv - jend)

         dims(1) = nbasis
         dims(2) = nbsuse
         if (nri .eq. 1) then
            call outoctfmt (iout, filenm, 2, 1, 'LCA', 2, dims, v(jend))
         elseif (nri .eq. 2) then
            call outoctfmt (iout, filenm, 2, 2, 'LCA', 2, dims, v(jend))
         endif

         if (nspblk .eq. 2) then
            call aclear (nri*nbasis*nbasis, v(jend))
            call fileio (2, -irwcb, nri*nbasis*nbsuse, v(jend), 0)

c$$$            call trspmo(.true.,.true.,nri,nbasis,nbasis,nbsuse,
c$$$     &           v(jend),v(jend + nri*nbasis*nbasis))
c$$$            call adjmo2(iout,iprint,.true.,I5DB1,I7FB1,nbasis,
c$$$     &           nbsuse*nri,NShell,ShellT,ShellC,v(jend))
c$$$            call trspmo(.false.,.true.,nri,nbasis,nbasis,
c$$$     &           nbsuse,v(jend),v(jend + nri*nbasis*nbasis))

            call localize(iout, iprint, 0, ilocal, methloc, .false.,
     &        natoms, nae, nbe, nbsuse, nbasis, v(jend), XX, v(jatmap),
     &        v(jend+nbasis*nbasis*nri), lenv - jend)

            dims(1) = nbasis
            dims(2) = nbsuse
            if (nri .eq. 1) then
               call outoctfmt (iout, filenm, 2, 1, 'LCB', 2, dims,
     &          v(jend))
            elseif (nri .eq. 2) then
               call outoctfmt (iout, filenm, 2, 2, 'LCB', 2, dims,
     &          v(jend))
            endif
         endif
      endif

C---------------------------------------------------------------------------
C     Output T2 amplitudes
C---------------------------------------------------------------------------

      if (outamp) then
C     load /orb/ block
         lenorb = intowp(28)
         call fileio(2,-irworb,lenorb,ispect,0)
C     Do some ugly additional size calculation
         nvava = nva*nva
         nvbvb = nvb*nvb
         novabt = noa*nob*nva*nvb
         novaat = noa*(noa-1)/2*nva*(nva-1)/2
         novbbt = nob*(nob-1)/2*nvb*(nvb-1)/2
         lenscr = max(noa*noa*nva*nva, nob*nob*nvb*nvb)
         jscr   = jend
         jend   = jscr + lenscr

         lenx = itqry(irwamp)
         write (iout, '(A,I20)') ' size of amp file = ', lenx
         lenx = 0
C     Here something strange happens, as my code coincides
C     with reading T1 in the (nv,no) form.
C     Apparently, Tom's code supports my point, and the actual order is (nv,no) (C-type)

C     dump T1 alpha part
         dims(1) = nva
         dims(2) = noa
         call fileio(2,-irwamp,noava,v(jend),lenx)
         lenx = lenx + noava
         call outoctfmt (iout, filenm, 2, 1, 't1a', 2, dims, v(jend))

         if(iopcl.ne.0) then
C     dump T1 beta part
            dims(1) = nvb
            dims(2) = nob
            call fileio(2,-irwamp,nobvb,v(jend),lenx)
            lenx = lenx + nobvb
            call outoctfmt (iout, filenm, 2, 1, 't1b', 2, dims,
     &           v(jend))

C     dump T2(alpha,alpha,alpha,alpha) - Warning: all T2 related code will
C     fail if any spatial symmetry was used for amplitude storage.
            if (novaat .gt. 0) then
c$$$                  dims(1) = noa
c$$$                  dims(2) = noa
c$$$                  dims(3) = nva
c$$$                  dims(4) = nva
               call fileio(2,-irwamp,novaat,v(jscr),lenx)
               lenx = lenx + novaat
C     unpack lower triangular
c$$$            call unpack_t2 (iout, noa, nva, v(jscr),
c$$$     &                 v(jend), lenv-jend)
c$$$                  call outoctfmt (iout, filenm, 2, 1, 't2aaaa', 4, dims,
c$$$     &                 v(jscr))
               dims(1) = novaat
               call outoctfmt (iout, filenm, 2, 1, 't2aaaa', 1, dims,
     &              v(jscr))

            endif
         endif
C     Here something strange happens, as my code coincides
C     with reading T2 in the (nv,nv,no,no) form, and commenst in Gaussian say it is (no,no,nv,nv).
C     Apparently, Tom's code supports my point, and the actual order is (nv,nv,no,no) (C-type)

C     dump T2(alpha,beta,alpha,beta)
         dims(1) = nvb
         dims(2) = nva
         dims(3) = nob
         dims(4) = noa
         call fileio(2,-irwamp,novabt,v(jend),lenx)
         lenx = lenx + novabt
         call outoctfmt (iout, filenm, 2, 1, 't2abab', 4, dims,
     &        v(jend))
         if(iopcl.ne.0) then
C     dump T2(beta,beta,beta,beta)
            if (novbbt .gt. 0) then
c$$$  dims(1) = nob
c$$$                  dims(2) = nob
c$$$                  dims(3) = nvb
c$$$                  dims(4) = nvb
               call fileio(2,-irwamp,novbbt,v(jscr),lenx)
               lenx = lenx + novbbt
c$$$                  call unpack_t2 (iout, noa, nva, v(jscr),
c$$$     &                 v(jend), lenv-jend)
c$$$                  call outoctfmt (iout, filenm, 2, 1, 't2bbbb', 4, dims,
c$$$     &                 v(jend))
               dims(1) = novbbt
               call outoctfmt (iout, filenm, 2, 1, 't2bbbb', 1, dims,
     &              v(jend))
            endif
         endif
      endif

C-----------------------------------------------------------------------------+
C     Print some summary                                                      |
C-----------------------------------------------------------------------------+

      if (iprint .gt. 0) then
         write (iout, '(A,I20)') ' basis size = ', nbasis
         write (iout, '(A,I20)') ' grid size = ', ngrid
         if (ri) write (iout, '(A,I20)') ' ndbf = ', ndbf
      endif

      call chainx(0)

      end
