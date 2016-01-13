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
C |  iop(79) - requests to read external grid from grid.dat    |
C |                0 - build grid normally                     |
C |               >0 - read external grid                      |     
C |  iop(80) - what to output                                  |
C |                1 - dump quadrature X                       |
C |                2 - dump overlap and AO to MO transform     |
C |                4 - dump MOs and eigenvalues                |
C |               10 - dump DF integrals and DF overlap        |
C |               20 - dump full integrals                     |
C |         options can be added, example: 13 - DF integrals + |
C |                  quadrature + overlap                      |      
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
      logical     dopure
      
      integer     idump, idumpdef, idp
      logical     outquad, outaox, outmos, outint, outintdf

      parameter   (idumpdef = 21)
      
C     memory allocation variables

      integer     jend, jovoff, jcsr, lenx, next, lenv
      
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
      character*(*)   filenm
      parameter   ( filenm = 'dump_all.mat')

C     Needed for commonb2 to work with implicit none

      integer    nshelb, janb, maxtyb 

C     RWF files

      integer     irwfac, irwspo, irwz, irwe

      parameter   ( irwfac = 824, irwspo = 826 )
      parameter   ( irwe = 827, irwz = 828 )

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

C     Try to recover iradan from a previous CPD run

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
      outint = ( (idp / 2) .ge. 1 )
      idp = mod(idp,2)
      outintdf = (idp .eq. 1)
      
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

      
C     Clean the previous file if it existed
      open(unit=771, iostat=ios, file=filenm, status='old')
      if (ios .eq. 0) close(771, status='delete')

C-----------------------------------------------------------------------------+
C     Dump some scalars                                                       |
C-----------------------------------------------------------------------------+

      call fileio (2, -irwgen, lengen, dgen, 0)

C     SCF energy
      call outoctfmt (iout, filenm, 2, 1, 'e_scf', 1, 1, dgen(32))
C     CC/CI energy
      call outoctfmt (iout, filenm, 2, 1, 'e_cc', 1, 1, dgen(9))

      
C-----------------------------------------------------------------------------+
C     Compute the collocation matrix X(mu,G) and save it to the RWF irwfac    |
C-----------------------------------------------------------------------------+
      
      if (extgrid) then
         call frmextquad (iout, iprint, 'grid.dat', 0, 
     &        nbas6d, ngrid, irwfac, v(jend), lenv - jend)
      else
         call frmquad (iout, 2, iradan, iranwt, irangd, iopcl,
     &        .false., nbsuse, junk, ngrid, irwfac, v(jend),
     &        lenv - jend)
      endif
      
C     Dump quadrature
      
      if (outquad) then
         call tstcor(nbas6d*ngrid,lenv,'dump_all: quadrature')      
         call fileio (2, -irwfac, nbas6d*ngrid, v(jend), 0)
         
         dims(1) = nbas6d
         dims(2) = ngrid
         call outoctfmt (iout, filenm, 2, 1,'Xu', 2, dims,
     &        v(jend))

         call fileio (2, -irwfac, nbas6d*ngrid, v(jend), 0)

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

      if (outmos) then
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


