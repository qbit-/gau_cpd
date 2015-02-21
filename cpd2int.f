

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

      real*8      cnvdef, acurcy, accdes
      real*8      vsgoal, defvsg
      integer     iprint, iphf, itype
      integer     ngrda, ngrdadf, ngrdamx
      integer     ngrdb, ngrdbdf, ngrdbmx
      integer     ngrdg, ngrdgdf, ngrdgmx
      integer     sval, sval1
      integer     maxcyc, cycdef, cycmxd
      integer     itst, idiis, mxscfc, ndidef, idisgo
      integer     ivshft, istrat, ifsel
      integer     igues1, idnmat, idnmt1, idnmt2, idnmt3, idnmts
      logical     icmplx, ispuhf, ispghf, isint, diis
      logical     evldmt, savedm, prtdm, dotwom, useao

      parameter   ( ngrdadf = 9, ngrdamx = 99 )
      parameter   ( ngrdbdf = 8, ngrdbmx = 99 )
      parameter   ( ngrdgdf = 9, ngrdgmx = 99 )
      parameter   ( cycdef = 512, cycmxd = 8192 )
      parameter   ( ndidef = 20 )
      parameter   ( defvsg = 2.0d0 )


C     memory allocation variables

      integer     jstrt, jend, lenv, mdv
      integer     nbas6d, ntt, ntt2, ntto, ntt2o, nbsq, nosq
      integer     szp, szf, szvec, szval, sznos
      integer     jp, jf, jval


C     PHF RWF files

      integer     irwp, irwf, irwvec, irwval, irwnos
      integer     irwciv, irw1dm, irw2dm

      parameter   ( irwp = 801, irwf = 802 )
      parameter   ( irwvec = 803, irwval = 804, irwnos = 805 )
      parameter   ( irwciv = 806, irw1dm = 807, irw2dm = 808 )


C     Gaussian RWF files
C       - irwgen - general array (see Link 1 for details)
C       - irwibf - file with IBF array
C       - irweig - MO orbital energies
C       - irwc?  - MO (x-spin) coefficients
C       - irwp?  - density matrix (x-spin)
C       - irws   - overlap matrix
C       - irwx   - transformation matrix (=S^(-1/2))
C       - irwh   - core Hamiltonian matrix

      integer     irwgen, irwibf, irweig
      integer     irwca, irwcb, irwpa, irwpb, irws, irwx, irwh

      parameter   ( irwgen = 501, irwibf = 508, irweig = 522 )
      parameter   ( irwca = 524, irwcb = 526, irwpa = 528, irwpb = 530 )
      parameter   ( irws = 514, irwx = 685, irwh = 515 )

C       - neq, neqshl, neqsh2 - files useful for symmetry

      integer     neq, neqshl, neqsh2

      parameter   ( neq = 580, neqshl = 565, neqsh2 = 726 )


C     symmetry related variables

      integer     jsym2e
      integer     nop1, nop2, nopuse
      integer     lenneq, lennes, lenne2
      integer     ineq, ineqsh, ineqs2


C     other variables

      real*8      enr, convg
      real*8      hsp, hph, energy, s2
      integer     nel, iopcl, junk
      integer     irws1, irwx1, irwh1
      logical     inobas


C     functions called

      real*8      scfacc, decacc
      integer     intowp, itqry


C     format statements

 1001 format (' PHF SCF was invoked with the following options:')
 1010 format ('   method selection:')
 1011 format (' iprint = ', I2)
 1012 format (' iphf   =  1, [ RHF + complex projection ]')
 1013 format (' iphf   =  2, [ UHF + complex projection ]')
 1014 format (' iphf   =  3, [ GHF + complex projection ]')
 1015 format (' iphf   =  4, [ UHF + spin projection ]')
 1016 format (' iphf   =  5, [ UHF + spin + complex projection ]')
 1017 format (' iphf   =  6, [ GHF + spin projection ]')
 1018 format (' iphf   =  7, [ GHF + spin + complex projection ]')
 1020 format ('   general options:')
 1021 format (' maxcyc = ', I6)
 1022 format (' maxcyc = ', I6, ', S      = ', I6)
 1023 format (' maxcyc = ', I6, ', S      = ', I4, '/2')
 1024 format (' ngrda  = ', I6, ', ngrdb  = ', I6, ', ngrdg  = ', I6)
 1030 format ('   accuracy requested:')
 1031 format (' conv criterion (delta energy) = ', 1P, E14.4)
 1032 format (' conv criterion (RMSDP, RMSDK) = ', 1P, E14.4)
 1033 format (' conv criterion (MaxP, MaxK)   = ', 1P, E14.4)
 1034 format (' requested integral accuracy   = ', 1P, E14.4)
 1040 format ('   use of symmetry:')
 1041 format (' two-electron integral symmetry not used')
 1042 format (' two-electron integrals replicated using symmetry')
 1050 format ('   strategy options:')
 1051 format (' level shifting is turned off')
 1052 format (' constant level shifting of ', F7.3, ' hartree')
 1053 format (' dynamic level shifting to a gap of ', F7.3, ' hartree')
 1054 format (' use of DIIS extrapolation is turned on')
 1055 format (' mxscfc = ', I6, ', idisgo = ', I6)
 1056 format (' istrat = ', I6, ', continue as usual if energy rises')
 1057 format (' istrat = ', I6, ', reduce DIIS space if energy rises')
 1058 format (' istrat = ', I6, ', increase level shft if energy rises')
 1059 format (' istrat = ', I6, ', apply DIIS every mxscfc cycles')
 1061 format (' regular PHF effective Fock matrix used')
 1062 format (' modified PHF effective Fock matrix used')
 1070 format ('   initial guess options:')
 1071 format (' igues1 =  1, read GDV guess')
 1072 format (' igues1 =  2, read GDV guess + complex phase to HOMO')

 1100 format (1X, '***', ' SCF Done after ', I4, ' cycles ***')
 1101 format (5X, '   E = ', F20.12, 3X, ' conv = ', 1P, E20.4)
 1102 format (5X, ' Hsp = ', F20.12, 3X, '  Hph = ', F20.12)
 1103 format (3X, '<S**2> = ', F20.12)
 1104 format (1X, '**********************************')

 1900 format (' PHF SCF - preparing initial guess')
 1901 format (' PHF SCF - starting SCF procedure')
 1902 format (' PHF SCF - quitting the program')



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


C     Set up prism and FMM control flags.

C       The calculation will crash if a semi-empirical Hamiltonian was
C       requested. FMM is set to default, but we will also ignore it.

      ifmm = 0
      iseall = 0

      call setpfl (iout, iprint, ifmm, ipflag, allowp, fmm, fmflag,
     $     fmflg1, nfxflg, ihmeth, omega, iseall, lseall, jstrt, v,
     $     ngot)

      if ( jstrt .lt. 1 ) jstrt = 1

      jend = jstrt
      lenv = ngot - jstrt + 1


C     Retrieve nuclear repulsion energy (enr) from general array.

      call fileio (2, -irwgen, lengen, dgen, 0)

      enr = dgen(41)


C     Read iopcl from ILSW file.

C       iopcl = 0,  real RHF
C             = 1,  real UHF
C             = 2,  complex RHF
C             = 3,  complex UHF
C             = 6,  complex GHF (there is no real GHF)

      call ilsw (2, 1, iopcl)

      if ( iopcl .gt. 3 .and. iopcl .ne. 6 ) then
        call gauerr ('Incorrect iopcl in phfdrv.')
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Options for PHF calculation  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      write (iout, *)
      write (iout, 1001)
      write (iout, *)


C     Set printing level.

      iprint = iop(6)

      if ( iprint .lt. 0 .or. iprint .gt. 10 ) iprint = 0



C     Type of PHF calculation
C     =======================

C       iphf - type of PHF calculation

C       iphf = 0,  [defaults to 1, 2, or 3, depending on Gaussian guess]
C            = 1,  closed-shell HF (RHF) + complex conjugation restoration
C            = 2,  open-shell HF   (UHF) + complex conjugation restoration
C            = 3,  general HF      (GHF) + complex conjugation restoration
C            = 4,  open-shell HF   (UHF) + spin projection
C            = 5,  open-shell HF   (UHF) + spin + complex projection
C            = 6,  general HF      (GHF) + spin projection
C            = 7,  general HF      (GHF) + spin + complex projection

      iphf = iop(7)

      if ( iphf .eq. 0 ) then
        if ( iopcl .eq. 0 .or. iopcl .eq. 2 ) then
          iphf = 1
        elseif ( iopcl .eq. 1 .or. iopcl .eq. 3 ) then
          iphf = 2
        elseif ( iopcl .eq. 6 ) then
          iphf = 3
        endif
      endif

      if ( iphf .le. 0 .or. iphf .gt. 7 ) then
        call gauerr ('Incorrect iphf in phfdrv.')
      endif


C     Decipher iphf.

C       icmplx - whether complex conjugation projection is on
C       ispuhf - whether UHF-based spin projection is used
C       ispghf - whether GHF-based spin projection is used
C       itype  - decides which spin blocks of transition density
C                matrices are stored

C       itype  = 1,  closed shell [alpha-alpha block]
C              = 2,  open shell [alpha-alpha, beta-beta blocks]
C              = 3,  general [all spin blocks are active]

      icmplx = iphf .le. 3 .or. iphf .eq. 5 .or. iphf .eq. 7
      ispuhf = iphf .eq. 4 .or. iphf .eq. 5
      ispghf = iphf .eq. 6 .or. iphf .eq. 7

      itype = 1

      if ( iphf .eq. 2 ) itype = 2
      if ( iphf .ge. 3 ) itype = 3


      write (iout, 1010)
      write (iout, 1011) iprint

      if ( iphf .eq. 1 ) write (iout, 1012)
      if ( iphf .eq. 2 ) write (iout, 1013)
      if ( iphf .eq. 3 ) write (iout, 1014)
      if ( iphf .eq. 4 ) write (iout, 1015)
      if ( iphf .eq. 5 ) write (iout, 1016)
      if ( iphf .eq. 6 ) write (iout, 1017)
      if ( iphf .eq. 7 ) write (iout, 1018)

      write (iout, *)



C     Integration grid
C     ================

C       ngrda - number of grid points to use for alpha integration
C       ngrdb - number of grid points to use for beta integration
C       ngrdg - number of grid points to use for gamma integration

      ngrda = iop(8)
      ngrdb = iop(9)
      ngrdg = iop(10)

      if ( .not. ( ispghf ) ) then
        ngrda = 1
        ngrdg = 1
      endif

      if ( .not. ( ispuhf .or. ispghf ) ) ngrdb = 1

      if ( ngrda .le. 0 ) ngrda = ngrdadf
      if ( ngrdb .le. 0 ) ngrdb = ngrdbdf
      if ( ngrdg .le. 0 ) ngrdg = ngrdgdf

      ngrda = min (ngrda, ngrdamx)
      ngrdb = min (ngrdb, ngrdbmx)
      ngrdg = min (ngrdg, ngrdgmx)



C     Spin projection
C     ===============

C       sval  - quantum number s on which to project
C       isint - .true. for integer s, .false. for half-integer s

C       sval is determined by iop(11), which selects the value of s
C       desired. Positive numbers should be used for integer s, negative
C       numbers for half-integer s. Note that s >= Sz = (Na - Nb) / 2,
C       where Na, Nb are the number of alpha and beta electrons,
C       respectively.

C       If iop (11) = 0, then sval = Sz = (Na - Nb) / 2.


      if ( abs (iop(11)) .gt. 0 ) sval1 = iop(11)

      if ( iop(11) .eq. 0 .or. .not. ( ispuhf .or. ispghf ) ) then
        if ( mod (nae-nbe, 2) .eq. 0 ) then
          sval = (nae-nbe)/2
          isint = .true.
        else
          sval = nae-nbe
          isint = .false.
        endif

      elseif ( mod (nae-nbe, 2) .eq. 0 ) then
        if ( sval1 .ge. 0 .and. sval1 .ge. (nae-nbe)/2 ) then
          sval = sval1
          isint = .true.
        else
          call gauerr ('Inconsistent multiplicity with iop (11)')
        endif

      else
        if ( sval1 .lt. 0 .and. -sval1 .ge. nae-nbe ) then
          sval = -sval1
          isint = .false.
        else
          call gauerr ('Inconsistent multiplicity with iop (11)')
        endif
      endif



C     Number of SCF iterations
C     ========================

C       maxcyc - maximum number of SCF cycles

      maxcyc = iop(12)

      if ( maxcyc .le. 0 ) maxcyc = cycdef

      maxcyc = min (maxcyc, cycmxd)

      write (iout, 1020)

      if ( ispuhf .or. ispghf ) then
        if ( isint ) write (iout, 1022) maxcyc, sval
        if ( .not. isint ) write (iout, 1023) maxcyc, sval

        write (iout, 1024) ngrda, ngrdb, ngrdg
      else
        write (iout, 1021) maxcyc
      endif

      write (iout, *)



C     Convergence criterion
C     =====================

C       itst - determines which variable should be used to test for
C              convergence

C       itst = 0,  [ same as 2 ]
C            = 1,  check for convergence in energy
C            = 2,  check for convergence in RMSDP
C            = 3,  check for convergence in MaxP

C       acurcy - requested convergence criterion

C       acurcy = 10^(-N),  where N is set by iop(14)  [ default = 8 ]


      itst = iop(13)

      if ( itst .gt. 3 .or. itst .lt. 1 ) itst = 2


      cnvdef = scfacc (0)

      if ( iop(14) .ge. 0 ) then
        acurcy = scfacc (iop(14))
      else
        acurcy = scfacc (0)
      endif

      accdes = decacc (0, 0.0d0, .false.)

      if ( acurcy .lt. cnvdef ) then
        accdes = accdes * acurcy / cnvdef
      endif

      write (iout, 1030)

      if ( itst .eq. 1 ) then
        write (iout, 1031) acurcy
      elseif ( itst .eq. 2 ) then
        write (iout, 1032) acurcy
      elseif ( itst .eq. 3 ) then
        write (iout, 1033) acurcy
      endif

      write (iout, 1034) accdes
      write (iout, *)



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

      if ( jsym2e .eq. 0 ) write (iout, 1041)
      if ( jsym2e .eq. 2 ) write (iout, 1042)

      write (iout, *)



C     Strategy options
C     ================

C       ivshft - controls level shifting

C       ivshft = -N,  dynamic level shifting to achieve a gap
C                     of -0.001*N
C              = -2,  dynamic level shifting to a default goal
C                     (same as -2000)
C              = -1,  no level shifting
C              =  0,  [same as -1]
C              =  N,  shift by 0.001*N

      ivshft = iop(15)

      if ( ivshft .eq. 0 ) ivshft = -1

      if ( ivshft .eq. -2 ) then
        vsgoal = defvsg
      else
        vsgoal = abs (dble (ivshft) / dble (1000))
      endif

      write (iout, 1050)

      if ( ivshft .eq. -1 ) then
        write (iout, 1051)
      elseif ( ivshft .gt. 0 ) then
        write (iout, 1052) vsgoal
      elseif ( ivshft .lt. -1 ) then
        write (iout, 1053) vsgoal
      endif


C       diis - whether to use DIIS or not

      idiis = iop(16)

      if ( idiis .eq. 1 ) then
        diis = .true.
      else
        diis = .false.
      endif

      if ( diis ) write (iout, 1054)


C       mxscfc - number of DIIS vectors to use

      mxscfc = iop(17)

      if ( diis .and. mxscfc .eq. 0 ) mxscfc = min (ndidef, maxcyc)
      
      diis = diis .and. mxscfc .gt. 1


C       idisgo - cycle at which DIIS should start

      idisgo = iop(18)

      if ( diis .and. idisgo .eq. 0 ) idisgo = 1

      if ( diis ) write (iout, 1055) mxscfc, idisgo


C       istrat - strategy options

C     with DIIS

C       istrat = 0,  [same as 1]
C              = 1,  continue as usual if energy rises
C              = 2,  reduce DIIS space if energy rises from previous
C                    cycle
C              = 3,  apply DIIS every mxscfc cycles

C     without DIIS, with level shifting

C       istrat = 0,  [same as 1]
C              = 1,  continue as usual if energy rises
C              = 2,  increase constant level shift or target gap by
C                    a factor of 1.5 if energy rises from previous cycle

      istrat = iop(19)

      if ( istrat .eq. 0 ) istrat = 1

      if ( istrat .eq. 1 ) then
        write (iout, 1056) istrat

      elseif ( istrat .eq. 2 ) then
        if ( diis ) then
          write (iout, 1057) istrat
        else
          write (iout, 1058) istrat
        endif

      elseif ( istrat .eq. 3 ) then
        write (iout, 1059) istrat

        if ( .not. diis ) then
          call gauerr ('Cannot use istrat = 3 without DIIS.')
        endif
      endif



C     Control of Fock matrix structure
C     ================================

C       ifsel - controls the structure of Fock matrix to be diagonalized
C               at every iteration

C       ifsel = 0,  [defaults to 1]
C             = 1,  use PHF effective Fock matrix
C             = 2,  use PHF effective Fock matrix for ov and vo blocks;
C                   use regular HF Fock matrix for oo and vv blocks

      ifsel = iop(22)

      if ( ifsel .ne. 1 ) ifsel = 2

      if ( ifsel .eq. 1 ) then
        write (iout, 1061)
      elseif ( ifsel .eq. 2 ) then
        write (iout, 1062)
      endif

      write (iout, *)



C     Initial guess options
C     =====================

C       igues1 - determines how to prepare initial guess

C       igues1 = 0,  [ same as 1 ]
C              = 1,  read Gaussian initial guess
C              = 2,  read Gaussian initial guess +
C                    apply complex phase factor to some of the elements
C                    of HOMO orbital (breaks cc symmetry)

      igues1 = iop(20)

      if ( igues1 .eq. 0 ) igues1 = 1

      write (iout, 1070)

      if ( igues1 .eq. 1 ) then
        write (iout, 1071)
      elseif ( igues1 .eq. 2 ) then
        write (iout, 1072)
      endif

      write (iout, *)



C     Density matrix options
C     ======================

C       idnmat - controls the evaluation of density matrices associated
C                with the PHF state

C       idnmat =   0,  [ same as 1 ]
C              =   1,  do not evaluate density matrices
C              =   2,  build 1PDM just to obtain natural orbital occupations
C              =   3,  build density matrices and save in RWF file
C              =   4,  build density matrices and output to external files
C                      ( 1pdm.dat, 2pdm.dat )

C              =  00,  [ same as 10 ]
C              =  10,  build only the 1PDM
C              =  20,  build both the 1PDM and the 2PDM

C              = 000,  [ same as 100 ]
C              = 100,  build density matrices in orthonormal AO basis
C              = 200,  build density matrices in regular AO basis

      idnmat = iop(21)

      idnmt1 = mod (idnmat, 10)
      idnmt2 = mod (idnmat, 100) / 10
      idnmt3 = mod (idnmat, 1000) / 100

      if ( idnmt1 .eq. 0 .or. idnmt1 .gt. 4 ) idnmt1 = 1
      if ( idnmt2 .eq. 0 .or. idnmt2 .gt. 2 ) idnmt2 = 1
      if ( idnmt3 .eq. 0 .or. idnmt3 .gt. 2 ) idnmt3 = 1

      evldmt = .false.
      savedm = .false.
      prtdm  = .false.

      if ( idnmt1 .eq. 1 ) then
      elseif ( idnmt1 .eq. 2 ) then
        evldmt = .true.
      elseif ( idnmt1 .eq. 3 ) then
        evldmt = .true.
        savedm = .true.
      elseif ( idnmt1 .eq. 4 ) then
        evldmt = .true.
        savedm = .true.
        prtdm  = .true.
      endif

      dotwom = .false.

      if ( idnmt2 .eq. 1 ) then
      elseif ( idnmt2 .eq. 2 ) then
        if ( savedm ) dotwom = .true.
      endif

      useao = .false.

      if ( idnmt3 .eq. 1 ) then
      elseif ( idnmt3 .eq. 2 ) then
        if ( savedm ) useao = .true.
      endif

      
C       idnmts - controls which density matrix is stored in the
C                appropriate Gaussian RWF file

C       idnmts = 0, [defaults to 1]
C              = 1, store the deformed HF state density matrix
C              = 2, store the PHF density matrix
C                   ( requires mod (iop(21), 10) >= 3 and
C                              mod (iop(21), 1000) / 100 = 1 )**

C         ** in words, the PHF density matrix has to be stored
C            in its own RWF file AND it has to be constructed in
C            the ORTHONORMAL AO basis

      idnmts = 1

      if ( iop(23) .eq. 2 .and. savedm .and. .not. useao ) then
        idnmts = 2
      endif



C     %%%%%%%%%%%%%%%%%%%%%%%
C     %  Memory allocation  %
C     %%%%%%%%%%%%%%%%%%%%%%%


C     Useful quantities.

      nel = nae + nbe

      call getnb6 (nbas6d)

      ntt   = nbasis * (nbasis + 1) / 2
      ntt2  = nbasis * (2*nbasis + 1)

      ntto  = nbsuse * (nbsuse + 1) / 2
      ntt2o = nbsuse * (2*nbsuse + 1)

      nbsq = nbasis * nbasis
      nosq = nbsuse * nbsuse


C     Define some quantities.
C       - size of density matrix (szp)
C       - size of Fock matrix (szf)
C       - size of eigenvectors of Fock (szvec)
C       - size of eigenvalues of Fock (sval)
C       - size of transformation matrix: orthonormal AO => NO (sznos)

      if ( itype .eq. 1 ) then
        szp = 2*ntto
        szf = 2*ntto
        szvec = 2*nosq
        szval = nbsuse
        sznos = 2*nosq

      elseif ( itype .eq. 2 ) then
        szp = 4*ntto
        szf = 4*ntto
        szvec = 4*nosq
        szval = 2*nbsuse
        sznos = 4*nosq

      elseif ( itype .eq. 3 .and. ispuhf ) then
        szp = 4*ntto
        szf = 4*ntto
        szvec = 4*nosq
        szval = 2*nbsuse
        sznos = 8*nosq

      elseif ( itype .eq. 3 ) then
        szp = 2*ntt2o
        szf = 2*ntt2o
        szvec = 8*nosq
        szval = 2*nbsuse
        sznos = 8*nosq
      endif
        

C     Allocate space for:
C       - HF density matrix (jp)
C       - effective Fock matrix (jf)

      jp   = jend
      jf   = jp   + szp
      jend = jf   + szf

      mdv = lenv - jend + 1

      call tstcor (jend-1, lenv, 'phfdrv')


C     Initialize PHF RWF files.

      call conddf (irwp, szp)
      call conddf (irwf, szf)
      call conddf (irwvec, szvec)
      call conddf (irwval, szval)
      call conddf (irwnos, sznos)



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Retrieve symmetry information  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C     Get number of symmetry operations.
C       - nop1,    order of the concise abelian subgroup
C       - nop2,    order of the full abelian subgroup
C       - nopuse,  actual oder to be used

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

      call tstcor (jend-1, lenv, 'phfdrv')


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


      call chainX (0)


      return
      end


