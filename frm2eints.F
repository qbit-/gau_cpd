      subroutine frm2eints (iout, iprint, jsym2e, neq, neqshl, 
     &     neqsh2, keep, accdes, next, irw2efull, v, lenv)
      
      implicit none

C     This is a simple routine to call PRISM and
C     store 2e integrals to an RWF (needs to be replaced
C     by a smarter code).
C     The integrals are stored in the full form without 8-fold simmetry
C     FixMe: 6d basis is not properly accounted for

C     Arguments:
C     [in]  jsym2e             - controls symmetry.           
C     [in]  neq, neqshl, neqs2 - symmetry related RWFs from Gaussian
C     [in]  keep               - if we need to keep integrals in memory after return
C     [in]  accdes             - desired accuracy from PRISM
C
C     [out] next               - offset to the next free location in the memory
C                                if integrals are kept there
C     [out] irw2efull          - RWF with integrals

#include "commonmol.inc"
#include "commonb.inc"
#include "repall.inc"

C     arguments

      integer    lenv
      integer    iout, iprint
      real*8     v(*)

      integer    jsym2e, neq, neqshl, neqsh2
      integer    next, irw2efull
      logical    keep

      real*8     accdes

C     symmetry related variables

      integer     nop1, nop2, nopuse
      integer     lenneq, lennes, lenne2
      integer     ineq, ineqsh, ineqs2

C     PBC variables

      integer     npdir, jpbcsta, jpbcend
      real*8      celvec(3,3)

C     prism variables

      integer     ifmm, iopcl, ipflag, fmflag, fmflg1, nfxflg, ihmeth
      integer     iseall
      integer     lseall, momega, nomega
      logical     allowp(50), fmm
      parameter   ( momega = DEFMOMEGA, nomega = 6 )
      real*8      omega(momega,nomega)

      integer     jpflsta
      integer     icntrl, i1cent, iraf, nbti
      integer     iprtf

C     other

      integer     jend, j2eint, jscr, junk
      integer     len2eint
      integer     ls, ntt, ntt2, nbsq, nbas6d   

      character*20 chap
      real*8      xx
      complex*16  cxx

C     functions called

      integer     intowp, itqry

C     constants

      real*8 two, half
      complex*16  cone
      parameter ( two = 2.0d0, half = 0.5d0 )
      parameter ( cone = (1.0d0, 0.0d0) )

C     debug
      
      integer     i, j, k, l
      logical     dbg

C     contraction

      

C     Memory allocation

      jend = 1

      ntt  = (nbasis*(nbasis + 1))/2
C      ntt2 = (ntt*(ntt+1))/2
      nbsq = nbasis*nbasis
      call getnb6(nbas6d)

C     allocate memory for integrals. For an easy setup for PRISM see l906

      len2eint = nbsq*nbsq
      j2eint = jend
      jscr   = j2eint + len2eint
      jend   = jscr + nbsq
      
      call tstcor (jend-1, lenv, 'frm2eints: integrals')

      if (keep) then
         next = jscr
      else
         next = 0
      endif

C     Initialize RWF for integrals

      call conddf (irw2efull, nbsq*nbsq)

C     Retrieve symmetry information for PRISM (symmetry is ignored for now)
      
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

      if ( jsym2e .eq. 1 ) then
        lenneq = intowp (nop1*nbasis)
        lennes = max (itqry (neqshl), 0)
        lenne2 = max (itqry (neqsh2), 0)
      endif

      ineq   = jend
      ineqsh = ineq   + lenneq
      ineqs2 = ineqsh + lennes
      jend   = ineqs2 + lenne2

      call tstcor (jend-1, lenv, 'frm2eints: symmetry')

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

C     set up PBC information (not needed, but leave here for future)
      
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
     $     fmflg1, nfxflg, ihmeth, omega, iseall,lseall,jpflsta,v(jend),
     $     lenv-jend+1)


C     Setup options for FoFcou
C     Request to return integrals in canonical order

      iopcl  = 0
      icntrl = 600 
      i1cent = 0

C     Supress excessive output
      iprtf = iprint - 3
      if (iprtf .lt. 0) iprtf = 0
      
C     Rafenetti control according to l906
      iraf   = 0

C     turn off symmetry
      jsym2e = 0
      nopuse = 1
      nop1   = 1
      
C     no idea what it is (set according l906)
      nbti   = 0
      
        call fofcou(iout,iprtf,ihmeth,iopcl,icntrl,iraf,nomega,
     &  omega,fmm,fmflag,fmflg1,nfxflg,ipflag,allowp,lseall,lseall,
     &  junk,junk,.false.,.false.,.false.,.false.,.false.,i1cent,accdes,
     &  cone,0,1,1,0,1,0,0,nbasis,
     &  jsym2e,nopuse,nop1,v(ineq),v(ineqsh),v(ineqs2),rotall,neqall,
     &  nbti,junk,xx,
     &  xx,xx,xx,xx,xx,xx,xx,xx,junk,xx,xx,junk,
     &  xx,xx,natoms,ian,atmchg,c,iattyp,iatfrg,xx,xx,v(jpbcsta),
     &  celvec,cxx,cone,v(j2eint),junk,junk,xx,xx,junk,
     &  junk,xx,v(jend),v(jend),lenv-jend+1)

C     Put integrals in a square order and also restore 
C     1/2 factors which PRISM put for I != J and K != L

        call square (v(j2eint), v(j2eint), nbsq, ntt, 0)

        do 100 ls = 1, ntt
           call square (v(j2eint+(ls-1)*nbsq),v(j2eint+(ls-1)*nbsq),
     &          nbasis,nbasis,0)
           call dscal (nbsq, two, v(j2eint+(ls-1)*nbsq), 1)
           call dscal (nbasis, half, v(j2eint+(ls-1)*nbsq), nbasis+1)
 100    continue

        call trspn2 (nbsq,v(j2eint))
        
        do 200 ls = 1, nbsq
           call square (v(j2eint+(ls-1)*nbsq),v(j2eint+(ls-1)*nbsq),
     &          nbasis,nbasis,0)
           call dscal (nbsq, two, v(j2eint+(ls-1)*nbsq), 1)
           call dscal (nbasis, half, v(j2eint+(ls-1)*nbsq), nbasis+1)
 200    continue

C     save integrals

      if (.not. keep) then
         call fileio (1, -irw2efull, len2eint, v(j2eint), 0) 
      endif

C     Debug printing

      dbg = .false.

      if (dbg) then

 55     format ('(', I4,',',I4,'|',I4,',',I4,') = ',F6.4)
        do i = 0, nbasis-1
           do j = 0, nbasis-1
              do k = 0, nbasis-1
                 do l = 0, nbasis-1
                 write (iout,55) i,j,k,l,v(j2eint+i*nbsq*nbasis+
     & j*nbsq + k*nbasis + l)
              enddo
           enddo
        enddo
      enddo

      endif

      end
