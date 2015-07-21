      subroutine frm2eri (iout, iprint, jsym2e, neq, neqshl, 
     &     neqsh2, ir12, keepov, keepints, itol, accdes, jovoff,
     &     next, irwria, irwrib, v, lenv)
      implicit none
C     
C     This is a wrapper routine to form RI overlap and integrals and place them 
C     into memory or in rwfs. Symmetry is not currently used. Daugther 
C     routines are based on rifrma and rifrmb from l906.
C     FixMe: 6d basis is not properly accounted for 
     
C     Arguments:
C     [in]  jsym2e             - controls symmetry.           
C     [in]  neq, neqshl, neqs2 - symmetry related RWFs from Gaussian
C     [in]  ir12               - metric to use in RI (1 - 1/r, 2 - overlap)
C     [in]  keepov             - if we need to keep the overlap inverse in memory
C     [in]  keepints           - if we need to keep integrals in memory after return
C     [in]  itol               - pseudoinversion tolerance
C     [in]  accdes             - desired accuracy from PRISM
C
C     [out] jovoff            - offset to RI ints. If overlap hasn't been kept in memory this
C                                is 0 
C     [out] next              - offset to the next free location in the memory. If integrals neither overlap
C                               are kept this is 0
C     [out] irwria            - RWF with DBF overlap
C     [out] irwrib            - RWF with RI integrals

#include "commonmol.inc"      
      
C     arguments

      integer    lenv
      integer    iout, iprint
      real*8     v(*)

      integer    jsym2e, neq, neqshl, neqsh2, ir12
      integer    irwria, irwrib
      integer    jovoff, next
      logical    keepov, keepints

      real*8     accdes

C     prism variables

      integer     ifmm, iopcl, ipflag, fmflag, fmflg1, nfxflg, ihmeth
      integer     iseall
      integer     lseall, momega, nomega
      logical     allowp(50), fmm
      parameter   ( momega = DEFMOMEGA, nomega = 6 )
      real*8      omega(momega,nomega)

      integer     jpflsta

C     other

      integer     jend, jria, jrints, junk
      integer     nbas6d, ndbf, ndbf6d   

      integer     itol
C     
      jend    = 1
      jovoff  = 0
      next    = 0

      call getnb6(nbas6d)
      
C     set up PRISM control flags

      ifmm = 0
      iseall = -1

      jpflsta = jend

      call setpfl (iout, iprint, ifmm, ipflag, allowp, fmm, fmflag,
     &     fmflg1, nfxflg, ihmeth, omega, iseall,lseall,jpflsta,v(jend),
     &     lenv-jend+1)

C     form the overlap. The sizes of auxiliarry basis are also returned here

      call frmria (iout, iprint, jovoff, irwria, ndbf, ndbf6d, keepov,
     &     ir12, itol, ipflag, allowp, natoms, ian, atmchg, c, iattyp,
     &     v(jend), lenv - jend)

      if ( keepov ) jend = jend + jovoff
      
      call frmrib (iout, iprint, next, irwrib,
     &     keepints, ipflag, allowp, accdes, ndbf, ndbf6d,
     &     nbas6d, natoms, ian, atmchg, c, iattyp, v(jend), lenv - jend)

      if ( keepints ) then
         jend = jend + next
         next = next + jovoff
      endif

      end
