      subroutine cpdfkmem (iout, iprint, iopcl, initf, nbasis, ngrid,
     & x, z, nmat, nmats, nmatt, pa, pb, fa, fb, v, lenv)

      implicit none
      
C-----------------------------------------------------------------------------+
C     This routine contracts densities with canonically decomposed integrals
C     in a way similar to FoFCou operation.
C     The integrals are Mulliken order, stored as:
C
C     (mn|rs) = X_{I,m} X_{I,n} Z_{I,J} X_{J,r} X_{J,s}
C     (repeated indices are summed)
C      
C     Densities are supplied in lower triangular 
C     format in vectors pa and pb, fock matrices are returned lower triangular in 
C     fa and fb.
C     Symmetric matrices are contracted with Coulomb and exchange:
C     ((mn|rs) + (ms|rn)) P_{s,r}
C     while triplet and
C     antisymmetric (singlet and triplet) are contracted with exchange
C     (ms|rn) P_{s,r}
C      
C     IMPORTANT: Z matrix is assumed symmetric
C     (but the code is not optimized for that), otherwise we would need to modify the
C     code to contract singlet antisymmetric densities with Coulomb part
C     
C     The formulas for contraction can be found in ctrcou and ctrexc
C      
C     Arguments:
C     [in]  iopcl   - type of the contraction (open/closed shell)
C     [in]  initf   - if we need to initialize the Fock matrix or to update current value
C     [in]  nbasis  - basis size
C     [in]  ngrid   - quadrature size
C     [in]  x       - nbasis*ngrid quadrature
C     [in]  z       - ngrid*ngrid Z matrix
C     [in]  nmat    - total number of lower triangular matrices in PA and PB 
C     [in]  nmats   - number of singlets stored in 1..nmats first lower triangular matrices in P
C                     (are contracted with both Coulomb and exchange)
C     [in]  nmatt   - number of triplets stored in nmats+1..nmatt positions in P
C                     (are contracted with exchange only)
C     [in]  p       - vector of nmat lower triangular matrices density matrices
C     [out] f       - vector of nmat lower triangular Fock matrices
C     
C-----------------------------------------------------------------------------+
C     arguments

      logical     initf
      integer     iopcl
      integer     iout, iprint
      integer     nbasis, ngrid
      integer     nmat, nmats, nmatt
      real*8      x(*), z(*)
      real*8      pa(*), pb(*), fa(*), fb(*)
      real*8      v(*)
      integer     lenv

C     pointers & sizes
      
      integer     lenx, lenz
      integer     lenscrco, lenscrex
      
C     other variables

      integer     ntt, nbsq, ngsq
      integer     ii
      integer     ioff, nmata
      
C     define sizes

      ntt  = nbasis*(nbasis+1)/2
      nbsq = nbasis*nbasis
      ngsq = ngrid*ngrid

      lenx = nbasis*ngrid
      lenz = ngsq
      
C     Prepare for parallel run (NYI)
      
      lenscrco = nbsq + lenx + ngrid
      lenscrex = lenx + ngsq
      
      call tstcor (max(lenscrco,lenscrex), lenv, 'cpdfkmem')

C     Do Coulomb part
      ioff = 1
      do 100 ii = 0, nmats - 1
         call ctrcou (iout, iprint, initf, 0, nbasis, ngrid, x,
     &        z, pa(1+ii*ntt), fa(1+ii*ntt), v)
 100  continue

C     multiply Coulomb by 2 if we are doing closed shell or
C     add beta contribution for open shell
      
      if (iopcl .eq. 0) then
         do 110 ii = 0, nmats - 1
            call ascale(ntt, 2.0d0, fa(1+ii*ntt),fa(1+ii*ntt))
 110     continue
      else
         do 120 ii = 0, nmats - 1
            call ctrcou (iout, iprint, .false., 0, nbasis, ngrid, x,
     &           z, pb(1+ii*ntt), fa(1+ii*ntt), v)
 120     continue
         if (initf) then
            do 130 ii = 0, nmats - 1
               call amove (ntt,fa(1+ii*ntt),fb(1+ii*ntt))
 130        continue
         else
            do 140 ii = 0, nmats - 1
            call ctrcou (iout, iprint, .false., 0, nbasis, ngrid, x,
     &           z, pb(1+ii*ntt), fb(1+ii*ntt), v)
            call ctrcou (iout, iprint, .false., 0, nbasis, ngrid, x,
     &           z, pa(1+ii*ntt), fb(1+ii*ntt), v)
 140        continue
         endif
      endif

C     Add exchange to PA and PB

      if (iopcl .eq. 0) then
         call cpdexmem (iout, iprint, initf, nbasis, ngrid,
     &        x, z, nmat, nmats, nmatt, pa, fa, v, lenv)
      else
         call cpdexmem (iout, iprint, initf, nbasis, ngrid,
     &        x, z, nmat, nmats, nmatt, pa, fa, v, lenv)
         call cpdexmem (iout, iprint, initf, nbasis, ngrid,
     &        x, z, nmat, nmats, nmatt, pb, fb, v, lenv)
         
      endif
      
      end
