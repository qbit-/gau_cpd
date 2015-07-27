      subroutine cpdfock (iout, iprint, iopcl, initf, nbasis, ngrid,
     & irwfac, irwz, nmat, nmats, nmatt, pa, pb, fa, fb, v, lenv)

      implicit none
      
C-----------------------------------------------------------------------------+
C     This routine contracts densities with canonically decomposed integrals
C     in a way similar to FoFCou operation. The factors of the integrals are
C     coming from RWFs irwx and irwz. The integrals are Mulliken order, stored as:
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
C     [in]  initf   - if we need to initialize the Fock matrix or to update current value
C     [in]  iopcl   - reminiscent to Gaussian's iopcl, used to specify if
C                     beta matrices are contracted. 0 - use only PA, 1 - use both PA and PB 
C     [in]  nbasis  - basis size
C     [in]  ngrid   - quadrature size
C     [in]  irwfac  - RWF which expected to contain nbasis*ngrid quadrature
C     [in]  irwz    - RWF which expected to contain ngrid*ngrid Z matrix
C     [in]  nmat    - total number of lower triangular matrices in PA and PB 
C     [in]  nmats   - number of singlets stored in 1..nmats first lower triangular matrices in P
C                     (are contracted with both Coulomb and exchange)
C     [in]  nmatt   - number of triplets stored in nmats+1..nmatt positions in P
C                     (are contracted with exchange only)
C     [in]  pa, pb  - vectors of nmat alpha/beta lower triangular matrices
C     [out] fa, fb  - vectors of nmat alpha/beta Fock matrices
C     
C-----------------------------------------------------------------------------+
C     arguments

      logical     initf 
      integer     iout, iprint
      integer     iopcl, nbasis, ngrid
      integer     nmat, nmats, nmatt
      integer     irwfac, irwz
      real*8      pa(*), pb(*)
      real*8      fa(*), fb(*)
      real*8      v(*)
      integer     lenv

C     pointers
      
      integer     jend, jx, jz
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
      
      jend = 1
      jx   = jend
      jz   = jx + lenx
      jend = jz + lenz

      call tstcor (jend - 1, lenv, 'cpdfock')
      
C     Load matrices of CPD decomposed 2e integrals

      call fileio (2, -irwfac, lenx, v(jx), 0)
      call fileio (2, -irwz, lenz, v(jz), 0)

C     input checks
      
      if ((iopcl .ne. 0) .and. (iopcl .ne. 1)) 
     &     call gauerr ('cpdfock: unknown value of iopcl')      
      
      if ((nmats + nmatt) .gt. nmat) 
     &     call gauerr ('cpdfock: nmats + nmatt > nmat')      

      call cpdfkmem (iout, iprint, iopcl, initf, nbasis, ngrid,
     &        v(jx), v(jz), nmat, nmats, nmatt,
     &        pa, fa, v(jend), lenv - jend)

      if (iopcl .eq. 1) then 
         call cpdfkmem (iout, iprint, iopcl, initf, nbasis,
     &        ngrid,
     &        v(jx), v(jz), nmat, nmats, nmatt,
     &        pb, fb, v(jend), lenv - jend)
      endif

      end
