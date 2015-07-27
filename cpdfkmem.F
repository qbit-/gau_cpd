      subroutine cpdfkmem (iout, iprint, initf, nbasis, ngrid,
     & x, z, nmat, nmats, nmatt, p, f, v, lenv)

      implicit none
      
C-----------------------------------------------------------------------------+
C     This routine contracts densities with canonically decomposed integrals
C     in a way similar to FoFCou operation. It is intended to be called from cpdfock,
C     and works with one spin case per run.
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
      integer     iout, iprint
      integer     nbasis, ngrid
      integer     nmat, nmats, nmatt
      real*8      x(*), z(*)
      real*8      p(*), f(*)
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

C     do full contraction with the singlet symmetric part
      
      do 100 ii = 0, nmats - 1
         call ctrcou (iout, iprint, initf, 0, nbasis, ngrid, x,
     &        z, p(1+ii*ntt), f(1+ii*ntt), v)
         call ctrexc (iout, iprint, .false., 0, nbasis, ngrid, x,
     &        z, p(1+ii*ntt), f(1+ii*ntt), v)
 100  continue
      
C     do exchange part for triplets (symmetric)
      
      ioff = ntt*nmats + 1
      do 200 ii = 0, nmatt - 1
         call ctrexc (iout, iprint, initf, 0, nbasis, ngrid, x,
     &        z, p(ioff+ii*ntt), f(ioff+ii*ntt), v)
 200  continue

C     do exchange for antisymmetric part

      ioff  = ntt*(nmats + nmatt) + 1
      nmata = nmat - nmats - nmatt
      do 300 ii = 0, nmata - 1
         call ctrexc (iout, iprint, initf, 1, nbasis, ngrid, x,
     &        z, p(ioff+ii*ntt), f(ioff+ii*ntt), v)
 300  continue

      end
